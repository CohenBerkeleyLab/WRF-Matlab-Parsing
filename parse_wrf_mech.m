function [ J, species, isfixed, photo_calls ] = parse_wrf_mech( mech_name )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('mech_name','var')
    mech_name = 'r2smh-simple';
end

species_file = fullfile('/Users/Josh/Documents/MATLAB/BEHR/WRF_Utils/Models/R2SMH',strcat(mech_name,'.spc'));
eqn_file = fullfile('/Users/Josh/Documents/MATLAB/BEHR/WRF_Utils/Models/R2SMH',strcat(mech_name,'.eqn'));

% The mechanism will be constructed as a vector of anonymous functions that
% accept a concentration vector and photolysis rate vector.
%
% An additional emissions vector will be added at each timestep as well.

[species, isfixed] = parse_species(species_file);
[J, photo_calls] = construct_mechanism(eqn_file, species);
J(isfixed) = {0};
end

function [species, isfixed] = parse_species(species_file)
% This function will return the list of species defined for the mechanism
% as well as a logical vector indicating if those species are to be fixed
% by some physical process rather than the chemistry.
%
% First read in the species names from the species file, as well as if it
% is meant to be solved by chemistry or fixed by physical processes.
species = cell(1,1000); % will cut down in the end
isfixed = false(1,1000);
fixed_i = 0;
i=1;
fid = fopen(species_file);
tline = fgetl(fid);
while ischar(tline)
    tline = strtrim(tline);
    if ismember('#DEFVAR',tline) % the following species are variable
        fixed_i = 0;
    elseif ismember('#DEFFIX',tline) % the following species are fixed
        fixed_i = 1;
    elseif ismember('#',tline) % unanticipated definition bloc
        fprintf('Unknown definition block: %s\n',tline);
    else
        % Some lines have curly braces. Since I don't know exactly what
        % they mean, and I'm only interested in what species are needed,
        % they can be removed.
        parsed_line = strrep(tline,'{','');
        parsed_line = strrep(parsed_line,'}','');
        parsed_line = strsplit(parsed_line,'='); % species comes before the = sign
        spec = strtrim(parsed_line{1});
        xx = strcmp(spec, species);
        if sum(xx) > 0 && isfixed(xx) ~= fixed_i
            error('species_definition:species_fixed_and_var','Specie %s is defined as both a fixed and variable specie',spec)
        elseif sum(xx) > 0
            fprintf('Redundant definition of %s, skipping\n',spec);
        else
            species{i} = spec;
            isfixed(i) = fixed_i;
            i=i+1;
        end
    end
    tline = fgetl(fid);
end
species(i:end) = [];
isfixed(i:end) = [];
fclose(fid);
end

function [J, photo_calls] = construct_mechanism(eqn_file, species)
J = cell(numel(species),1);
J(:) = {0};

fid = fopen(eqn_file);
tline = fgetl(fid);
photo_calls = {};
while ischar(tline)
    % skip comment lines or "dummy" reactions
    if ismember('#',tline) || ~isempty(regexp(tline,'JUNK','once'))
    else
        [products, product_coeff, reactants, reactant_coeff, k, isphoto] = read_wrf_mech_line(tline);
        
        % Construct the derivative function and add it to any existing such
        % functions in the mechanism for this species.
        rr = ismember(species, reactants);
        if sum(rr) ~= numel(reactants)
            nf = ~ismember(reactants, species);
            error('mechanism_parse:unknown_reactant','The reactant %s cannot be identified in the species list',strjoin(reactants(nf),', '));
        end
        
        if isa(k,'function_handle')
            k = func2str(k);
        elseif isphoto
            % We'll need this to know which function handles to construct
            % in the real mechanism.
            photo_calls{end+1} = k;
            k = sprintf('j(%d)',numel(photo_calls)); % This will call the proper element of vector j which must be passed to the jacobian and contains the photolysis rates. 
        end
        
        f = sprintf('%s * prod(c(%s) .^ %s)', k, mat2str(find(rr)), mat2str(reactant_coeff));
        
        % Add this for each reactant
        xx = find(ismember(species, reactants));
        if numel(xx) ~= numel(reactants)
            nf = ~ismember(reactants, species);
            error('mechanism_parse:unknown_reactant','The reactant %s cannot be identified in the species list',strjoin(reactants(nf),', '));
        end
        for j=1:numel(xx)
            if isa(J{xx(j)},'function_handle')
                J{xx(j)} = eval(sprintf('@(c,j,TEMP,C_M) %f * %s + %s', -reactant_coeff(j), f, func2str(J{xx(j)})));
            else
                J{xx(j)} = eval(sprintf('@(c,j,TEMP,C_M) %f * %s', -reactant_coeff(j), f));
            end
        end
        
        % Add this for each product
        xx = find(ismember(species, products));
        if numel(xx) ~= numel(products)
            nf = ~ismember(products, species);
            error('mechanism_parse:unknown_reactant','The product %s cannot be identified in the species list',strjoin(products(nf),', '));
        end
        for j=1:numel(xx)
            if isa(J{xx(j)},'function_handle')
                curr_fun = func2str(J{xx(j)});
                curr_fun = strrep(curr_fun,'@(c,j,TEMP,C_M)','');
                J{xx(j)} = eval(sprintf('@(c,j,TEMP,C_M) %f * %s + %s', product_coeff(j), f, curr_fun));
            else
                J{xx(j)} = eval(sprintf('@(c,j,TEMP,C_M) %f * %s', product_coeff(j), f));
            end
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
end