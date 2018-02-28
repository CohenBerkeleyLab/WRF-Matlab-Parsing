classdef ReactionNetwork < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        reactions;
    end
    
    methods
        function obj = ReactionNetwork(mechanism_file)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.reactions = obj.ParseMechFile(mechanism_file);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods(Static, Access = private)
        function reactions = ParseMechFile(mech_file)
            % This should eventually also parse the species file, if only
            % to make sure that no erroneous species are expected in the
            % equation file.
            reactions = cell(1000,1); % will clean up extra entries at the end;
            fid = fopen(mech_file);
            try
                tline = fgetl(fid);
                i_rxn = 1;
                while ischar(tline)
                    if strcmp(tline(1),'#')
                        % Skip comment lines
                        tline = fgetl(fid);
                        continue
                    end
                    fprintf('Parsing line: %s\n', tline);
                    reactions{i_rxn} = KPP_OOP.Reaction(tline);
                    i_rxn = i_rxn + 1;
                    tline = fgetl(fid);
                end
            catch err
                fclose(fid);
                rethrow(err);
            end
            fclose(fid);
            reactions(i_rxn:end) = [];
        end
    end
end

