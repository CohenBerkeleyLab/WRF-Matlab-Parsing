function [ conc_out  ] = run_wrf_mech( init_cond )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% Initial conditions. Concentrations should be in molec. cm^{-3}. Each row
% of the cell array should have the name of the species (case insensitive)
% and the concentration.

Nair = 2e19; % number density of air at surface
T = 298;

if ~exist('init_cond','var')
    init_cond = {   'NO2', 0;...
                    'NO', 1e-9 * Nair;...
                    'M', Nair};
end

[ J, species, isfixed, photo_calls ] = parse_wrf_mech('r2smh-simple'); %#ok<ASGLU>

lon = -84.39;
lat = 33.775;
date_in = '2013-06-01';

nt = 10800; % number timesteps
dt = 1; % seconds
save_freq = 60; % how frequently to save out the concentrations
start_utc_time = 19; % hr utc
run_mode = 'normal';
%%%%%%%%%%%%%%%%%%
%%% Initialize %%%
%%%%%%%%%%%%%%%%%%

C = zeros(size(species)); % concentration vector
for a=1:size(init_cond,1)
    xx = strcmpi(init_cond{a,1},species);
    C(xx) = init_cond{a,2};
end

conc_out = nan(numel(C), ceil(nt/save_freq));

tuv_hr = nan;

%%%%%%%%%%%
%%% RUN %%%
%%%%%%%%%%%

switch run_mode
    case 'normal'
        % Run in regular mode
        for t=1:nt
            utc_time = start_utc_time + t*dt/(3600);
            mech_timestep;
        end
        
    case 'steady-state'
        % Run in SS mode
        utc_time = start_utc_time;
        while true
            mech_timestep;
        end
end


function mech_timestep
    % Should we call TUV again? Do so if we are closer to a new hour than
    % we were before.
    if round(utc_time) ~= tuv_hr
        tuv_hr = round(utc_time);
        j = call_tuv(photo_calls, date_in, tuv_hr, lon, lat, true);
    end
    dC = zeros(size(C));
    for i=1:numel(C)
        dC(i) = J{i}(C,j,T,Nair)*dt;
    end
    C = C + dC;
    
    % Save if requested
    if mod(t,save_freq) == 0
        conc_out(:,ceil(t/save_freq)) = C;
    end
end

end

