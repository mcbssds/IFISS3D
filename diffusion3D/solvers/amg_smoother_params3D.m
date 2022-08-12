function [smoother_params] = amg_smoother_params3D(amg_grid, smoother_type, no_sweeps)
%AMG_SMOOTHER_PARAMS3D generates structure specifying AMG smoother parameters
%
% smoother_params = amg_smoother_params3D(amg_grid, smoother_type, no_sweeps)
%
% Need to call amg_smoother_setup first!
%
% Inputs:
%
%   amg_grid        structure generated using 'amg_grids_setup'
%   smoother_type   string to specify smoother type; one of 
%
%                   'PGS' - point Gauss Seidel
%                   'PDJ' - point damped Jacobi
%
%   no_sweeps       no. of pre- and post-smoothing steps to apply
%
% Output:
% 
% smoother_params:  structure array where smoother_params(n) contains
%                   smoother parameters for level==n. Components are:
%
%      smoother_params(n).type      smoother type
%      smoother_params(n).nsweeps   no. of times smoother is applied (pre and post)
%      smoother_params(n).pre       pre-sweep direction ('forward' or 'backward')
%      smoother_params(n).post      post-sweep direction ('forward' or 'backward')
%      smoother_params(n).damping   damping for PDJ, default value is 0.5
%
%  IFISS function: CP; 12th August 2022. 
%  Copyright (c) 2007 by J. Boyle

no_levels = length(amg_grid);
smoother_params = [];
damping = 0.5;  % damping for PDJ - NEEDS FIXING IN 3D/USER DEFINED INPUT??

for level = 1:no_levels
    smoother_params(level).nsweeps = no_sweeps;
    smoother_params(level).type = smoother_type;
    smoother_params(level).pre = 'forward';    
    smoother_params(level).post = 'backward';        
    if  strcmp(smoother_type,'PDJ')
        smoother_params(level).damping = damping;
    end
end


