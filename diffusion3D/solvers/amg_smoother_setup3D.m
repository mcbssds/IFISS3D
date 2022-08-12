function [smoother_data] = amg_smoother_setup3D(grid_data, smoother_params)
%AMG_SMOOTHER_SETUP3D generates smoother data for AMG
%
% smoother_data = amg_smoother_setup3D(grid_data, smoother_params)
%
% For use with 'amg_v_cycle'. Takes grid_data generated using amg_setup 
% and creates smoother_data according to structure smoother_params.
% 
% Inputs:
%
%   grid_data        structure array generated using 'amg_grids_setup'
%   smoother_params  structure array generated using 'amg_smoother_params'
%                    where smoother_params(n) contains smoother parameters
%                    for level==n
%
% Output:
%
%  smoother_data     structure array containing smoother information
%
% IFISS function: CP 12th August 2022
% Copyright (c) 2007 by J. Boyle

size = length(grid_data);

for level = 1:size
    % get the smoother information for this level 
    smoother = smoother_params(level);
    smoother_data(level).smoother_params = smoother;
    % point Gauss Seidel   
    if strcmp(smoother.type,'PGS')
        % set up for forward sweep
        if strcmp(smoother.pre, 'forward') | strcmp(smoother.post, 'forward')
            smoother_data(level).LD = tril(grid_data(level).A);     % part containing D & L
        end
        % set up backward sweep
        if strcmp(smoother.pre, 'backward') | strcmp(smoother.post, 'backward')     
            smoother_data(level).UD = triu(grid_data(level).A);     % part containing D & U
        end
    % point damped Jacobi sweep 
    elseif strcmp(smoother.type,'PDJ')
        omega = smoother.damping;   % relaxation factor for damped Jacobi
        n = length(grid_data(level).A);
        smoother_data(level).D = (1/omega)*spdiags(diag(grid_data(level).A), 0, n, n);    
     end
end
