function [c_G,c_A,c_S,c_1] = amg_grid_complexity(amg_grid)
%AMG_GRID_COMPLEXITY checks grid/operator complexity and stencil size
%
% [c_G,c_A,c_S,c_1] = amg_grid_complexity(amg_grid);
%
% Input:
%           amg_grid    structure array generated by amg_grids_setup3D
%
% Outputs:
%       c_G       Grid complexity (sum of n_i / n_1)
%       c_A       Operator complexity (sum of nnz(A_i) / nnz(A_1)
%       c_S       Average stencil size (sum (nnz(A_i)/n_i) / L)
%       c_1       Average stencil size original problem (finest level)
%
% IFISS function: CP; 12th August 2022; DJS 24 August 2022
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

% Total no. of levels, including finest
L=max(size(amg_grid));
fprintf('\n total number of levels (L) is %i',L)

s=zeros(L,1); y=s; z=s;
for i=1:L
    s(i)=max(size(amg_grid(i).A));
    t(i)=nnz(amg_grid(i).A);
    z(i)=t(i)/s(i);
end

c_G=sum(s)/s(1);
c_A=sum(t)/t(1);
c_S=sum(z)/L;
c_1=z(1);

fprintf('\nstarting stencil size (c_1) is %2.2f', c_1)
fprintf('\n      grid complexity (c_G) is %2.2f', c_G)
fprintf('\n  operator complexity (c_A) is %2.2f', c_A)
fprintf('\n average stencil size (c_S) is %2.2f \n\n', c_S)



