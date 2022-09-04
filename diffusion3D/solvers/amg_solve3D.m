%AMG_SOLVE3D HSL_MI20 solution of predefined 3D diffusion problem
%IFISS scriptfile: DJS 4 September 2022
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

if exist('pde','var')==0
    error('Oops ... you need to set up a specific discrete problem first!'),
end

if exist('hsl_mi20_setup')~=2
    error('Oops ... you need to install hsl_mi20 functionality!'),
end

if pde==1   % 3D diffusion problem
    fprintf('solving discretised diffusion problem...\n')
    % select Krylov subspace method
    itmeth = 2; %-------- MINRES
    
    % set parameters
    tol = 1e-10;
    maxit = 100;
    
    % select preconditioner and construct it
    fprintf('fast AMG preconditioner\n');
    tic
    control = hsl_mi20_control;
    control.st_parameter=2/3
    inform = hsl_mi20_setup(Agal,control)
    fprintf('\nAMG setup done.\n'), toc

    % zero initial guess
    x0=zeros(size(fgal));
    tic
    [x_amg,flag,relres,iter,resvec] = ...
     minres(Agal,fgal,tol,maxit,'hsl_mi20_precondition');
     hsl_mi20_finalize;
     etoc = toc;
else
error('Oops ... undefined PDE problem'),
end

% print and plot results
if flag ==0
    % successful convergence
    fprintf('convergence in %3i iterations\n',iter)
    nr0=resvec(1);
    fprintf('\n    k  log10(||r_k||/||r_0||)   \n')
    for its=1:iter+1,
        fprintf('%5i %16.4f \n', its-1, log10(resvec(its)/nr0));
    end
    fprintf('Bingo!\n')
    fprintf('\n  %9.4e seconds\n\n\n',etoc)
    %%% plot residuals
resplot(resvec)
else
    nr0=resvec(1);
    fprintf('\n    k  log10(||r_k||/||r_0||)   \n')
    for its=1:iter+1,
        fprintf('%5i %16.4f \n', its-1, log10(resvec(its)/nr0));
    end
    fprintf('iteration aborted! Iteration returned with flag equal to  %2i \n',flag)
    %%% plot residuals
resplot(resvec)
end
