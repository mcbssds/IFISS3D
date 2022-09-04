%IT_SOLVE3D iterative solution of predefined 3D diffusion problem
%IFISS scriptfile: CP 12 August 2022; DJS  24 August 2022
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester
global amg_grid amg_smoother   
% Declare global variables
if exist('pde','var')==0
    error('Oops ... you need to set up a specific discrete problem first!'),
end

if pde==1   % 3D diffusion problem
    fprintf('solving discretised 3D diffusion problem...\n')
    % select Krylov subspace method
    itmeth = default('PCG/MINRES? 1/2 (default PCG)',1);
    
    % set parameters
    tol = 1e-10;
    maxit = default('maximum number of iterations? (default 100)',100);
    
    % select preconditioner and construct it
    fprintf('preconditioner:\n');
    fprintf('   0  none\n');
    fprintf('   1  diagonal\n');
    fprintf('   2  incomplete Cholesky\n');
    fprintf('   3  algebraic multigrid\n');
    precon = default('default is incomplete Cholesky ',2);
    if precon==0     % none
        M1=[]; M2=[];
    elseif precon==1 % diagonal
        D=diag(diag(Agal)); M1=sqrt(D); M2=M1;
    elseif precon==2 % incomplete Cholesky
        tic,M1 = ichol(Agal); M2=M1'; fprintf('\nsetup done.\n'), toc
    elseif precon==3 % AMG
        % uses global variables amg_grid amg_smoother
        regrid = default('reuse AMG grid sequence? yes/no 1/2 (default no)',2);
        if regrid~=1, tic, amg_grid = amg_grids_setup3D(Agal,1);
        fprintf('\nsetup done.\n'), toc
        amg_grid_complexity(amg_grid); end
        plot_mg = default('plot AMG grid sequence? yes/no 1/2 (default no)',2);
        if plot_mg==1, amg_coarsen_plot3D(amg_grid, xyz); end
        smoothopt = default('PDJ/PGS smoother? 1/2 (point Gauss Seidel)',2);
        if smoothopt==1
            fprintf('point damped Jacobi smoothing ..\n')
            smoother_params = amg_smoother_params3D(amg_grid, 'PDJ', 2);
        else
            fprintf('point Gauss-Seidel smoothing ..\n')
            smoother_params = amg_smoother_params3D(amg_grid, 'PGS', 2);
        end
        amg_smoother = amg_smoother_setup3D(amg_grid, smoother_params);
    else
        error('Oops ... invalid preconditioner')
    end


    % zero initial guess
    x0=zeros(size(fgal));
    tic
    if itmeth==1 %---- PCG
        if precon < 3
            fprintf('\nPCG iteration ...\n');
            [x_it,flag,relres,iter,resvec] = pcg(Agal,fgal,tol,maxit,M1,M2,x0);
        else
            [x_it,flag,relres,iter,resvec] = ...
             pcg(Agal,fgal,tol,maxit, @amg_v_cycle, [], x0, amg_grid, amg_smoother); 
        end
    elseif itmeth==2 %---- MINRES
        if precon < 3
            fprintf('\nMINRES iteration ...\n');
            [x_it,flag,relres,iter,resvec] = minres(Agal,fgal,tol,maxit,M1,M2,x0);
        else
            [x_it,flag,relres,iter,resvec] = ...
              minres(Agal,fgal,tol,maxit,@amg_v_cycle,[],x0, amg_grid,amg_smoother);
        end
    else
        error('Oops ... invalid iterative method!')
    end
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
