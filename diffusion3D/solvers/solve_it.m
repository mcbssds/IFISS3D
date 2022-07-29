%SOLVE_IT iterative solver for high-dimensional Galerkin system
% IFISS scriptfile: DJS 29 July 2022
% Copyright (c) 2022 G. Papanikos C.E. Powell, D.J. Silvester

if exist('pde','var')==0,
   error('You need to set up a specific discrete problem first!'), 
end

%----- preset parameters
tol = 1e-10;
maxit = 199;

if pde==1,  %------ Poisson problem
itmeth = 2; %------ MINRES
%------ incomplete cholesky
M1 = ichol(Agal); M2=M1';
   x0=zeros(size(fgal));
   tic
   if itmeth==1,
      fprintf('\nPCG iteration ...');
      [x_it,flag,relres,iter,resvec] = pcg(Agal,fgal,tol,maxit,M1,M2,x0);
   elseif itmeth==2,
      fprintf('\nMINRES iteration ...');
      [x_it,flag,relres,iter,resvec] = minres(Agal,fgal,tol,maxit,M1,M2,x0);
   end
  etoc = toc;
end

%------ print and plot results
if flag ==0,
   % successful convergence
   fprintf('convergence in %3i iterations\n',iter)
else
   % aborted iteration
   nr0=resvec(1);
   fprintf('\n    k  log10(||r_k||/||r_0||)   \n')
   for its=1:iter+1,
      fprintf('%5i %16.4f \n', its-1, log10(resvec(its)/nr0));
   end
   fprintf('iteration aborted! Iteration returned with flag equal to  %2i \n',flag)
   resplot(resvec)
end
