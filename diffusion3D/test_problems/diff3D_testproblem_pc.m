% DIFF3D_TESTPROBLEM_PC (for Windows users) driver for 3D test problem setup
% IFISS scriptfile: DJS; 27 July 2022.
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

gohome
close all
clear variables, clear functions
fprintf('\nspecification of reference Poisson problem.')
fprintf('\nchoose specific example');
fprintf('\n     1  Cube domain, constant source')
fprintf('\n     2  Stair-shaped domain, constant source')
fprintf('\n     3  Borehole domain, constant source\n')
fprintf('\n    11  Cube domain, analytic trilinear    solution ')
fprintf('\n    12  Cube domain, analytic triquadratic solution')
fprintf('\n    13  Cube domain, analytic tricubic     solution\n')


sn = default('',1);

%
if sn==1
system('copy .\diffusion3D\test_problems\unit_rhs3D.m .\diffusion3D\specific_rhs3D.m');
system('copy .\diffusion3D\test_problems\zero_bc3D.m .\diffusion3D\specific_bc3D.m');
cube_diff
    
energy_norm_squared_approx = x_it'*Agal*x_it;
energy_norm_squared_ref   = 0.645391924647045;
fprintf('Energy norm squared of approximation is %4.4f\n', energy_norm_squared_approx);
fprintf('Energy norm squared of reference solution is %4.4f\n', energy_norm_squared_ref);
e = sqrt(energy_norm_squared_ref - energy_norm_squared_approx);
fprintf('  ||u_{ref}-u_{FEM}||_{A} = %6.3e\n',e);
    if qmethod == 1
        effectivity_index = errorest./e;
        fprintf('Effectivity index = %4.4f\n',effectivity_index);
    end
    
elseif sn==2
system('copy .\diffusion3D\test_problems\unit_rhs3D.m .\diffusion3D\specific_rhs3D.m');
system('copy .\diffusion3D\test_problems\zero_bc3D.m .\diffusion3D\specific_bc3D.m');
stair_diff
  
energy_norm_squared_approx = x_it'*Agal*x_it;
energy_norm_squared_ref = 0.296720625256252;  % energy norm of the reference solution
fprintf('Energy norm squared of approximation is %4.4f\n', energy_norm_squared_approx);
fprintf('Energy norm squared of reference solution is %4.4f\n', energy_norm_squared_ref);
e = sqrt(energy_norm_squared_ref - energy_norm_squared_approx);
fprintf('  ||u_{ref}-u_{FEM}||_{A} = %6.3e\n',e);
    if qmethod == 1
        effectivity_index = errorest./e;
        fprintf('Effectivity index = %4.4f\n',effectivity_index);
    end

elseif sn==3
system('copy .\diffusion3D\test_problems\unit_rhs3D.m .\diffusion3D\specific_rhs3D.m');
system('copy .\diffusion3D\test_problems\zero_bc3D.m .\diffusion3D\specific_bc3D.m');
bore_diff

elseif sn == 11
system('copy .\diffusion3D\test_problems\zero_rhs3D.m .\diffusion3D\specific_rhs3D.m');
system('copy .\diffusion3D\test_problems\trilinear_bc.m .\diffusion3D\specific_bc3D.m');
cube_diff
    
energy_norm_squared_approx = x_it'*A*x_it;
fprintf('Energy norm squared of approximation is  %4.4f\n', energy_norm_squared_approx);
true_energy_norm_squared = 128/3;
fprintf('Energy norm squared of exact solution is %4.4f\n', true_energy_norm_squared);
if qmethod == 1
    e = energy_norm_error(xyz,ev,x_it,sn,qmethod);
   else
   e = sqrt(energy_norm_squared_approx - true_energy_norm_squared);
   end
fprintf('  ||u_{ref}-u_{FEM}||_{A} = %6.3e\n',e);

elseif sn == 12
system('copy .\diffusion3D\test_problems\triquadratic_rhs.m .\diffusion3D\specific_rhs3D.m');
system('copy .\diffusion3D\test_problems\zero_bc3D.m .\diffusion3D\specific_bc3D.m');
cube_diff
        
energy_norm_squared_approx = x_it'*Agal*x_it;
fprintf('Energy norm squared of approximation is %4.4f\n', energy_norm_squared_approx);
true_energy_norm_squared = 2048/225;
fprintf('Energy norm squared of exact solution is %4.4f\n', true_energy_norm_squared);
e = sqrt(true_energy_norm_squared - energy_norm_squared_approx);
fprintf('  ||u_{ref}-u_{FEM}||_{A} = %6.3e\n',e);
    if qmethod == 1
        effectivity_index = errorest./e;
        fprintf('Effectivity index = %4.4f\n',effectivity_index);
    end

elseif sn == 13
system('copy .\diffusion3D\test_problems\tricubic_rhs.m .\diffusion3D\specific_rhs3D.m');
system('copy .\diffusion3D\test_problems\zero_bc3D.m .\diffusion3D\specific_bc3D.m');
cube_diff
energy_norm_squared_approx = x_it'*A*x_it;
fprintf('Energy norm squared of approximation is %4.4f\n', energy_norm_squared_approx);
true_energy_norm_squared = 2048/18375;
fprintf('Energy norm squared of exact solution is %4.4f\n', true_energy_norm_squared);
e = sqrt(true_energy_norm_squared - energy_norm_squared_approx);
fprintf('  ||u_{ref}-u_{FEM}||_{A} = %6.3e\n',e);
    if qmethod == 1
        effectivity_index = errorest./e;
        fprintf('Effectivity index = %4.4f\n',effectivity_index);
    end

else
    error('Reference problem datafile not found!')
end


