%INSTALL_PC3D sets up IFISS 3D on non-UNIX computer
%   IFISS scriptfile: DJS; 19 August 2023
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

if strncmp(computer,'MAC2',4)
   fprintf('\nInstalling MAC-OS specific files must be done manually:\n')
%
elseif strncmp(computer,'PCWIN',3)
   gohome;
   fprintf('\nInstalling PC-OS specific files.\n')
%
delete .\diffusion3D\test_problems\diff_testproblem.m
system('copy .\diffusion3D\test_problems\diff_testproblem_pc.m .\diffusion3D\test_problems\diff_testproblem.m');
end
fprintf('\nAll installed!\n')
