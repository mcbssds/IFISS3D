function batchmode3D(testproblem)
%BATCHMODE3D enables batch processing for IFISS testproblem in 3D
%   batchmode3D(testproblem);
%   input
%          testproblem  character string must have the form
%                        "*_batch".m where "*" begins with "3D"
%   side effect
%          If batchmode terminates prematurely because of an error or
%          execution of cntl-C, interactive input with IFISS may not
%          work correctly.  This is fixed by typing "activemode".
%
%
%   IFISS function: DJS; 18 August 2020.
% Copyright (c)  2022  G.Papanikos, C.E. Powell, D.J. Silvester

global BATCH FID

% file containing input data 
batchfile=[testproblem,'_batch.m'];
[FID,message]=fopen(batchfile,'r');
% Error return on nonexistent or misnamed input file
if strcmp(message,'')~=1
   error(['INPUT FILE ERROR: ' message])
else
   disp(['Working in batch mode from data file ' batchfile])
end
if ~strncmp(testproblem,'3D',2),
    errmsg = 'INPUT FILE ERROR:\n';
    errmsg = [errmsg,'   Batch input filenames must have the form "*_batch.m"'];
    errmsg = [errmsg,' where "*" begins with\n'];
    errmsg = [errmsg,'   "3D" for generation of a 3D problem'];
    error('BATCH:err',errmsg);    
end  

% batch run
% switch to activate batch mode (off/on 0/1) (see "default.m")
BATCH=1;
% run driver
diff3D_testproblem
gohome, cd datafiles, save batchrun.mat

% switch back to interactive mode
activemode
return
