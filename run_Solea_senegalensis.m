%close all; 
global pets

pets = {'Solea_senegalensis'};
check_my_pet(pets); 

% Read about how to set estimation and output options (estim_options) on the online
% manual: http://www.debtheory.org/wiki/index.php?title=Run_file

estim_options('default'); 
estim_options('max_step_number',5e3); 
estim_options('max_fun_evals',5e3);  

estim_options('pars_init_method', 2); %1 to continue 2 read from pars_init
estim_options('results_output', 1); % -1 in the command window 
estim_options('method', 'nm');
estim_pars; 



 % see http://www.debtheory.org/wiki/index.php?title=Estimation_options
% filter 
% 1 - use filter (default); 
% 0 - do not; 
% 
% pars_init_method: 
% 0 - get initial values from automatized computation (default) 
% 1 - read initial estimates from .mat file (for continuation) 
% 2 - read initial estimates from pars_init file 
% 
% results_output 
% 0 only saves data to .mat - use this for (automatic) continuations
% 1 no saving to .mat file, prints results to html (1) or screen (-1) 
% 2 saves to .mat file, prints results to html (2) or screen (-2)
% 3 saves data to .mat file and graphs to .png files, prints results to html or screen
% 4 saves data to .mat file and graphs to .png files; prints results to html or screen;
%   prints html with implied properties
% 5 saves data to .mat file and graphs to .png files; prints results to html or screen; 
%   prints html with implied properties including related species
% method 

