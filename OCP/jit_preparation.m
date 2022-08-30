function [S] = jit_preparation(S)
% --------------------------------------------------------------------------
% jit_preparation
%   This functions makes preparations for using just-in-time compiling
%   functionality of CasADi
%   (https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F)
%   1) add settings for jit 
%   2) copy required files to dedicated folder
% 
% INPUT:
%   - S -
%   * setting structure S
%
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: Lars D'Hondt
% Original date: 29/August/2022
%
% Last edit by:
% Last edit date: 
% --------------------------------------------------------------------------



%% Prepare working dir
% Create working dir for jit
work_dir_name = ['tmp_jit_' S.post_process.result_filename];
work_dir_path = fullfile(S.misc.subject_path,work_dir_name);
S.solver.jit_work_dir_path = work_dir_path;

if ~isfolder(work_dir_path)
    mkdir(work_dir_path);
else
    warning(['"' work_dir_path '" already exists, just-in-time compiling could have issues'])
end

% copy .dll to working dir
copyfile(fullfile(S.misc.subject_path,S.misc.external_function),work_dir_path);

% backwards compatibility .lib 
base_name = S.misc.external_function(1:end-4);
lib_name = [base_name '.lib'];
lib_path = fullfile(S.misc.subject_path,lib_name);
if ~isfile(lib_path)
    pathLib = fullfile(S.misc.main_path,'opensimAD','install-ExternalFunction',base_name,'lib',lib_name);
    copyfile(pathLib,lib_path);
end

% copy .lib to working dir
copyfile(lib_path,work_dir_path);

S.solver.jit_libs{1} = lib_name;

%% Settings
% parallel computing
if strcmp(S.solver.parallel_mode,'openmp')
    % to use openmp, it needs to be added to the compiler
    if ~any(strcmp(S.solver.jit_compiler_flags(:),{'/openmp'}))
        S.solver.jit_compiler_flags{end+1} = '/openmp';
    end
    % set the max number of threads openmp can use
    setenv('OMP_NUM_THREADS',num2str(S.solver.N_threads));
end


cd(S.solver.jit_work_dir_path)







