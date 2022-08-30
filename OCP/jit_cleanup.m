function [] = jit_cleanup(S)
% --------------------------------------------------------------------------
% jit_cleanup
%   This functions removes temporary files using just-in-time compiling
%   functionality of CasADi
%   (https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F)
% 
% INPUT:
%   - S -
%   * setting structure S
% 
% OUTPUT:
%   - this function has no outputs -
% 
% Original author: Lars D'Hondt
% Original date: 29/August/2022
%
% Last edit by:
% Last edit date: 
% --------------------------------------------------------------------------

cd(S.misc.main_path)

% Remove working dir
work_dir_name = ['tmp_jit_' S.post_process.result_filename];
work_dir_path = fullfile(S.misc.subject_path,work_dir_name);
if isfolder(work_dir_path)
    dir_is_removed = 0;
    try
        dir_is_removed = rmdir(work_dir_path,'s');
    catch err
        if ~dir_is_removed
            warning(['Unable to remove "' work_dir_name '" from "',...
                S.misc.subject_path '", because:' err.message])
        end
    end
end

