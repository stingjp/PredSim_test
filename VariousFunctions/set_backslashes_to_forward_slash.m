function [path_with_forward_slash] = set_backslashes_to_forward_slash(path_with_backslash)

path_tmp = path_with_backslash;
idx_bs = strfind(path_tmp,'\');
path_tmp(idx_bs) = '/';
path_with_forward_slash = path_tmp;









