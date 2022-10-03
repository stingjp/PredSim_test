function [func, lib] = compile_CasADI_Function(fun,work_dir,libs,varargin)
% --------------------------------------------------------------------------
% compile_CasADI_Function
%   This functions compiles a CasADi Function into a .dll, and loads it
%   back as external function. 
% 
% INPUT:
%   - fun -
%   * CasADi Function
%
%   - work_dir -
%   * path to folder where the generated files are stored
%
%   - libs -
%   % .lib dependencies
%
% OUTPUT:
%   - func -
%   * external function
%
%   - lib -
%   * .lib file, needed for linking the code
% 
% Original author: Lars D'Hondt
% Original date: 29/August/2022
%
% Last edit by:
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*

% name of function
if length(varargin)>=1
    fun_name = varargin{1};
else
    fun_name = fun.name;
end

%% Generate C code
% create code generator
CG = CodeGenerator([fun_name '.c']);
% add fun
CG.add(fun);
% add Jacobian of fun
CG.add(fun.jacobian())
% generate C code
CG.generate([work_dir '/']);

%% Compile C code into dll

options.flags = {'/Ox'};
options.folder = work_dir;
if ~isempty(libs)
    options.linker_flags = libs;
end
options.cleanup = false;
options.name = fun_name;
options.temp_suffix = false;
Ip = Importer([fun_name '.c'],'shell',options);

func = external(fun.name,Ip);

%% Get name of .lib file
LibFiles = dir(fullfile(work_dir,'*.lib'));
lib = LibFiles(1).name;







end