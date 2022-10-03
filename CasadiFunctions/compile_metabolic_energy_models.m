function [] = compile_metabolic_energy_models(S)
% --------------------------------------------------------------------------
% compile_metabolic_energy_models
%   This function compiles the metabolic energy model functions into a .dll.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%
% OUTPUT:
%   - (this fuction returns no outputs) -
%   * 
% 
% Original author: Lars D'Hondt
% Original date: 28/Sept/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------



import casadi.*

%% f_getMetabolicEnergySmooth2004all
% L. J. Bhargava, M. G. Pandy, and F. C. Anderson, "A phenomenological model
%   for estimating metabolic energy consumption in muscle contraction,” Journal 
%   of biomechanics, vol. 37, no. 1, pp. 81–88, 2004.

% parameters
act_SX          = SX.sym('act_SX',1,1); % Muscle activations
exc_SX          = SX.sym('exc_SX',1,1); % Muscle excitations
lMtilde_SX      = SX.sym('lMtilde_SX',1,1); % N muscle fiber lengths
vM_SX           = SX.sym('vM_SX',1,1); % Muscle fiber velocities
Fce_SX          = SX.sym('FT_SX',1,1); % Contractile element forces
Fpass_SX        = SX.sym('FT_SX',1,1); % Passive element forces
Fiso_SX         = SX.sym('Fiso_SX',1,1); % N forces (F-L curve)
musclemass_SX   = SX.sym('musclemass_SX',1,1); % Muscle mass
pctst_SX        = SX.sym('pctst_SX',1,1); % Slow twitch ratio
modelmass_SX    = SX.sym('modelmass_SX',1); % Model mass
b_SX            = SX.sym('b_SX',1); % Parameter determining tanh smoothness
FMo_SX          = SX.sym('FMo_SX',1,1); % Optimal fibre force

% Bhargava et al. (2004)
[energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2004all(exc_SX,act_SX,...
    lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,pctst_SX,Fiso_SX,...
    FMo_SX,modelmass_SX,b_SX);

SX_in = {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,pctst_SX,...
    Fiso_SX,FMo_SX,modelmass_SX,b_SX};
SX_out = {energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,energy_model_sm_SX};
for i=1:length(SX_in)
    name_in{i} = SX_in{i}.name;
end

f_getMetabolicEnergySmooth2004all = Function('f_getMetabolicEnergySmooth2004all',...
    SX_in,SX_out,name_in,{'Edot','Adot','Mdot','Sdot','Wdot','Edot+basal'});


%% Generate C code
tmp_dir = pwd;
cd(fullfile(S.misc.main_path,'CasadiFunctions'))

% create code generator
CG = CodeGenerator('metabolic_energy_models.c');
% add functions add Jacobian
CG.add(f_getMetabolicEnergySmooth2004all);
CG.add(f_getMetabolicEnergySmooth2004all.jacobian());

% generate C code
CG.generate([fullfile(S.misc.main_path,'CasadiFunctions') '/']);

%% Compile
options.flags = {'/O2'};
options.folder = pwd;
options.cleanup = false;
options.name = 'metabolic_energy_models';
options.temp_suffix = false;
Ip = Importer('metabolic_energy_models.c','shell',options);


cd(tmp_dir)






