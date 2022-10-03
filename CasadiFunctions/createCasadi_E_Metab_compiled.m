function [f_getMetabolicEnergySmooth2004all] = createCasadi_E_Metab_compiled(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_E_Metab
%   Function to create Casadi functions for metabolic energy.
%   Uses compiled metabolic energy model (.dll).
%   
% INPUT:
%   - S -
%   * setting structure S
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_getMetabolicEnergySmooth2004all -
%   * Casadi functions for metabolic energy based on
%   L. J. Bhargava, M. G. Pandy, and F. C. Anderson, "A phenomenological model
%   for estimating metabolic energy consumption in muscle contraction,” Journal 
%   of biomechanics, vol. 37, no. 1, pp. 81–88, 2004.
% 
% Original author: Ines Vandekerckhove, KU Leuven
% Original date: 01-12-2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*
NMuscle = model_info.muscle_info.NMuscle;
FMo = struct_array_to_double_array(model_info.muscle_info.parameters,'FMo');

% Compile functions if needed
dll_path = fullfile(S.misc.main_path,'CasadiFunctions','metabolic_energy_models.dll');
test_dll = isfile(dll_path);
test_lib = isfile(fullfile(S.misc.main_path,'CasadiFunctions','metabolic_energy_models.lib'));
if ~(test_lib && test_dll)
    compile_metabolic_energy_models(S);
end

% Load functions
F_getMetabolicEnergySmooth2004all = external('f_getMetabolicEnergySmooth2004all',dll_path);

% Map to number of muscles
F_getMetabolicEnergySmooth2004all_map = F_getMetabolicEnergySmooth2004all.map(NMuscle);


%% Metabolic energy models
act_MX          = MX.sym('act_MX',NMuscle); % Muscle activations
exc_MX          = MX.sym('exc_MX',NMuscle); % Muscle excitations
lMtilde_MX      = MX.sym('lMtilde_MX',NMuscle); % N muscle fiber lengths
vM_MX           = MX.sym('vM_MX',NMuscle); % Muscle fiber velocities
Fce_MX          = MX.sym('Fce_MX',NMuscle); % Contractile element forces
Fpass_MX        = MX.sym('Fpass_MX',NMuscle); % Passive element forces
Fiso_MX         = MX.sym('Fiso_MX',NMuscle); % N forces (F-L curve)
musclemass_MX   = MX.sym('musclemass_MX',NMuscle); % Muscle mass
pctst_MX        = MX.sym('pctst_MX',NMuscle); % Slow twitch ratio
modelmass_MX    = MX.sym('modelmass_MX',1); % Model mass
b_MX            = MX.sym('b_MX',1); % Parameter determining tanh smoothness
% Bhargava et al. (2004)

[energy_total_sm_MX,Adot_sm_MX,Mdot_sm_MX,Sdot_sm_MX,Wdot_sm_MX,energy_model_sm_MX] =...
    F_getMetabolicEnergySmooth2004all_map(exc_MX,act_MX,...
    lMtilde_MX,vM_MX,Fce_MX,Fpass_MX,musclemass_MX,pctst_MX,Fiso_MX,...
    FMo,repmat(modelmass_MX,NMuscle,1),repmat(b_MX,NMuscle,1));

% adapt energy model to compensate evaluating function 1 muscle at a time
energy_model_sm_MX = energy_model_sm_MX(1)-energy_total_sm_MX(1) + sum(energy_total_sm_MX);

f_getMetabolicEnergySmooth2004all = ...
    Function('f_getMetabolicEnergySmooth2004all_all_muscles',...
    {exc_MX,act_MX,lMtilde_MX,vM_MX,Fce_MX,Fpass_MX,musclemass_MX,...
    pctst_MX,Fiso_MX,modelmass_MX,b_MX},{energy_total_sm_MX,...
    Adot_sm_MX,Mdot_sm_MX,Sdot_sm_MX,Wdot_sm_MX,sum(energy_model_sm_MX)});

end