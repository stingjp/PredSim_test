function [f_ForceEquilibrium_FtildeState_all_tendon,f_FiberLength_TendonForce_tendon,...
    f_FiberVelocity_TendonForce_tendon,f_lT_vT] = createCasadi_ContractDynam_compiled(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_ContractDynam_compiled
%   Function to create Casadi functions for muscle contraction dynamics.
%   Uses compiled muscle contraction dynamics (.dll).
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_ForceEquilibrium_FtildeState_all_tendon -
%   * function for force equilibrium between muscle fiber and tendon
%
%   - f_FiberLength_TendonForce_tendon -
%   * function to compute fiber length of a muscle
%
%   - f_FiberVelocity_TendonForce_tendon -
%   * function to compute fiber lenghtening velocity of a muscle
%
%   - f_lT_vT -
%   * function to compute tendon length and lengthening velocity of a tendon
% 
% Original author: Lars D'Hondt
% Original date: 27/Sept/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*
N_muscles = model_info.muscle_info.NMuscle;

% Compile functions if needed
dll_path = fullfile(S.misc.main_path,'CasadiFunctions','muscle_contraction_dynamics.dll');
test_dll = isfile(dll_path);
test_lib = isfile(fullfile(S.misc.main_path,'CasadiFunctions','muscle_contraction_dynamics.lib'));
if ~(test_lib && test_dll)
    compile_muscle_contraction_dynamics(S);
end

% Load functions
F_ForceEquilibrium_FtildeState_all_tendon = external('f_ForceEquilibrium_FtildeState_all_tendon',dll_path);
F_FiberLength_TendonForce_tendon = external('f_FiberLength_TendonForce_tendon',dll_path);
F_FiberVelocity_TendonForce_tendon = external('f_FiberVelocity_TendonForce_tendon',dll_path);
F_lT_vT = external('f_lT_vT',dll_path);

% Map to number of muscles
F_ForceEquilibrium_FtildeState_all_tendon_map = F_ForceEquilibrium_FtildeState_all_tendon.map(N_muscles);
F_FiberLength_TendonForce_tendon_map = F_FiberLength_TendonForce_tendon.map(N_muscles);
F_FiberVelocity_TendonForce_tendon_map = F_FiberVelocity_TendonForce_tendon.map(N_muscles);
F_lT_vT_map = F_lT_vT.map(N_muscles);


% tensions = struct_array_to_double_array(model_info.muscle_info.parameters,'specific_tension');
FMo = struct_array_to_double_array(model_info.muscle_info.parameters,'FMo')';
lMo = struct_array_to_double_array(model_info.muscle_info.parameters,'lMo')';
vMmax = struct_array_to_double_array(model_info.muscle_info.parameters,'vMmax')';
lTs = struct_array_to_double_array(model_info.muscle_info.parameters,'lTs')';
alphao = struct_array_to_double_array(model_info.muscle_info.parameters,'alphao')';
tendon_stiff = struct_array_to_double_array(model_info.muscle_info.parameters,'tendon_stiff')';
tendon_stiff_shift = struct_array_to_double_array(model_info.muscle_info.parameters,'tendon_stiff_shift')';
dampingCoefficient = ones(N_muscles,1)*S.misc.dampingCoefficient;
muscle_pass_stiff_shift = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_pass_stiff_shift')';
muscle_pass_stiff_scale = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_pass_stiff_scale')';
muscle_strength = struct_array_to_double_array(model_info.muscle_info.parameters,'muscle_strength')';

%% Muscle contraction dynamics
% Function for Hill-equilibrium
FTtilde     = MX.sym('FTtilde',N_muscles); % Normalized tendon forces
a           = MX.sym('a',N_muscles); % Muscle activations
dFTtilde    = MX.sym('dFTtilde',N_muscles); % Time derivative tendon forces
lMT         = MX.sym('lMT',N_muscles); % Muscle-tendon lengths
vMT         = MX.sym('vMT',N_muscles); % Muscle-tendon velocities
tensions    = MX.sym('tension',N_muscles); % Tensions

% Hilldiff    = SX(N_muscles,1); % Hill-equilibrium
% FT          = SX(N_muscles,1); % Tendon forces
% Fce         = SX(N_muscles,1); % Contractile element forces
% Fiso        = SX(N_muscles,1); % Normalized forces from force-length curve
% vMmax       = SX(N_muscles,1); % Maximum contraction velocities
% massM       = SX(N_muscles,1); % Muscle mass
% Fpass       = SX(N_muscles,1); % Passive element forces


[Hilldiff,FT,Fce,Fpass,Fiso,vMmax_MX,massM] = ...
    F_ForceEquilibrium_FtildeState_all_tendon_map(a,FTtilde,dFTtilde,lMT,vMT,FMo,...
    lMo,lTs,alphao,vMmax,tensions,tendon_stiff,tendon_stiff_shift,...
    dampingCoefficient,muscle_pass_stiff_shift,muscle_pass_stiff_scale,muscle_strength);

f_ForceEquilibrium_FtildeState_all_tendon = ...
    Function('f_forceEquilibrium_FtildeState_all_tendon_all_muscles',{a,FTtilde,...
    dFTtilde,lMT,vMT,tensions},{Hilldiff',FT',Fce',Fpass',Fiso',vMmax_MX',massM'},...
    {'a','FTtilde','dFTtilde','lMT','vMT','tension_SX'},...
    {'Hilldiff','FT','Fce','Fpass','Fiso','vMmax','massM'});

% Function to get (normalized) muscle fiber lengths

[lM,lMtilde] = F_FiberLength_TendonForce_tendon_map(FTtilde,...
    lMo,lTs,alphao,lMT,tendon_stiff,tendon_stiff_shift);

f_FiberLength_TendonForce_tendon = Function(...
    'f_FiberLength_Ftilde_tendon_all_muscles',{FTtilde,lMT},{lM',lMtilde'},...
    {'FTtilde','lMT'},{'lM','lMtilde'});

% Function to get (normalized) muscle fiber velocities
[vM,vMtilde] = F_FiberVelocity_TendonForce_tendon_map(FTtilde,...
    dFTtilde,lMo,lTs,alphao,vMmax,lMT,vMT,tendon_stiff,tendon_stiff_shift);

f_FiberVelocity_TendonForce_tendon = Function(...
    'f_FiberVelocity_Ftilde_tendon_all_muscles',{FTtilde,dFTtilde,lMT,vMT},...
    {vM',vMtilde'},{'FTtilde','dFTtilde','lMT','vMT'},{'vM','vMtilde'});

% Function to get tendon kinematics
[lT,vT] = F_lT_vT_map(FTtilde,...
    dFTtilde,lMo,lTs,alphao,vMmax,lMT,vMT,tendon_stiff,tendon_stiff_shift);

f_lT_vT = Function('f_lT_vT_all_muscles',{FTtilde,dFTtilde,lMT,vMT},...
    {lT',vT'},{'FTtilde','dFTtilde','lMT','vMT'},{'lT','vT'});









end