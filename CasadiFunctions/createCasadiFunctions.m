function [f_casadi] = createCasadiFunctions(S,model_info)
% --------------------------------------------------------------------------
%createCasadiFunctions.m
%   Overview function from which al casadi functions are created
% 
% INPUT:
%   - S -
%   * setting structure S
%  
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% OUTPUT:
%   - f_casadi -
%   * Struct containing all casadi functions.
% 
% Original author: Tom Buurke
% Original date: 02/12/2021
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

%% Create generic casadi functions
f_casadi = createCasadi_GenHelper(S,model_info);

%% Create Casadi functions for musculoskeletal geometry
f_casadi.lMT_vMT_dM = createCasadi_MSKGeometry(S,model_info);

%% Create Casadi functions for muscle contraction dynamics
if S.solver.mus_dyn_codegen == 0 || S.solver.mus_dyn_codegen == 2
    [forceEquilibrium_FtildeState_all_tendon, FiberLength_TendonForce_tendon,...
        FiberVelocity_TendonForce_tendon,lT_vT] = createCasadi_ContractDynam(S,model_info);
    
    f_casadi.forceEquilibrium_FtildeState_all_tendon = forceEquilibrium_FtildeState_all_tendon;
    f_casadi.FiberLength_TendonForce_tendon = FiberLength_TendonForce_tendon;
    f_casadi.FiberVelocity_TendonForce_tendon = FiberVelocity_TendonForce_tendon;
    f_casadi.lT_vT = lT_vT;

end

% new functions, uses codegen, compiling, mapping
if S.solver.mus_dyn_codegen
    [f_ForceEquilibrium_FtildeState_all_tendon,f_FiberLength_TendonForce_tendon,...
        f_FiberVelocity_TendonForce_tendon,f_lT_vT] = createCasadi_ContractDynam_compiled(S,model_info);
    
    % test
    if S.solver.mus_dyn_codegen == 2
        addpath(fullfile(S.misc.main_path,'Tests'));
        compareCasADiFunctions(f_casadi.forceEquilibrium_FtildeState_all_tendon,f_forceEquilibrium_FtildeState_all_tendon,10);
        compareCasADiFunctions(f_casadi.FiberLength_TendonForce_tendon,f_FiberLength_TendonForce_tendon,10);
        compareCasADiFunctions(f_casadi.FiberVelocity_TendonForce_tendon,f_FiberVelocity_TendonForce_tendon,10);
        compareCasADiFunctions(f_casadi.lT_vT,f_lT_vT,10);
    end

    % use new function
    f_casadi.forceEquilibrium_FtildeState_all_tendon = f_ForceEquilibrium_FtildeState_all_tendon;
    f_casadi.FiberLength_TendonForce_tendon = f_FiberLength_TendonForce_tendon;
    f_casadi.FiberVelocity_TendonForce_tendon = f_FiberVelocity_TendonForce_tendon;
    f_casadi.lT_vT = f_lT_vT;

end

%% Create Casadi functions for passive torques
[f_casadi.PassiveStiffnessMoments,f_casadi.PassiveDampingMoments,f_casadi.LimitTorques,...
    f_casadi.AllPassiveTorques,f_casadi.AllPassiveTorques_cost] = createCasadi_PassiveMoments(S,model_info);

%% Create Casadi functions for activation dynamics
if model_info.ExtFunIO.jointi.nq.torqAct > 0
    [f_casadi.ActuatorActivationDynamics] = createCasadi_ActDynam(S,model_info);
end

%% Create Casadi functions for metabolic energy.
if S.solver.metab_codegen == 0 || S.solver.metab_codegen == 2
    [f_casadi.getMetabolicEnergySmooth2004all] = createCasadi_E_Metab(S,model_info);
end

if S.solver.metab_codegen
    [f_getMetabolicEnergySmooth2004all] = createCasadi_E_Metab_compiled(S,model_info);

    if S.solver.metab_codegen==2
        addpath(fullfile(S.misc.main_path,'Tests'));
        compareCasADiFunctions(f_casadi.getMetabolicEnergySmooth2004all,f_getMetabolicEnergySmooth2004all,10);
    end

    f_casadi.getMetabolicEnergySmooth2004all = f_getMetabolicEnergySmooth2004all;
end

%% Create Casadi function to get step length
if isfield(model_info.ExtFunIO,'origin') && ...
        isfield(model_info.ExtFunIO.origin,'calcn_r') && isfield(model_info.ExtFunIO.origin,'calcn_l') && ...
        ~isempty(model_info.ExtFunIO.origin.calcn_r) &&  ~isempty(model_info.ExtFunIO.origin.calcn_l)
    
    [f_casadi.f_getCalcnOriginInWorldFrame,f_casadi.f_getStepLength] = createCasadi_StepLength(S,model_info);
end


end