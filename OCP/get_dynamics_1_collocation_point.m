function [f_coll] = get_dynamics_1_collocation_point(S,model_info,f_casadi)

%% OCP: collocation equations
% Define CasADi variables for static parameters
tfk         = MX.sym('tfk'); % MX variable for final time
% Define CasADi variables for states
ak          = MX.sym('ak',NMuscle);
aj          = MX.sym('akmesh',NMuscle,d);
akj         = [ak aj];
FTtildek    = MX.sym('FTtildek',NMuscle);
FTtildej    = MX.sym('FTtildej',NMuscle,d);
FTtildekj   = [FTtildek FTtildej];
Qsk         = MX.sym('Qsk',nq.all);
Qsj         = MX.sym('Qsj',nq.all,d);
Qskj        = [Qsk Qsj];
Qdotsk      = MX.sym('Qdotsk',nq.all);
Qdotsj      = MX.sym('Qdotsj',nq.all,d);
Qdotskj     = [Qdotsk Qdotsj];
if nq.torqAct > 0
    a_ak        = MX.sym('a_ak',nq.torqAct);
    a_aj        = MX.sym('a_akmesh',nq.torqAct,d);
    a_akj       = [a_ak a_aj];
end
% Define CasADi variables for controls
vAk     = MX.sym('vAk',NMuscle);
if nq.torqAct > 0
    e_ak    = MX.sym('e_ak',nq.torqAct);
end

% Define CasADi variables for "slack" controls
dFTtildej   = MX.sym('dFTtildej',NMuscle,d);
Aj          = MX.sym('Aj',nq.all,d);
J           = 0; % Initialize cost function
eq_constr   = {}; % Initialize equality constraint vector
ineq_constr1   = {}; % Initialize inequality constraint vector
ineq_constr2   = {}; % Initialize inequality constraint vector
ineq_constr3   = {}; % Initialize inequality constraint vector
ineq_constr4   = {}; % Initialize inequality constraint vector
ineq_constr5   = {}; % Initialize inequality constraint vector
ineq_constr6   = {}; % Initialize inequality constraint vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step
h = tfk/N;
% Loop over collocation points
for j=1:d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unscale variables
    Qskj_nsc = Qskj.*(scaling.Qs'*ones(1,size(Qskj,2)));
    Qdotskj_nsc = Qdotskj.*(scaling.Qdots'*ones(1,size(Qdotskj,2)));
    FTtildekj_nsc = FTtildekj.*(scaling.FTtilde'*ones(1,size(FTtildekj,2)));
    dFTtildej_nsc = dFTtildej.*scaling.dFTtilde;
    Aj_nsc = Aj.*(scaling.Qdotdots'*ones(1,size(Aj,2)));
    vAk_nsc = vAk.*scaling.vA;
    
    QsQdotskj_nsc = MX(nq.all*2, d+1);
    QsQdotskj_nsc(1:2:end,:) = Qskj_nsc;
    QsQdotskj_nsc(2:2:end,:) = Qdotskj_nsc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get muscle-tendon lengths, velocities, and moment arms
    qinj    = Qskj_nsc(:, j+1);
    qdotinj = Qdotskj_nsc(:, j+1);
    [lMTj,vMTj,MAj] =  f_casadi.lMT_vMT_dM(qinj',qdotinj');
    % Derive the moment arms of all the muscles crossing each joint
    for i=1:nq.musAct
        MA_j.(model_info.ExtFunIO.coord_names.muscleActuated{i}) = ...
            MAj(mai(model_info.ExtFunIO.jointi.muscleActuated(i)).mus',...
            model_info.ExtFunIO.jointi.muscleActuated(i));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffj,FTj,Fcej,Fpassj,Fisoj] = ...
        f_casadi.forceEquilibrium_FtildeState_all_tendon(akj(:,j+1),...
        FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
        lMTj,vMTj,tensions);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get metabolic energy rate if in the cost function
    % Get muscle fiber lengths
    [~,lMtildej] = f_casadi.FiberLength_TendonForce_tendon(...
        FTtildekj_nsc(:,j+1),lMTj);
    % Get muscle fiber velocities
    [vMj,~] = f_casadi.FiberVelocity_TendonForce_tendon(...
        FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
        lMTj,vMTj);
    if strcmp(S.metabolicE.model,'Bhargava2004')
        % Get metabolic energy rate Bhargava et al. (2004)
        [e_totj,~,~,~,~,~] = f_casadi.getMetabolicEnergySmooth2004all(...
            akj(:,j+1),akj(:,j+1),lMtildej,vMj,Fcej,Fpassj,...
            MuscleMass',pctsts,Fisoj,model_info.mass,S.metabolicE.tanh_b);
    else
        error('No energy model selected');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get passive joint torques for dynamics
    Tau_passj = f_casadi.AllPassiveTorques(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
    % Get passive joint torques for cost function
    Tau_passj_cost = f_casadi.AllPassiveTorques_cost(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
    
    % Expression for the state derivatives at the collocation points
    Qsp_nsc      = Qskj_nsc*C(:,j+1);
    Qdotsp_nsc   = Qdotskj_nsc*C(:,j+1);
    FTtildep_nsc = FTtildekj_nsc*C(:,j+1);
    ap           = akj*C(:,j+1);
    if nq.torqAct > 0
        a_ap         = a_akj*C(:,j+1);
    end
    % Append collocation equations
    % Dynamic constraints are scaled using the same scale
    % factors as the ones used to scale the states
    % Activation dynamics (implicit formulation)
    eq_constr{end+1} = (h*vAk_nsc - ap)./scaling.a;
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = (h*dFTtildej_nsc(:,j) - FTtildep_nsc)./scaling.FTtilde';
    % Skeleton dynamics (implicit formulation)
    qdotj_nsc = Qdotskj_nsc(:,j+1); % velocity
    eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.Qs';
    eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./scaling.Qdots';
    % Torque actuator activation dynamics (explicit formulation)
    if nq.torqAct > 0
        da_adtj = f_casadi.ActuatorActivationDynamics(e_ak,a_akj(:,j+1));
        eq_constr{end+1} = (h*da_adtj - a_ap)./scaling.a_a;
    end

    % Add contribution to the cost function
    J = J + ...
        W.E          * B(j+1) *(f_casadi.J_muscles_exp(e_totj,W.E_exp))/model_info.mass*h + ...
        W.a          * B(j+1) *(f_casadi.J_muscles(akj(:,j+1)'))*h + ...
        W.q_dotdot   * B(j+1) *(f_casadi.J_not_arms_dof(Aj(model_info.ExtFunIO.jointi.noarmsi,j)))*h + ...
        W.pass_torq  * B(j+1) *(f_casadi.J_lim_torq(Tau_passj_cost))*h + ...
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(vAk))*h + ...
        W.slack_ctrl * B(j+1) *(f_casadi.J_muscles(dFTtildej(:,j)))*h;

    if nq.torqAct > 0
        J = J + W.e_arm      * B(j+1) *(f_casadi.J_torq_act(e_ak))*h;
    end
    if nq.arms > 0
        J = J + W.slack_ctrl * B(j+1) *(f_casadi.J_arms_dof(Aj(model_info.ExtFunIO.jointi.armsi,j)))*h;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call external function (run inverse dynamics)
    [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add path constraints
    for i=1:nq.all
        % total coordinate torque
        Ti = 0;

        % muscle moment
        cond_special_mtp = strcmp(S.subject.mtp_type,'2022paper') &&...
            contains(model_info.ExtFunIO.coord_names.all{i},'mtp');

        if ismember(i,model_info.ExtFunIO.jointi.muscleActuated) && ~cond_special_mtp
            % muscle forces
            FTj_coord_i = FTj(mai(i).mus',1);
            % total muscle moment
            M_mus_i = f_casadi.(['musc_cross_' num2str(sumCross(i))])...
            (MA_j.(model_info.ExtFunIO.coord_names.all{i}),FTj_coord_i);
            % add to total moment
            Ti = Ti + M_mus_i;
        end
        
        % torque actuator
        if nq.torqAct > 0 && ismember(i,model_info.ExtFunIO.jointi.torqueActuated)
            idx_act_i = find(model_info.ExtFunIO.jointi.torqueActuated(:)==i);
            T_act_i = a_akj(idx_act_i,j+1).*scaling.ActuatorTorque(idx_act_i);
            Ti = Ti + T_act_i;
        end

        % passive moment
        if ~ismember(i,model_info.ExtFunIO.jointi.floating_base)
            Ti = Ti + Tau_passj(i);
        end
        
        % total coordinate torque equals inverse dynamics torque
        eq_constr{end+1} = Tj(i,1) - Ti;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Activation dynamics (implicit formulation)
    act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
    act2 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tact);
    ineq_constr1{end+1} = act1;
    ineq_constr2{end+1} = act2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = Hilldiffj;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constraints to prevent parts of the skeleton to penetrate each other.
    if isfield(model_info.ExtFunIO,'origin')
        % Origins calcaneus (transv plane) at minimum 9 cm from each other.
        if isfield(model_info.ExtFunIO.origin,'calcn_r') && isfield(model_info.ExtFunIO.origin,'calcn_l') && ...
                ~isempty(model_info.ExtFunIO.origin.calcn_r) &&  ~isempty(model_info.ExtFunIO.origin.calcn_l)
            Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.calcn_r([1 3]),1) - ...
                Tj(model_info.ExtFunIO.origin.calcn_l([1 3]),1));
            ineq_constr3{end+1} = Qconstr;
        end
        % Constraint to prevent the arms to penetrate the skeleton
        % Origins femurs and ipsilateral hands (transv plane) at minimum
        % 18 cm from each other.
        if isfield(model_info.ExtFunIO.origin,'femur_r') && isfield(model_info.ExtFunIO.origin,'hand_r') && ...
                ~isempty(model_info.ExtFunIO.origin.femur_r) &&  ~isempty(model_info.ExtFunIO.origin.hand_r)
            Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.femur_r([1 3]),1) - ...
                Tj(model_info.ExtFunIO.origin.hand_r([1 3]),1));
            ineq_constr4{end+1} = Qconstr;
        end
        if isfield(model_info.ExtFunIO.origin,'femur_l') && isfield(model_info.ExtFunIO.origin,'hand_l') && ...
                ~isempty(model_info.ExtFunIO.origin.femur_l) &&  ~isempty(model_info.ExtFunIO.origin.hand_l)
            Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.femur_l([1 3]),1) - ...
                Tj(model_info.ExtFunIO.origin.hand_l([1 3]),1));
            ineq_constr4{end+1} = Qconstr;
        end
        % Origins tibia (transv plane) at minimum 11 cm from each other.
        if isfield(model_info.ExtFunIO.origin,'tibia_r') && isfield(model_info.ExtFunIO.origin,'tibia_l') && ...
                ~isempty(model_info.ExtFunIO.origin.tibia_r) &&  ~isempty(model_info.ExtFunIO.origin.tibia_l)
            Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.tibia_r([1 3]),1) - ...
                Tj(model_info.ExtFunIO.origin.tibia_l([1 3]),1));
            ineq_constr5{end+1} = Qconstr;
        end
        % Origins toes (transv plane) at minimum 10 cm from each other.
        if isfield(model_info.ExtFunIO.origin,'toes_r') && isfield(model_info.ExtFunIO.origin,'toes_l') && ...
                ~isempty(model_info.ExtFunIO.origin.toes_r) &&  ~isempty(model_info.ExtFunIO.origin.toes_l)
            Qconstr = f_casadi.J_nn_2(Tj(model_info.ExtFunIO.origin.toes_r([1 3]),1) - ...
                Tj(model_info.ExtFunIO.origin.toes_l([1 3]),1));
            ineq_constr6{end+1} = Qconstr;
        end
    end

end % End loop over collocation points

eq_constr = vertcat(eq_constr{:});
ineq_constr1 = vertcat(ineq_constr1{:});
ineq_constr2 = vertcat(ineq_constr2{:});
ineq_constr3 = vertcat(ineq_constr3{:});
ineq_constr4 = vertcat(ineq_constr4{:});
ineq_constr5 = vertcat(ineq_constr5{:});
ineq_constr6 = vertcat(ineq_constr6{:});


% Casadi function to get constraints and objective
coll_input_vars_def = {tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,Qdotsj,vAk,dFTtildej,Aj};
if nq.torqAct > 0
    coll_input_vars_def = [coll_input_vars_def,{a_ak,a_aj,e_ak}];
end
f_coll = Function('f_coll',coll_input_vars_def,...
    {eq_constr,ineq_constr1,ineq_constr2,ineq_constr3,...
    ineq_constr4,ineq_constr5,ineq_constr6,J});