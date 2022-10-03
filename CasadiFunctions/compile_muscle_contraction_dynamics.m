function [] = compile_muscle_contraction_dynamics(S)
% --------------------------------------------------------------------------
% compile_muscle_contraction_dynamics
%   This function compiles the muscle contraction dynamics functions into a
%   .dll.
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
% Original date: 27/Sept/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------




import casadi.*


%% Muscle contraction dynamics
% Function for Hill-equilibrium
a           = SX.sym('a',1); % Muscle activations
FTtilde     = SX.sym('FTtilde',1); % Normalized tendon forces
dFTtilde    = SX.sym('dFTtilde',1); % Time derivative tendon forces
lMT         = SX.sym('lMT',1); % Muscle-tendon lengths
vMT         = SX.sym('vMT',1); % Muscle-tendon velocities
FMo         = SX.sym('FMo',1); % optimal force
lMo         = SX.sym('lMo',1); % optimal fibre length
lTs         = SX.sym('lTs',1); % tendon slack langth
alphao      = SX.sym('alphao',1); % optimal pennation angle
vMmax       = SX.sym('vMmax',1); % max contraction velocity
tension     = SX.sym('tension',1); % Tensions
tendon_stiff = SX.sym('tendon_stiff',1); % Tendon stiffness
tendon_stiff_shift = SX.sym('tendon_stiff_shift',1); % shift curve
dampingCoefficient = SX.sym('d',1); % damping
stiffness_shift = SX.sym('shift',1); % shift passive curve
stiffness_scale = SX.sym('scale',1); % scale passive curve
muscle_strength = SX.sym('shift',1); % scale active curve

load('Fvparam.mat','Fvparam');
load('Fpparam.mat','Fpparam');
load('Faparam.mat','Faparam');

[err, FT, Fce, Fpass, Fiso, vMmax1, massM] = ...
    ForceEquilibrium_FtildeState_all_tendon(a,FTtilde,dFTtilde,lMT,vMT,FMo,lMo,...
    lTs,alphao,vMmax,Fvparam,Fpparam,Faparam,tension,tendon_stiff,tendon_stiff_shift,...
    0,dampingCoefficient,stiffness_shift,stiffness_scale,muscle_strength);


SX_in = {a,FTtilde,dFTtilde,lMT,vMT,FMo,lMo,lTs,alphao,vMmax,tension,tendon_stiff,...
    tendon_stiff_shift,dampingCoefficient,stiffness_shift,stiffness_scale,muscle_strength};
SX_out = {err, FT, Fce, Fpass, Fiso, vMmax1, massM};
for i=1:length(SX_in)
    name_in{i} = SX_in{i}.name;
end

f_ForceEquilibrium_FtildeState_all_tendon = Function('f_ForceEquilibrium_FtildeState_all_tendon',...
    SX_in,SX_out,name_in,{'Hilldiff','FT','Fce','Fpass','Fiso','vMmax','massM'});

%
[lM,lMtilde,lT] = FiberLength_TendonForce_tendon(FTtilde,lMo,lTs,alphao,...
    lMT,tendon_stiff,tendon_stiff_shift,0);

SX_in = {FTtilde,lMo,lTs,alphao,lMT,tendon_stiff,tendon_stiff_shift};
SX_out = {lM,lMtilde};
clearvars('name_in')
for i=1:length(SX_in)
    name_in{i} = SX_in{i}.name;
end

f_FiberLength_TendonForce_tendon = Function('f_FiberLength_TendonForce_tendon',...
    SX_in,SX_out,name_in,{'lM','lMtilde'});

%
[vM,vMtilde,vT] = FiberVelocity_TendonForce_tendon(FTtilde,dFTtilde,lMo,...
    lTs,alphao,vMmax,lMT,vMT,tendon_stiff,tendon_stiff_shift,0);

SX_in = {FTtilde,dFTtilde,lMo,lTs,alphao,vMmax,lMT,vMT,tendon_stiff,tendon_stiff_shift};
SX_out = {vM,vMtilde};
clearvars('name_in')
for i=1:length(SX_in)
    name_in{i} = SX_in{i}.name;
end

f_FiberVelocity_TendonForce_tendon = Function('f_FiberVelocity_TendonForce_tendon',...
    SX_in,SX_out,name_in,{'vM','vMtilde'});

SX_in = {FTtilde,dFTtilde,lMo,lTs,alphao,vMmax,lMT,vMT,tendon_stiff,tendon_stiff_shift};
SX_out = {lT,vT};
clearvars('name_in')
for i=1:length(SX_in)
    name_in{i} = SX_in{i}.name;
end

f_lT_vT = Function('f_lT_vT',SX_in,SX_out,name_in,{'lT','vT'});


%% Generate C code
tmp_dir = pwd;
cd(fullfile(S.misc.main_path,'CasadiFunctions'))

% create code generator
CG = CodeGenerator('muscle_contraction_dynamics.c');
% add functions add Jacobian
CG.add(f_ForceEquilibrium_FtildeState_all_tendon);
CG.add(f_FiberLength_TendonForce_tendon);
CG.add(f_FiberVelocity_TendonForce_tendon);
CG.add(f_lT_vT);

CG.add(f_ForceEquilibrium_FtildeState_all_tendon.jacobian());
CG.add(f_FiberLength_TendonForce_tendon.jacobian());
CG.add(f_FiberVelocity_TendonForce_tendon.jacobian());
CG.add(f_lT_vT.jacobian());

% generate C code
CG.generate([fullfile(S.misc.main_path,'CasadiFunctions') '/']);

%% Compile
options.flags = {'/O2'};
options.folder = pwd;
options.cleanup = false;
options.name = 'muscle_contraction_dynamics';
options.temp_suffix = false;
Ip = Importer('muscle_contraction_dynamics.c','shell',options);


cd(tmp_dir)




