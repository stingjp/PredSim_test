
close all
clear
clc

import casadi.*

% add repository to workspace
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/PreProcessing'])

% old dll
name_1 = 'PredSim_2D';
F1  = external('F',fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '.dll'])); 
load(fullfile(pathRepo,'Subjects',name_1,['F_' name_1 '_IO.mat'])); 
IO1 = IO;

% new osim
name_2 = 'subject1a_2D';
F2  = external('F',fullfile(pathRepo,'Subjects',name_2,['F_' name_2 '.dll'])); 
load(fullfile(pathRepo,'Subjects',name_2,['F_' name_2 '_IO.mat'])); 
IO2 = IO;

% new osim
name_3 = 'subject1a_2D_v2';
F3  = external('F',fullfile(pathRepo,'Subjects',name_3,['F_' name_3 '.dll'])); 
load(fullfile(pathRepo,'Subjects',name_3,['F_' name_3 '_IO.mat'])); 
IO3 = IO;

clearvars('IO')

Fs = {F1,F2,F3};
IOs = {IO1,IO2,IO3};
% legnames = {name_1,name_2,name_3};
legnames = {'old .dll','new dll','new dll with 3 spheres'};
markers = {'*','o','.'};

% Fs = {F1,F2};
% IOs = {IO1,IO2};
% legnames = {name_1,name_2};
% markers = {'.','o'};

figure(1)

figure(2)

for ji=1:length(Fs)

    F = Fs{ji};
    IO = IOs{ji};
    
    coord_names = fieldnames(IO.coordi);
    n_coord = length(coord_names);
    
    Qs = zeros(n_coord,1);
    Qdots = Qs;
    Qddots = Qs;
    
    Qs(IO.coordi.pelvis_ty) = 0.95;
    
    T1 = calcID(Qs,Qdots,Qddots,F);
    
    x = linspace(-10,10,40);
    
    Tall = zeros(length(T1),length(x));
    
    for i=1:length(x)
        Qs(IO.coordi.pelvis_tx) = x(i);
        T2 = calcID(Qs,Qdots,Qddots,F);
        Tall(:,i) = T2;
    end
    
    % diff2 = T2-T1;
    % disp(diff2(1:n_coord))
    
    figure(1)

    subplot(3,2,2)
    hold on
    p1=plot(x,Tall(IO.GRFs.right_foot(1),:),markers{ji});
    ylabel('right GRF x')
    xlabel('pelvis tx (m)')
    title('right')

    subplot(3,2,4)
    hold on
    plot(x,Tall(IO.GRFs.right_foot(2),:),markers{ji},'Color',p1.Color)
    ylabel('right GRF y')
    xlabel('pelvis tx (m)')

    subplot(3,2,6)
    if IO.GRFs.right_foot(3) ~= IO.GRFs.right_foot(1)
        hold on
        plot(x,Tall(IO.GRFs.right_foot(3),:),markers{ji},'Color',p1.Color)
        ylabel('right GRF z')
        xlabel('pelvis tx (m)')
    end

    subplot(3,2,1)
    hold on
    plot(x,Tall(IO.GRFs.left_foot(1),:),markers{ji},'Color',p1.Color)
    ylabel('left GRF x')
    xlabel('pelvis tx (m)')
    title('left')

    if ji==length(Fs)
        legend(legnames,'Interpreter','none','Location','best')
    end

    subplot(3,2,3)
    hold on
    plot(x,Tall(IO.GRFs.left_foot(2),:),markers{ji},'Color',p1.Color)
    ylabel('left GRF y')
    xlabel('pelvis tx (m)')

    subplot(3,2,5)
    if IO.GRFs.right_foot(3) ~= IO.GRFs.right_foot(1)
        hold on
        plot(x,Tall(IO.GRFs.left_foot(3),:),markers{ji},'Color',p1.Color)
        ylabel('left GRF z')
        xlabel('pelvis tx (m)')
    end
    

    if isfield(IO.GRFs,'contact_sphere_1')
        figure(2)
        for i=1:4
            idxsi = IO.GRFs.(['contact_sphere_' num2str(i)]);
            for j=1:3
                subplot(3,4,i+4*(j-1))
                hold on
                plot(x,Tall(idxsi(j),:),markers{ji},'Color',p1.Color)
                if j==1
                    title(['contact_sphere_' num2str(i)],'Interpreter','none')
                end
            end
            xlabel('pelvis tx (m)')
        end
        subplot(3,4,1)
        ylabel('GRF x')
        subplot(3,4,5)
        ylabel('GRF y')
        subplot(3,4,9)
        ylabel('GRF z')

    end

end

%%
function [Ts] = calcID(qs,dqs,ddqs,Fext)
import casadi.*
qsdqs = zeros(length(qs)*2,1);
qsdqs(1:2:end) = qs;
qsdqs(2:2:end) = dqs;

res = Fext([qsdqs; ddqs]);
Ts = full(res);

end