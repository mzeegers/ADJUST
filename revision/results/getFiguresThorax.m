
clc; clearvars; close all;

% figure directory
figDir = [pwd '/figures/exampleThorax/'];
if ~exist(figDir,'dir'), mkdir(figDir); end

%% load data

load('example_Thorax_nonInversCrimeFinal.mat');

%% save results

% GT
for i=1:3
    if i<3
        imwrite(rescale(reshape(strOut.strProb.A(:,i),[512 512])),[figDir 'Thorax_GT_' num2str(i) '.png']);
    else
        imwrite(rescale(reshape(strOut.strProb.A(:,3)+strOut.strProb.A(:,4)+strOut.strProb.A(:,5),[512 512])),[figDir 'Thorax_GT_' num2str(i) '.png']);
    end
end

% RU
for i=1:size(strOut.RU.A,2)
    imwrite(rescale(reshape(strOut.RU.A(:,i),[512 512])),[figDir 'Thorax_RU_' num2str(i) '.png']);
end

% UR
for i=1:size(strOut.UR.A,2)
    imwrite(rescale(reshape(strOut.UR.A(:,i),[512 512])),[figDir 'Thorax_UR_' num2str(i) '.png']);
end

% cJoint
for i=1:size(strOut.cJoint.A,2)
    imwrite(rescale(reshape(strOut.cJoint.A(:,i),[512 512])),[figDir 'Thorax_cJoint_' num2str(i) '.png']);
end

% ADJUST
for i=1:size(strOut.ADJUST.A,2)
    imwrite(rescale(reshape(strOut.ADJUST.A(:,i),[512 512])),[figDir 'Thorax_ADJUST_' num2str(i) '.png']);
end

%%

F     = readmatrix('Thorax_MatrixFXray100chan_17matBoneBloodIodineSofttissueBloodLungsIodine_20to80kV_CORRECTED.csv');
F  = F(1:3,:);

xV = linspace(20,80,100); plot(xV,F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf, [figDir 'Thorax_GT_F'],'epsc');
saveas(gcf, [figDir 'Thorax_GT_F'],'png');

% xV = linspace(20,80,100); plot(xV,strOut.strProb.F','LineWidth',2); pbaspect([3 1 1]); ylabel('attenuation'); xlabel('frequency (keV)'); set(gca,'FontSize',12)
% saveas(gcf, [figDir 'Thorax_GT_F'],'epsc');
% saveas(gcf, [figDir 'Thorax_GT_F'],'png');

xV = linspace(20,80,100); plot(xV,strOut.RU.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf, [figDir 'Thorax_RU_F'],'epsc');
saveas(gcf, [figDir 'Thorax_RU_F'],'png');

xV = linspace(20,80,100); plot(xV,strOut.UR.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf, [figDir 'Thorax_UR_F'],'epsc');
saveas(gcf, [figDir 'Thorax_UR_F'],'png');

xV = linspace(20,80,100); plot(xV,strOut.cJoint.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf, [figDir 'Thorax_cJoint_F'],'epsc');
saveas(gcf, [figDir 'Thorax_cJoint_F'],'png');

xV = linspace(20,80,100); plot(xV,strOut.ADJUST.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf, [figDir 'Thorax_ADJUST_F'],'epsc');
saveas(gcf, [figDir 'Thorax_ADJUST_F'],'png');