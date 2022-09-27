
%clc; clearvars; close all;

% figure directory
figDir = [pwd '/figures/exampleMixedDisk/'];
if ~exist(figDir,'dir'), mkdir(figDir); end

%% load data

load('example_MixedDisk_nonInversCrime_fromArsenicOnward.mat');

%% save results

% GT
for i=1:size(strOut.RU.A,2)
    imwrite(rescale(reshape(strOut.GT.A(:,i),[512 512])),[figDir 'MixedDisk_GT_' num2str(i) '.png']);
end

% RU
for i=1:size(strOut.RU.A,2)
    imwrite(rescale(reshape(strOut.RU.A(:,i),[512 512])),[figDir 'MixedDisk_RU_' num2str(i) '.png']);
end

% UR
for i=1:size(strOut.UR.A,2)
    imwrite(rescale(reshape(strOut.UR.A(:,i),[512 512])),[figDir 'MixedDisk_UR_' num2str(i) '.png']);
end

% cJoint
for i=1:size(strOut.cJoint.A,2)
    imwrite(rescale(reshape(strOut.cJoint.A(:,i),[512 512])),[figDir 'MixedDisk_cJoint_' num2str(i) '.png']);
end

% ADJUST
for i=1:size(strOut.ADJUST.A,2)
    imwrite(rescale(reshape(strOut.ADJUST.A(:,i),[512 512])),[figDir 'MixedDisk_ADJUST_' num2str(i) '.png']);
end

%%
xV = linspace(20,80,100); plot(xV,strOut.GT.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf,[figDir 'MixedDisk_GT_f'],'epsc');
xV = linspace(20,80,100); plot(xV,strOut.RU.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf,[figDir 'MixedDisk_RU_f'],'epsc');
xV = linspace(20,80,100); plot(xV,strOut.UR.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf,[figDir 'MixedDisk_UR_f'],'epsc');
xV = linspace(20,80,100); plot(xV,strOut.cJoint.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf,[figDir 'MixedDisk_cJoint_f'],'epsc');
xV = linspace(20,80,100); plot(xV,strOut.ADJUST.F','LineWidth',2); pbaspect([3 1 1]); ylabel('Attenuation'); xlabel('Energy (keV)'); set(gca,'FontSize',12)
saveas(gcf,[figDir 'MixedDisk_ADJUST_f'],'epsc');