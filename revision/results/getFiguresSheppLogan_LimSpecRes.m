
clc; clearvars; close all;

% figure directory
figDir = [pwd '/figures/exampleDisk_LimSpecRes/'];
if ~exist(figDir,'dir'), mkdir(figDir); end

%% load data

load('example_SheppLogan_nonInversCrime__LimSpecResFinal.mat');

%% save results

% GT
for i=1:size(strOut.GT.A,2)
    imwrite(rescale(reshape(strOut.GT.A(:,i),[512 512])),[figDir 'Disk_LimSpecRes_GT_' num2str(i) '.png']);
end

% RU
for i=1:size(strOut.RU.A,2)
    imwrite(rescale(reshape(strOut.RU.A(:,i),[512 512])),[figDir 'Disk_LimSpecRes_RU_' num2str(i) '.png']);
end

% UR
for i=1:size(strOut.UR.A,2)
    imwrite(rescale(reshape(strOut.UR.A(:,i),[512 512])),[figDir 'Disk_LimSpecRes_UR_' num2str(i) '.png']);
end

% cJoint
for i=1:size(strOut.cJoint.A,2)
    imwrite(rescale(reshape(strOut.cJoint.A(:,i),[512 512])),[figDir 'Disk_LimSpecRes_cJoint_' num2str(i) '.png']);
end

% ADJUST
for i=1:size(strOut.ADJUST.A,2)
    imwrite(rescale(reshape(strOut.ADJUST.A(:,i),[512 512])),[figDir 'Disk_LimSpecRes_ADJUST_' num2str(i) '.png']);
end
