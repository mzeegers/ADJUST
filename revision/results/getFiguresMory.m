
%clc; clearvars; close all;

% figure directory
figDir = [pwd '/figures/exampleMory/'];
if ~exist(figDir,'dir'), mkdir(figDir); end

%% load data

load('example_ComparisonsMory_nonICFinal.mat');

%% save results

% GT
for i=1:size(str.True.A,2)
    imwrite(rescale(reshape(str.True.A(:,i),[128 128])),[figDir 'Mory_GT_' num2str(i) '.png']);
end

% RU
for i=1:size(str.RU.A,2)
    imwrite(rescale(reshape(str.RU.A(:,i),[128 128])),[figDir 'Mory_RU_' num2str(i) '.png']);
end

% UR
for i=1:size(str.UR.A,2)
    imwrite(rescale(reshape(str.UR.A(:,i),[128 128])),[figDir 'Mory_UR_' num2str(i) '.png']);
end

% cJoint
for i=1:size(str.Joint.A,2)
    imwrite(rescale(reshape(str.Joint.A(:,i),[128 128])),[figDir 'Mory_cJoint_' num2str(i) '.png']);
end

% ADJUST
for i=1:size(str.ADJUST.A,2)
    imwrite(rescale(reshape(str.ADJUST.A(:,i),[128 128])),[figDir 'Mory_ADJUST_' num2str(i) '.png']);
end

% Cai
for i=1:size(str.Cai.A,2)
    imwrite(rescale(reshape(str.Cai.A(:,i),[128 128])),[figDir 'Mory_Cai_' num2str(i) '.png']);
end

% Weidinger
for i=1:size(str.Weidinger.A,2)
    imwrite(rescale(reshape(str.Weidinger.A(:,i),[128 128])),[figDir 'Mory_Weidinger_' num2str(i) '.png']);
end

% Long
for i=1:size(str.Long.A,2)
    imwrite(rescale(reshape(str.Long.A(:,i),[128 128])),[figDir 'Mory_Long_' num2str(i) '.png']);
end

% Mechlem
for i=1:size(str.Mechlem.A,2)
    imwrite(rescale(reshape(str.Mechlem.A(:,i),[128 128])),[figDir 'Mory_Mechlem_' num2str(i) '.png']);
end

% Barber
for i=1:size(str.Barber.A,2)
    imwrite(rescale(reshape(str.Barber.A(:,i),[128 128])),[figDir 'Mory_Barber_' num2str(i) '.png']);
end