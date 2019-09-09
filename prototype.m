clc;
clear;
tic
 %% ================================Set Parameters=======================%%
 gridSpacing = 16;

 %% =====================================================================%%

 out=regexp(pwd,'\','split');
 setDir = '';
 for i=1:length(out)-1
     setDir = fullfile(setDir,out(i));
 end
 setDir = char(fullfile(setDir,'img'));
 imds = imageDatastore(setDir,'IncludeSubfolders',true,'LabelSource',...
     'foldernames');
 clear i out setDir
 %% =====================Calculate Sift Features=======================%%

 %createHistogram(imds.Files{17,1});
 desc = denseSIFT(imds.Files{21,1}, gridSpacing);
%load('defineOrient.mat');
toc;