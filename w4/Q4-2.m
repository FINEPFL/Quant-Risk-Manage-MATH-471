%% Homework 4-2 - Quantitative Risk Management
% By Tianxiao MA
% Date: March 20th, 2017
clear;
clc;
load('Workspace.mat');
%%
X = [IBM,MCD,MMM,WMT];
F = [SP500];

X = price2ret(X); % Continiously compounded by default, not need to specify
F = price2ret(F);
n = length(X);
F = [ones(n,1),F];

B = inv(F'*F)*F'*X;

Epsilon = X - F*B;
SampleCorrelationEpsilon = corrcoef(Epsilon);
SampleCorrelationX = corrcoef(X);