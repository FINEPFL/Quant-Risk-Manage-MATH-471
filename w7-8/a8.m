% Quantitative Risk Management
%   Assignment 7
%   Authors:
%               Mengjie Zhao
%               Tianxiao Ma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc

raw_data = xlsread('data.xls');
x1 = raw_data(:, 1);
x2 = raw_data(:, 2);
N  = size(raw_data, 1);

% scatterhist of real data
figure
scatterhist(x1, x2);
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')

% generate pseudo data
x1_pseu = zeros(size(x1));
x2_pseu = zeros(size(x2));

for i = 1:N
    x1_pseu(i) = 1/(N+1) * length(find(x1<=x1(i)));
    x2_pseu(i) = 1/(N+1) * length(find(x2<=x2(i)));
end

PSE = [x1_pseu x2_pseu];

figure
scatterhist(x1_pseu, x2_pseu);
xlabel('$x_1^{pseu}$', 'interpreter', 'latex')
ylabel('$x_2^{pseu}$', 'interpreter', 'latex')

% MLE 
% theta_x is the parameter estimated and vc is the negative loglikelihood
[theta_c, vc] = fmincon(@(x) mll(PSE, x, 'Clayton'), 0, 0, 0);
[theta_g, vg] = fmincon(@(x) mll(PSE, x, 'Gumbel'), 1, 0, 0);
[theta_f, vf] = fmincon(@(x) mll(PSE, x, 'Frank'), 0, 0, 0);

