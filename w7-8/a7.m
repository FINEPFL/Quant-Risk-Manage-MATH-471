% Quantitative Risk Management
%   Assignment 7
%   Authors:
%               Mengjie Zhao
%               Tianxiao Ma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc

% Q1 
sigma = 0:0.1:10;
rho_max = (exp(sigma) - 1)  ./ (sqrt((exp(1)-1) .* (exp(sigma.^2) - 1)));
rho_min = (exp(-sigma) - 1) ./ (sqrt((exp(1)-1) .* (exp(sigma.^2) - 1)));
plot(sigma, rho_max, 'r', 'linewidth', 1.5); hold on; plot(sigma, rho_min, 'b-o', 'linewidth', 1.5, 'markersize', 4);
xlabel('$\sigma$', 'interpreter', 'latex')
ylabel('$\rho$', 'interpreter', 'latex')
legend('\rho_{max}', '\rho_{min}');
axis([0 5 -0.91 1.01]); grid on;
set(gca, 'fontsize', 15)
