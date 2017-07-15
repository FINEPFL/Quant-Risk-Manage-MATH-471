% QRM Assignment10
%   Authors:
%            Mengjie Zhao
%            Tianxiao Ma
%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

% params
N   = 1000;  % N bonds
M   = 10000; % MC simulation samples
PI  = 0.01;  % PI is the probability a single bond default
RHO = 0.01;

% Q1
rho_correlated  = repmat(RHO, N, N);
rho_correlated(1:N+1:end) = 1; % let diag be 1
 
% Gaussian
X_corr = copularnd('Gaussian', rho_correlated , M);
for i = 1:M
L_g_corr(i) = length(find(X_corr(i,:)<=PI));
end
Percentile_g = quantile(L_g_corr,[0.95,0.99,0.999]);
corr_g = corrcoef(double(X_corr(:,1)<=0.01),double(X_corr(:,2)<=0.01));
corr_g(2);

% subplot(1, 3, 1)
% histogram(L_g_corr,'Normalization','pdf'); hold on;
% plot(binopdf((0:50), N, PI), 'linewidth', 3);
% axis([0 50 0 0.14])
% legend('Dependent', 'Independent')
% xlabel('Defaults','interpreter','latex'); 
% ylabel('Probability','interpreter','latex')
% title('Gaussian Copula','interpreter','latex')
% set(gca, 'fontsize', 15)

% % t student 3
X_corr_t3 = copularnd('t', rho_correlated, 3, M);
for i = 1:M
L_t3_corr(i) = length(find(X_corr_t3(i,:)<=PI));
end
Percentile_t3 = quantile(L_t3_corr,[0.95,0.99,0.999]);
corr_t3 = corrcoef(double(X_corr_t3(:,1)<=0.01),double(X_corr_t3(:,2)<=0.01));
corr_t3(2);

% subplot(1, 3, 2)
% histogram(L_t3_corr, 420, 'Normalization','pdf'); hold on;
% plot(binopdf((0:50), N, PI), 'linewidth', 3);
% axis([0 50 0 0.7])
% legend('Dependent', 'Independent')
% xlabel('Defaults','interpreter','latex'); 
% ylabel('Probability','interpreter','latex')
% title('Student-t Copula, 3','interpreter','latex')
% set(gca, 'fontsize', 15)

% % t student 10
X_corr_t10 = copularnd('t', rho_correlated, 10, M);
for i = 1:M
L_t10_corr(i) = length(find(X_corr_t10(i,:)<=PI));
end
Percentile_t10 = quantile(L_t10_corr,[0.95,0.99,0.999]);
corr_t10 = corrcoef(double(X_corr_t10(:,1)<=0.01),double(X_corr_t10(:,2)<=0.01));
corr_t10(2);

% subplot(1, 3, 3)
% histogram(L_t10_corr, 120, 'Normalization','pdf'); hold on;
% plot(binopdf((0:50), N, PI), 'linewidth', 3);
% axis([0 50 0 0.3])
% legend('Dependent', 'Independent')
% xlabel('Defaults','interpreter','latex');
% ylabel('Probability','interpreter','latex')
% title('Student-t Copula, 10','interpreter','latex')
% set(gca, 'fontsize', 15)

% Q2
% standard MC
clearvars
var_mc = atanh(0.99 * 2 - 1);
for i = 1:1000
    sample   = atanh(rand(1000, 1) .* 2 - 1);
    ES_mc(i) = mean(sample(sample>var_mc));
end
mean(ES_mc);
sqrt(var(ES_mc)/length(ES_mc));

% importance sampling
x  = atanh(0.99 * 2 - 1);
x_ = x - log(1-rand(1000, 1000));
for i = 1:1000
    sample = x_(:, i);
    sample = sample(sample>x);
    d = 2 * ((exp(sample)+exp(-sample)).^2 .* exp(-sample+x)).^(-1);
    ES_is(i) = mean(sample.*d/(1-0.99));
end
mean(ES_is);
sqrt(var(ES_is)/length(ES_is));
    