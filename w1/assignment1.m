% This is a script for solving assignment 1 of Quantitative Risk Management 
% course at EPFL 2017. If you meet any problem when evalutating or running
% the code, please feel free to drop an email to the authors.
% *************************************************************************
% Author:
%         Mengjie Zhao, mengjie.zhao@epfl.ch
%         Tianxiao Ma,  tianxiao.ma@epfl.ch
% *************************************************************************
% Here we acknowlege that all following contents are working results of
% ourselves, w.r.t. EPFL honor code.

clearvars; close all; clc

%% Question 2
% 1.
S_t = 100; r_t = 0.05; sigma_t = 0.2;

u_1_tpd   = 0;
std_1_tpd = 0.01;

u_2_tpd   = 0;
std_2_tpd = 0.0001;

u_3_tpd   = 0;
std_3_tpd = 0.001;

corr_1_3 = -0.5;

% generating data for calculating loss -- L and L_delta
N = 10000;
MU = [u_1_tpd; u_2_tpd; u_3_tpd]';

% X2_tpd is independent with X1_tpd and X3_tpd, cor(X1_tpd, X3_tpd) = -0.5,
% then we build the correlation matrix, then map it to covariace matrix
% since it is required by the mvnrnd function
cor_mat = [   1   0  -0.5;
              0   1     0;
           -0.5   0     1];
var_map = [std_1_tpd; std_2_tpd; std_3_tpd] *...
                            [std_1_tpd; std_2_tpd; std_3_tpd]';
cov_mat = cor_mat .* var_map;
R = mvnrnd(MU, cov_mat, N);
corr(R(:, 1), R(:, 3)) % verify the correlation is -0.5
corr(R(:, 1), R(:, 2)) % verify the correlation is 0

% def parameters ànd helper functions for calculation
t = 0; T = 1; K = 100; Delta = 1/252;
get_d1 = @(t, T, K, S, r, sigma) (log(S./K) + (r + 0.5*sigma.^2) * (T-t))/...
                                                    (sigma * sqrt(T-t));
get_d2 = @(t, T, K, S, r, sigma) (get_d1(t, T, K, S, r, sigma) - ...
                                                     sigma * sqrt(T-t));
                                               
get_C_bs = @(t, T, K, S, r, sigma) (S*normcdf(get_d1(t, T, K, S, r, sigma)))...
                - exp(-r*(T-t)) * K * normcdf(get_d2(t, T, K, S, r, sigma));

d1_t = get_d1(t, T, K, S_t, r_t, sigma_t);
d2_t = get_d2(t, T, K, S_t, r_t, sigma_t);
C_bs_t = get_C_bs(t, T, K, S_t, r_t, sigma_t);

% in matrix R we have X1_tpd, X2_tpd and X3_tpd respectively. tpd means t
% plus delat
t = Delta;
S_tpd = exp(R(:, 1) + log(S_t));
r_tpd = R(:, 2) + r_t;
sigma_tpd = R(:, 3) + sigma_t;

for i=1:N
C_bs_tpd(i) = get_C_bs(t, T, K, S_tpd(i), r_tpd(i), sigma_tpd(i));
end
L_t_tpd = -(C_bs_tpd - C_bs_t);

% get handles from hist to make the hist to pdf and integral to 1
[N, X] = hist(L_t_tpd, 95);
bar(X, N/(sum(N)*diff(X(1:2))), 1)
grid on;
xlabel('$L(t, t+\Delta)$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('Probability')
set(gca, 'fontsize', 12)

delta = normcdf(d1_t)
rho = K * T * exp(-r_t * T) * normcdf(d2_t)
vega = S_t * normpdf(d1_t) * sqrt(T)

lized_L_t_tpd = -(delta * S_t * R(:, 1) + rho * R(:, 2) + vega * R(:, 3));
figure(2)
[N_, X_] = hist(lized_L_t_tpd, 95);
bar(X_, N_/(sum(N_)*diff(X_(1:2))), 1)
grid on
xlabel('$L^\delta(t, t+\Delta)$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('Probability')
set(gca, 'fontsize', 12)
