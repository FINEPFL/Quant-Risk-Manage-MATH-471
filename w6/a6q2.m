% QRM Assignment 6         %
% Authors:                 %
%           Mengjie Zhao   %
%           Tianxiao Ma    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc

%% Processing the data
% load full data
raw_data = load('msft.txt');

% previous two years and last three years
prep_price = raw_data(1:501, 3);
last_price = raw_data(502:end, 3);

%% 1. Fit GARCH(1, 1) to previous two years of returns
% (a) compute log returns Rt, we can use the price2ret function
mu = mean(prep_price);
Rt = price2ret(prep_price);
Xt = Rt - mean(Rt);
n  = length(Xt);

% (b) let sigma0 = std(Xt)
sigma0 = std(Xt);

% (c) use mle to estimate parameters
[params val_obj] = ...
    fmincon(...
            @(vec)neg_likelihood(vec, Xt, sigma0, n),...
            [0.01, 0.01, 0.01],...
            [],[],[],[],[0,0,0],[Inf,Inf,Inf]...
           );
       % [0.01, 0.01, 0.01] is the initial point of gradient descent
       % params = [alpha0, alpha1, beta1]
       % val_obj = residual of estimation - smaller than tolerance 
       
% 2. Use last three years for testing
% Gaussian innovation
mu_last = mean(last_price);
Rt_last = price2ret(last_price);
Xt_last = Rt_last - mean(Rt_last);

% calculating sigma, slide 15
sigma = zeros(length(last_price), 1);
sigma(1) = sigma0;
for t = 2:length(last_price)
    sigma(t) = sqrt(params(1) + params(2)*Xt_last(t-1)^2 + params(3)*sigma(t-1)^2);
end
       
alpha  = [0.95, 0.99];
VaR_95 = zeros(length(last_price), 1);
VaR_99 = zeros(length(last_price), 1);

for i=1:length(last_price)
    VaR_95(i) = last_price(i)*sigma(i)*norminv(0.95, 0, 1);
    VaR_99(i) = last_price(i)*sigma(i)*norminv(0.99, 0, 1);
end

plot(VaR_95); hold on; plot(VaR_99);
legend('Gaussian, VaR_{95}', 'Gaussian, VaR_{99}')

% 3
% ************* TO ADD IF THE RETURN ON NEXT DAY BREACHES **************

% 4. Now innovations have t-student distribution, degree of freedom v added
[params_t obj_val_t] = ...
        fmincon(...
                @(vec)neg_likelihood_t(vec,  Xt, sigma0, n),...
                [0.01, 0.01, 0.01, 3],... % note the initial degree of freedom is set to 3
                [],[],[],[],[0,0,0,0,0],[Inf,Inf,Inf,Inf]...
               );
               % params_t = [alpha0, alpha1, beta1, v]

% calculating sigma assuming t distribution, slide 15
sigma_T = zeros(length(last_price), 1);
sigma_T(1) = sigma0;
for t = 2:length(last_price)
    sigma_T(t) = sqrt(params_t(1) + params_t(2)*Xt_last(t-1)^2 + params_t(3)*sigma_T(t-1)^2);
end
       
alpha  = [0.95, 0.99];
VaR_95_T = zeros(length(last_price), 1);
VaR_99_T = zeros(length(last_price), 1);
reg = sqrt( (params_t(4)-2) / params_t(4) ); % sqrt((v-2)/v) where v is fitted degree of freedom
for i=1:length(last_price)
    VaR_95_T(i) = last_price(i)*sigma_T(i)*reg*tinv(0.95, params_t(4));
    VaR_99_T(i) = last_price(i)*sigma_T(i)*reg*tinv(0.99, params_t(4));
end
figure
plot(VaR_95_T); hold on; plot(VaR_99_T);
legend('t-Student, VaR_{95}', 't-Student, VaR_{99}')

% 5. Compare calculated and theoretical variance of Xt
% (a) for gaussian assumption
% slide 21
testing_g = params(2) + params(3)
theo_var_g = params(1)/(1 - params(2) - params(3))

% (b) for t student
testing_t = params_t(2) + params_t(3)
theo_var_t = params_t(1)/(1 - params_t(2) - params_t(3))








