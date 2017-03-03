clearvars; close all; clc

%% Question 1
S     = [100, 50, 25]';
share = [1, 3, 5]';
mean_ = [0, 0, 0]';
std_  = ([1, 2, 3] .* 0.001)';
alpha = 0.90:0.01:0.99;
nu    = [3, 10, 50]';
N     = 10000;
VaR   = zeros(length(alpha), 4);
VaR_m = zeros(length(alpha), 4);
ES    = zeros(length(alpha), 4);

% section 1 - (a)
X_    = trnd(nu(1), N, 3);
X     = X_ .* repmat((std_./sqrt(nu(1)/(nu(1)-2)))', N, 1);
L_3   = sum(repmat((S.*share)', N, 1) .* (1-exp(X)), 2);

VaR(:, 1)   = quantile(L_3, alpha);
VaR_m(:, 1) = VaR(:, 1) - mean(L_3);
for i=1:length(alpha)
    idx = L_3 > VaR(i, 1);
    ES(i, 1) = mean(L_3(idx));
end

% section 1 - (b)
X_    = trnd(nu(2), N, 3);
X     = X_ .* repmat((std_./sqrt(nu(2)/(nu(2)-2)))', N, 1);
L_10   = sum(repmat((S.*share)', N, 1) .* (1-exp(X)), 2);

VaR(:, 2)   = quantile(L_10, alpha);
VaR_m(:, 2) = VaR(:, 2) - mean(L_10);
for i=1:length(alpha)
    idx = L_10 > VaR(i, 2);
    ES(i, 2) = mean(L_10(idx));
end

% section 1 - (c)
X_    = trnd(nu(3), N, 3);
X     = X_ .* repmat((std_./sqrt(nu(3)/(nu(3)-2)))', N, 1);
L_50   = sum(repmat((S.*share)', N, 1) .* (1-exp(X)), 2);

VaR(:, 3)   = quantile(L_50, alpha);
VaR_m(:, 3) = VaR(:, 3) - mean(L_50);
for i=1:length(alpha)
    idx = L_50 > VaR(i, 3);
    ES(i, 3) = mean(L_50(idx));
end

% section 1 - (d)
X_    = randn(N, 3);
X     = X_ .* repmat(std_', N, 1);
L_N   = sum(repmat((S.*share)', N, 1) .* (1-exp(X)), 2);

VaR(:, 4)   = quantile(L_N, alpha);
VaR_m(:, 4) = VaR(:, 4) - mean(L_N);
for i=1:length(alpha)
    idx = L_N > VaR(i, 4);
    ES(i, 4) = mean(L_N(idx));
end

% plotting
figure(1)
plot(alpha, VaR, '-s', 'linewidth', 1.5); grid on;
legend('\nu = 3', '\nu = 10','\nu = 50', 'Normal')
xlabel('$\alpha$', 'interpreter', 'latex')
ylabel('$VaR_\alpha$', 'interpreter', 'latex')
set(gca, 'fontsize', 15)

figure(2)
plot(alpha, VaR_m, '-s', 'linewidth', 1.5); grid on;
legend('\nu = 3', '\nu = 10','\nu = 50', 'Normal')
xlabel('$\alpha$', 'interpreter', 'latex')
ylabel('$VaR_\alpha$', 'interpreter', 'latex')
set(gca, 'fontsize', 15)

figure(3)
plot(alpha, ES, '-s', 'linewidth', 1.5); grid on;
legend('\nu = 3', '\nu = 10','\nu = 50', 'Normal')
xlabel('$\alpha$', 'interpreter', 'latex')
ylabel('$VaR_\alpha$', 'interpreter', 'latex')
set(gca, 'fontsize', 15)

% section 2 -(a)
close all
figure(4)
plot(alpha, VaR-VaR_m, '-s', 'linewidth', 1.5); grid on;
legend('\nu = 3', '\nu = 10','\nu = 50', 'Normal')
xlabel('$\alpha$', 'interpreter', 'latex')
ylabel('$VaR_\alpha-VaR_{mean}$', 'interpreter', 'latex')
% axis([min(alpha) max(alpha) 0 0.03])
set(gca, 'fontsize', 15)