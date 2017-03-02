clear; clc

S_t = 100; lambda = 1; d = 1; mu_X = 0; SD_X = 0.01; f = [3, 10, 50]; N = 10000;

alpha =  SD_X ./ sqrt(f ./ (f-2));
T(:, 1) = trnd (f(1, 1), N, 1);
T(:, 2) = trnd (f(1, 2), N, 1);
T(:, 3) = trnd (f(1, 3), N, 1);
T(:, 4) = randn (N, 1);

X(:, 1) = alpha(:, 1) * T(:, 1);
X(:, 2) = alpha(:, 2) * T(:, 2);
X(:, 3) = alpha(:, 3) * T(:, 3);
X(:, 4) = SD_X * T(:, 4);

L_t = -lambda * S_t * (exp(X)-1);

figure(1)
[N_1, X] = hist(L_t(:, 1), 150);
bar(X, N_1/(sum(N_1)*diff(X(1:2))), 1);
x = -8:0.1:8;
y = normpdf(x, 0, 1);
hold on
plot(x, y, 'r', 'linewidth', 3);
axis([-8 8 -inf inf]);
hold off
grid on;
xlabel('$L(t, t+\Delta)$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('Probability density')
title('Freedom $\nu=3$', 'interpreter', 'latex', 'fontsize', 15)
set(gca, 'fontsize', 12)

figure(2)
[N_2, X] = hist(L_t(:, 2), 100);
bar(X, N_2/(sum(N_2)*diff(X(1:2))), 1);
hold on
plot(x, y, 'r', 'linewidth', 3);
hold off
grid on;
xlabel('$L(t, t+\Delta)$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('Probability density')
title('Freedom $\nu=10$', 'interpreter', 'latex', 'fontsize', 15)
set(gca, 'fontsize', 12)    

figure(3)
[N_3, X] = hist(L_t(:, 3), 100);
bar(X, N_3/(sum(N_3)*diff(X(1:2))), 1);
hold on
plot(x,y, 'r', 'linewidth', 3);
hold off
grid on;
xlabel('$L(t, t+\Delta)$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('Probability density')
title('Freedom $\nu=50$', 'interpreter', 'latex', 'fontsize', 15)
set(gca, 'fontsize', 12)

figure(4)
[N_4, X] = hist(L_t(:, 4), 100);
bar(X, N_4/(sum(N_4)*diff(X(1:2))), 1);
hold on
plot(x,y, 'r', 'linewidth', 3);
hold off
grid on;
xlabel('$L(t, t+\Delta)$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('Probability density')
title('Normal', 'interpreter', 'latex', 'fontsize', 15)
set(gca, 'fontsize', 12)