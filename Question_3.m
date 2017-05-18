%% Question 3 for Exerise 9
% Mengjie ZHAO
% Tianxiao MA
%%
clear
clc
Price = xlsread('MicrosoftPrice.csv');
%% Problem 1
negt_logret = - log(Price(2:end, 3) ./ Price(1:end -1, 3));
n = 0.04/0.0004;
x = repmat(0.0004,1,n);
threshold = cumsum(x);
mean_excess = [];
filtered_mean_excess = [];
for i = 1:n
    m = 0.0004*i;
    filtered_mean_excess = negt_logret(find(negt_logret > m)) - m;
    mean_excess(i) = mean(filtered_mean_excess);
end

figure(1)
plot(threshold, mean_excess, '.', 'markersize', 15);
xlabel('Threshold','interpreter','latex','fontsize',16)
ylabel('Mean Excess','interpreter','latex','fontsize',16)

%% Problem 2
w = [];
w = negt_logret(find(negt_logret>0.01)) - 0.01;
parm = gpfit(w);
p = 0.125/0.0001;
y = repmat(0.0001,1,p);
Threshold = cumsum(y);
fitted_pdf = gppdf(Threshold,parm(1),parm(2),0.01);
fitted_cdf = cumsum(fitted_pdf.*0.0001);
[f,v] = ecdf(w);

figure(2)
plot(Threshold(100:end),fitted_cdf(100:end),'LineWidth',3);
hold on;
plot(v+0.01,f,'.','markersize',15);
xlabel('x (negative return) ','interpreter','latex','fontsize',16);
ylabel('$F_u(x-u)$','interpreter','latex','fontsize',16);
h = legend('Fitted Distribution','Empirical Distribution','Location','northeast');
set(h,'Interpreter','latex');
hold off;