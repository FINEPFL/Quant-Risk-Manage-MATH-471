%% Assignment 9 - Quantitative Risk Management
%% Exercise 3
close all
clear all
clc

%% Loading Data of MSFT from May 6, 2010 to May 5, 2015.
raw_data = load('MicrosoftPrice.csv');
% msft = xlsread('MSFT.xlsx','D:D');
msft = raw_data(:, 3);
neg_log = -log(msft(2:end)./msft(1:end-1));      % negative log return of microsoft price

%% Mean Excess Plot - MSFT negative returns
n = 0.04/0.0004;
x = repmat(0.0004,1,n);
threshold = cumsum(x);
mean_excess_neg_log = [];
for i= 1:n
    filtered_excess_neg_log = [];
    t = 0.0004*i;
    filtered_excess_neg_log = neg_log(find(neg_log>t)) - t;
    mean_excess_neg_log(i) = mean(filtered_excess_neg_log);
end

% plot(threshold,mean_excess_neg_log,'^','markersize',10);
% xlabel('Threshold','interpreter','latex','fontsize',16)
% ylabel('Mean Excess','interpreter','latex','fontsize',16)

%% Fit a GPD model to the excess distribution over the threshold by maximum likelihood estimation
close all
w = [];
w = neg_log(find(neg_log>0.01)) - 0.01;
parameter = gpfit(w);
distribution = [];
m = 0.125/0.0001;
x = repmat(0.0001,1,m);
threshold = cumsum(x);
fitted_distribution = gpcdf(threshold,parameter(1),parameter(2),0.01);
[f,v] = ecdf(w);
figure
plot(threshold(100:end),fitted_distribution(100:end),'LineWidth',3);
hold on;
plot(v+0.01,f,'.','markersize',10);

xlabel('x (negative return) ','interpreter','latex','fontsize',15);
ylabel('$F_u(x-u)$','interpreter','latex','fontsize',15);
h = legend('Fitted Distribution','Empirical Distribution','Location','northeast');
set(h,'Interpreter','latex');



















