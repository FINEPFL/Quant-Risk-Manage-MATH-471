clear all
close all
clc

%% sample mean excess plot
%data = xlsread('msft.csv');
% T = readtable('msft.csv');
raw_data = load('MicrosoftPrice.csv');
p = raw_data(:, 3);
rr = -log(p(1:end-1,1)./p(2:end,1));
r = sort(rr,'ascend');
r=r(r>0); % negative returns (positive losses)
e = ones(length(r)-2,1);
for i=1:length(e)
e(i,1) = sum(r(r>r(i+1))-r(i+1))/length(r(r>r(i+1)));
end
 %% plotting
% figure(1);
% scatter(r(3:end,1),e,linspace(1,10,1),'filled');
% xlabel('Threshold');
% ylabel('Sample Mean Excess');
% title('Sample Mean Excess Plot ? MSFT negative returns');
figure(2)
t=find(r>=0.04);
scatter(r(3:t(1),1),e(1:t(1)-2,1),linspace(1,10,1),'filled');
xlim([0 0.04]);
xlabel('Threshold');
ylabel('Sample Mean Excess');
title('Sample Mean Excess Plot ? MSFT negative returns');
%% fitting GDP model
clear t; t=0.01
clear e;
e = r(r>0.01)-0.01 % exceedance vector
par = gpfit(e);
gdistr = gpcdf(r(r>0.01),par(1),par(2),0.01);
figure(3)
scatter(r(r>0.01),gdistr);
hold on
h = cdfplot(r(r>0.01));
set(h,'linewidth',3.5);
xlim([0 max(r)]);
xlabel('Negative return (>0.01)');
ylabel('Empirical and GP distr.');
title('Estimated Tail Distribution')