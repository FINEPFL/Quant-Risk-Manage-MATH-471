clearvars; close all; clc

handle = fopen('data.csv');
raw_data = textscan(handle, '%f %f %s %f', 'delimiter', ',');
fclose(handle);

i = 1;
numMSFT = 0;
numYHOO = 0;
numITEL = 0;
while i <= length(raw_data{1, 3})
    switch(raw_data{1, 3}{i})
        case 'MSFT'
            numMSFT = numMSFT + 1;
        case 'YHOO'
            numYHOO = numYHOO + 1;
        otherwise
            numITEL = numITEL + 1;
    end
    i = i + 1;
end

msft_price = raw_data{1, 4}(1:numMSFT);
itel_price = raw_data{1, 4}(1+numMSFT:numMSFT+numITEL);
yhoo_price = raw_data{1, 4}(1+numMSFT+numITEL:end);

Price = [msft_price itel_price yhoo_price];

%% march 11 2013 to march 10 2016 corresponds to 503 and 1259
getReturn = @(P) log(P(2:end, :)./P(1:end-1, :));
Return    = getReturn(Price);
b_ = [100;
      100*Price(503, 1)/ Price(503, 2);
      100*Price(503, 1)/ Price(503, 3);];
     
for i=1:(size(Return, 1)-503)
    mu    = mean(Return(i:i+503, :))';
    Sigma = cov(Return(i:i+503, :));
    b = b_ .* Price(i+503, :)';
    L(i) = -b' * Return(i+503, :)';
    VaR_95(i) = -b' * mu + sqrt(b' * Sigma * b) * norminv(0.95);
    VaR_99(i) = -b' * mu + sqrt(b' * Sigma * b) * norminv(0.99);
    
    breach_95(i) = VaR_95(i) < L(i);
    breach_99(i) = VaR_99(i) < L(i);

end
figure(1)
xaxis = 1:960;
plot(L, 'k', 'linewidth', 1); hold on;
plot(VaR_95, 'b', 'linewidth', 2); plot(VaR_99, 'r', 'linewidth', 2); 
legend('Loss', 'VaR_{95}', 'VaR_{99}')
xlabel('date', 'interpreter', 'latex')
ylabel('loss', 'interpreter', 'latex')
plot(xaxis(breach_95), VaR_95(breach_95), 'g.', 'markersize', 25)
plot(xaxis(breach_99), VaR_99(breach_99), 'y.', 'markersize', 25)
set(gca, 'fontsize', 15)
times95 = sum(double(breach_95))
times99 = sum(double(breach_99))

%% Q2
% 1
clearvars; close all; clc
p = 0.5;
VaR_95 = geoinv(0.95, p)
alpha  = 0.9:0.0001:0.99;
VaR_alpha = geoinv(alpha, p);
plot(alpha, VaR_alpha, 'linewidth', 1.8)
xlabel('$\alpha$', 'interpreter', 'latex')
ylabel('$VaR_{\alpha}$', 'interpreter', 'latex')
title('$L{\sim}geom(0.5)$', 'interpreter', 'latex')
set(gca, 'fontsize', 15)
% 2
clearvars; close all; clc
alpha  = 0.9:0.0001:0.99;
VaR_X = poissinv(alpha, 1);
VaR_Y = poissinv(alpha, 2);
VaR_Z = poissinv(alpha, 3);
plot(alpha, VaR_X, alpha, VaR_Y, alpha, VaR_Z, 'linewidth', 1.8)
legend('VaR_\alpha(X)', 'VaR_\alpha(Y)', 'VaR_\alpha(L)')
xlabel('$\alpha$', 'interpreter', 'latex')
ylabel('$VaR_{\alpha}$', 'interpreter', 'latex')
title('$L{\sim}Poisson(\lambda)$', 'interpreter', 'latex')
set(gca, 'fontsize', 15)

