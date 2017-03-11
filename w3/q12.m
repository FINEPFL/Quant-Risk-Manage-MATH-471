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

%% march 11 2013 to march 10 2016 corresbonds to 503 and 1259
getReturn = @(P) log(P(2:end, :)./P(1:end-1, :));
Return    = getReturn(Price);
b_ = [100;
      100*Price(505, 1)/ Price(505, 2);
      100*Price(505, 1)/ Price(505, 3);];
     
for i=1:(size(Return, 1)-503)
    mu    = mean(Return(i:i+503, :))';
    Sigma = cov(Return(i:i+503, :));
    b = b_ .* Price(i+504, :)';
    L(i) = -b' * Return(i+503, :)';
    VaR_95(i) = -b' * mu + sqrt(b' * Sigma * b) * norminv(0.95);
    VaR_99(i) = -b' * mu + sqrt(b' * Sigma * b) * norminv(0.99);
    breach_95(i) = VaR_95(i) > L(i);
    breach_99(i) = VaR_99(i) > L(i);
end
figure(1)
plot(L, 'k'); hold on;
plot(VaR_95, 'b'); plot(VaR_99, 'r'); 
