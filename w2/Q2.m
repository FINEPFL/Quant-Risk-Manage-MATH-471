%% Exercise 2 - Delta Hedge Call Option
clear 
clc
m = 10000;
mu =[0; 0];
Sigma = [(10^-3)^2, -0.5*10^-3*10^-4;
         -0.5*10^-3*10^-4, (10^-4)^2];
     
RV = mvnrnd(mu, Sigma, m); % Generate the multivarite normal RV
X1 = RV(:,1,end);
X2 = RV(:,2,end);

t = 0;
T = 0.5;
rt = 0.05;
sigmat = 0.2;
St = 100;
K = 100;
Delta = 1/252; % 252 trading day in 1 year

%% Monte Carlo Simulation
[Call, Put] = blsprice(St, K, rt, T-t, sigmat);
[CD,PD] = blsdelta(St,K,rt,T-t,sigmat,0);
D = -CD; % The humber of stocks that we need to buy (here since it is
% negative, we will short them).
Vt = Call + D*St; % The value of our position today
Z1t = log(St);
Z2t = sigmat;
Z1T = Z1t + X1;
Z2T = Z2t + X2;
ST = exp(Z1T);
sigmaT = Z2T;
[CallT, PutT] = blsprice(ST, K, rt, T-t-Delta, sigmaT);
VT = CallT + D*ST; % Value of the position after one day
L = -(VT - Vt);
L = sort(L,'descend');
meanL = mean(L);
% For alpha 0.95
alpha = 0.95;
N = m*0.05;
VaR_MC95 = L(N) 
VaRmean_MC95 = VaR_MC95 - mean(L)
ESa_95 = mean(L(1:N))
% For alpha 0.99
alpha = 0.99;
N = m*0.01;
VaR_MC99 = L(N) 
VaRmean_MC99 = VaR_MC99 - mean(L)
ESa_MC99 = mean(L(1:N))

%% Linearized Loss method
t = 0;
T = 0.5;
rt = 0.05;
sigmat = 0.2;
St = 100;
K = 100;
Delta = 1/252; % 252 trading day in 1 year

[CT,    PT] = blstheta(St,K,rt,T-t,sigmat,0);
[CD,PD] = blsdelta(St,K,rt,T-t,sigmat,0);
V = blsvega(St,K,rt,T-t,sigmat,0);

LinearizedL = - (CT*Delta + CD*St*X1 + V*X2 + D*St*X1);
LinearizedL = sort(LinearizedL,'descend');
meanLinearizedL = mean(LinearizedL);

% For alpha 0.95
alpha = 0.95;
N = m*0.05;
VaR_MC95L = LinearizedL(N) 
VaRmean_MC95L = VaR_MC95L - mean(LinearizedL)
ESa_MC95L = mean(LinearizedL(1:N))

% For alpha 0.99
alpha = 0.99;
N = m*0.01;
VaR_MC99L = LinearizedL(N) 
VaRmean_MC99L = VaR_MC99L - mean(LinearizedL)
ESa_MC99L = mean(LinearizedL(1:N))

%% Variance-Covariance Method
c = CT*Delta;
b = [0;
     V];
% For alpha = 0.95
alpha = 0.95;
VaR_VC95 = - c - b'*mu + sqrt(b'*Sigma*b)*norminv(alpha)
Varmean_VC95 = VaR_VC95 - (- c - b'*mu)
ESa_VC95 = - c - b'*mu + sqrt(b'*Sigma*b)*(normpdf((norminv(alpha)))/(1-alpha))
% For alpha = 0.99
alpha = 0.99;
VaR_VC99 = - c - b'*mu + sqrt(b'*Sigma*b)*norminv(alpha)
Varmean_VC99 = VaR_VC99 - (- c - b'*mu)
ESa_VC99 = - c - b'*mu + sqrt(b'*Sigma*b)*(normpdf((norminv(alpha)))/(1-alpha))




