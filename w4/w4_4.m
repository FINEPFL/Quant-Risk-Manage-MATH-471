clearvars; close all; clc

rng(1)
N = 10000;
A = [ 1  0 0 0;
      1  1 0 0;
     -1  2 3 0;
      1 -1 1 1];

x = trnd(5, N, 4);
X = (A * x')';

% calculate Loss, simply sum them up 
L = sum(X, 2);

% calculate VaR at 0.95
VaR = quantile(L, .95)
[vecs, vals] = eig(cov(X));
mu = mean(X);

Y = (vecs' * (X - ones(length(X), 1) * mu)')';

% approximate X using first two principle components
ap_X = ones(length(X), 1) * mu + ([vecs(:, 4) vecs(:, 3)] * [Y(:, 4) Y(:, 3)]')';
ap_L = sum(ap_X, 2);
ap_VaR = quantile(ap_L, 0.95)