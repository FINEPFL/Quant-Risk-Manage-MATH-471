function f = neg_likelihood(vec, x_t, sigma_0, n)
    % when the innovation is assumed to be normal
    % starting points
    x = vec(1);
    y = vec(2);
    z = vec(3);
    
    Xt = x_t;
    sigma0 = sigma_0;
    N = n;
    Sigma = sigma0 * ones(N+1, 1);
    
    for i=2:N
        Sigma(i) = sqrt(x + y.*(Xt(i-1).^2) + z.*(Sigma(i-1).^2));
    end
    % the likelihood function (to be sumed)
    % from lecture slide 15
    l = normpdf(Xt(2:end)./Sigma(2:end-1))./Sigma(2:end-1);
    f = - sum(log(l));
    
end

    