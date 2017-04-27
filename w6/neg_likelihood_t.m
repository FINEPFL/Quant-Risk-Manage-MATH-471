function f = neg_likelihood_t(vec, x_t, sigma_0, n)
    % when the innovation is assumed to be t student distributed
    % starting points
    x = vec(1);
    y = vec(2);
    z = vec(3);
    v = vec(4);
    
    Xt = x_t;
    sigma0 = sigma_0;
    N = n;
    Sigma = sigma0 * ones(N+1, 1);
    
    for  i = 2:N
         Sigma(i) = sqrt(x + y.*(Xt(i-1).^2) + z.*(Sigma(i-1).^2));
    end
    
    % from lecture slide 15
    % note the added reg
    reg = sqrt((v-2)/v);
    l = tpdf(Xt(2:end)./(reg.*Sigma(2:end-1)), v)./(reg.*Sigma(2:end-1));
    f = - sum(log(l));

end