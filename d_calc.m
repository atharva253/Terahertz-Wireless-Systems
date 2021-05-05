function d = d_calc(n,m,K,Delta)
    NCR = @(N,R) factorial(N)/(factorial(R) * factorial(N-R));
    d = 0;
    for k = 0:n
        suml = 0;
        for l = 0:k
            term_l1 = NCR(k,l)*(gamma(n+m+2*l-k)*exp(1j*pi*(2*l-k)/2)* ((m+K)^2 - (K*Delta)^2)^(-1*(n+m)/2));
            z = (m+K)/sqrt((m+K)^2 - (K*Delta)^2);
           
            Mu = k-2*l;
            Nu = n+m-1;
           
          P = my_asso_legendre_func(Mu,Nu,z); 
          %%keyboard;
          %  P = (1/gamma(1-Mu)) .* ((z+1)/(z-1)).^(Mu/2) .* hypergeom([-1*Nu, Nu+1],1-Mu,(1-z)/2);
            term_l2 = P;
            suml = suml + (term_l1 * term_l2);
        end
        term_k = NCR(n,k) * (Delta/2)^k * suml;
        d = d + term_k;
    end
end
            