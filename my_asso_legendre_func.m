function val = my_asso_legendre_func(m,v,z)
if(m >= 0)
if(m <= v)
term_1 = 1/pi;
for i = 1:m
    term_1 = term_1 * (v + i);
end
int_term = @(phi) ((z + sqrt(z.^2 - 1) .* cos(phi)).^v) .* cos(m.*phi);
term_2 = integral(int_term,0,pi);
val = (1i)^(-mod(m,4)) * term_1 * term_2;
else
val = 0;
end
%val = term_1 * term_2;
else 
mu = -m;
const_term = ((z^2 - 1)^(mu/2)) / (2^mu * sqrt(pi) * gamma(mu + 0.5));
int_term = @(t) ((1 - t.^2).^(mu-0.5))./((z + t.*sqrt(z^2 - 1)).^(mu-v));
val = integral(int_term,-1,1);
val = (1i)^(-mod(m,4)) * const_term * val;
end
end


