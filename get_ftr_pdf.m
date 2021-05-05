function pdf = get_ftr_pdf(x)
    m = 5; 
    K = 15;
    Delta = 0.3;
    Km = (m^m) / gamma(m); % See the calculations and the FTR Paper to get Km and Kj
    sigma_thz = 1/sqrt(2*(1+K));
    sum_p = 0;
    for p = 0:40
        fG = (x.^p) .* exp(-x./(2*sigma_thz^2)) ./ (gamma(p+1) * (2*sigma_thz^2)^(p+1));
        term = K^p * d_calc(p,m,K,Delta) .* fG ./ factorial(p);
        sum_p = sum_p + term;
    end
    pdf = sum_p * Km;
end

