function pdf = get_combined_pdf(x)
    m = 5; 
    K = 15;
    Delta = 0.3;
    Gamm = sqrt(1.8291);
    A0 = 0.5816;
    Km = (m^m) / gamma(m); % See the calculations and the FTR Paper to get Km and K_Gamm
    K_Gamm = (Gamm^2) / (2*A0^(Gamm^2));
    sigma_thz = 1/sqrt(2*(1+K));
    sum_p = 0;
    for p = 0:40
        fG = (x.^((Gamm^2/2)-1)) .* igamma(-Gamm^2/2+p+1,x/(A0^2 * 2*sigma_thz^2)) ./ (gamma(p+1) .* (2*sigma_thz.^2).^(Gamm.^2/2));
        term = K^p * d_calc(p,m,K,Delta) .* fG ./ factorial(p);
        sum_p = sum_p + term;
    end
    pdf = sum_p * Km * K_Gamm;
end

