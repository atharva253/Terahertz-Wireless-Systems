tic
clear all;
close all;
freq = 275*10^9;
eta = 4;
c= 3*10^8;
gt=10^5.5;
gr=gt;
d=30;
lambda = c/freq;
k=THz_Pathloss(freq);
path_gain_thz= ((c*sqrt(gt*gr))./(4*pi*freq*d)).*exp(-0.5*k*d); % hL = hAL * hFL
signal_bw = 10*10^9;
noise_psd = 3.8*10^(-17);
noise_power = noise_psd*signal_bw;
signal_power_dbm= -20:10:60;
signal_power_lin = 10.^(signal_power_dbm/10)*10^(-3);
power_noise_ratio_db = 20:10:60;    %signal_power_dbm-noise_power_dbm;
power_noise_ratio= 10.^(power_noise_ratio_db/10);
rate_sim_loop = [];
SNR_loop = [];
SNR_num_loop = [];
rate_num_loop = [];
iter=10^6;
    
    for i=power_noise_ratio
    % SNR_0 = actual_SNR * |hL|^2
    y_0_thz= i*abs(path_gain_thz)^2;
    %---------------------------------------------------------------------------------%
    % In this portion, we are modelling the the pointing error hp
    %---------------------------------------------------------------------------------%
    a = 0.05;                         %radius of the reception antenna effective area = 5cm
    wd_1 = 2.5*(d/1000);             %transmission beam footprint radius at distance d1
    u = (sqrt(pi)/sqrt(2))*(a/wd_1)  ;                   
    weq = sqrt(wd_1^2*((sqrt(pi)*erf(u))/(2*u*exp(-u^2))));
    sigma_s = 0.05;
    s =   abs(erf(u))^2;
    % In this code, phi = gamma^2 and s = A0.
    phi = weq^2/(2*sigma_s^2);     %0.2;% 2^(2.5); %0.006513 
    %fraction of the collected power when the transceivers antennas are fully-aligned aka A0 in RIS Paper    s =   abs(erf(u))^2;  %%0.0032;% 0.0032;%0.0032                           
    hp_val= 0:1/iter:s;
    pdf_hp = (phi./(s.^phi)).*hp_val.^(phi-1);
    hp = randpdf(pdf_hp,hp_val,[1,iter]);
    %------------------------------------------------------------------
    %  Now we begin the modelling hf. 
    %--------------------------------------------------------------------
     hf = get_hf_thz(5,15,0.3,iter);
    %%%% Simulation
    
    %g= abs(hp).^2 .* abs(hf).^2 .* y_0_thz;
    g= abs(hf).^2 .* y_0_thz;
    SNR = mean(g,'all');
    SNR_loop = [SNR_loop SNR];
    rate_sim = mean(log2(1+g));    
    rate_sim_loop = [rate_sim_loop rate_sim];           
    
    %%%%% Numerical 
    pdf_thz = @(x) get_ftr_pdf(x);
    SNR_term = @(x) x .* pdf_thz(x);
    SNR_num = integral(SNR_term,0,Inf) * y_0_thz;
    SNR_num_loop = [SNR_num_loop SNR_num];
    rate_num_term = @(x) log2(1+x*y_0_thz) .* pdf_thz(x);
    rate_num = abs(integral(rate_num_term,0,Inf));
    rate_num_loop = [rate_num_loop rate_num]; 
   end
 
 figure(1)
 grid on
 plot(power_noise_ratio_db,rate_sim_loop,'b--')
 xlabel('$P/N_0$ (dB)','FontWeight','normal','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
 %ylabel('Capacity (Bits/Sec/Hz)')
 ylabel('Capacity(Bits/Sec/Hz)','FontWeight','bold','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
%  legend('Simulation','Numerical','Lower Bound','Location','best')
% legend('Capacity(Simulation)','Location','best')
%  title('Capacity') 
hold on

grid on
 plot(power_noise_ratio_db,rate_num_loop,'r')
 xlabel('$P/N_0$ (dB)','FontWeight','normal','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
 %ylabel('Capacity (Bits/Sec/Hz)')
 ylabel('Capacity(Bits/Sec/Hz)','FontWeight','bold','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
%  legend('Simulation','Numerical','Lower Bound','Location','best')
 legend('Capacity(Simulation)','Capacity(Numerical)','Location','best')
%  title('Capacity') 
title('Capacity: FTR only') 

figure(2)
grid on
plot(power_noise_ratio_db,10*log10(abs(SNR_loop)),'b--')
xlabel('$P/N_0$ (dB)','FontWeight','normal','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
 %ylabel('Capacity (Bits/Sec/Hz)')
 ylabel('SNR (dB)','FontWeight','bold','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
%  legend('Simulation','Numerical','Lower Bound','Location','best')
 %legend('SNR(Simulation)','Location','best')
  

hold on
grid on
plot(power_noise_ratio_db,10*log10(abs(SNR_num_loop)),'r')
xlabel('$P/N_0$ (dB)','FontWeight','normal','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
 %ylabel('Capacity (Bits/Sec/Hz)')
 ylabel('SNR (dB)','FontWeight','bold','Color','k','FontSize',12,'Fontname', 'Arial','Interpreter', 'latex')
%  legend('Simulation','Numerical','Lower Bound','Location','best')
 legend('SNR(Simulation)','SNR(Numerical)','Location','best')
  title('SNR: FTR only') 
toc