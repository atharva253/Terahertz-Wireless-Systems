function k=THz_Pathloss(freq)
     
c = 3*10^8;
%%freq = 275*10^9
% lambda = c/freq;
% % eta = 2
% d1 = [1:1000]

% gt = 10;%10^(2.7)%dBi
% gr =10;% 10^(2.7)%
q1 = 0.2205;
q2 = 0.1303;
q3 = 0.0294;
q4 = 0.4093;
q5 = 0.0925;
q6 = 2.014;
q7 = 0.1702;
q8 = 0.0303;
q9 = 0.537;
q10 = 0.0956;
c1 = 5.54* 10^-37  ;               % Per Hz^3
c2 = -3.94*10^-25  ;               % Per Hz^2
c3 = 9.06*10^-14  ;                % per Hz
c4 = -6.36*10^-3  ;                % Per Hz^3
p1 = 10.835*10^-2;%*d1             % per cm
p2 = 12.664*10^-2;%*d1             % per cm
si = 0.5;
p = 101325  ;                      %Pa
g1 = 6.1121;
g2 = 1.0007;
g3 = 3.46*10^-6;
g4 = 17.502;
g5 = 273.15;
g6 = 32.18;
t = 296;
fi_h = 1013.25 ;                   % hPa (1 hPa = 100 Pa)
%pw = 28.1                        % Using Buck's Equation
pw = g1*(g2+g3*fi_h)*exp((g4*(t-g5))/(t-g6));
v = (si*pw)/(100*p);
x = (q1*v*(q2*v+q3))/((q4*v+q5)^2+((freq/(100*c)-p1)^2));
y = (q6*v*(q7*v+q8))/((q9*v+q10)^2+((freq/(100*c)-p2)^2));
z = (c1*freq^3)+(c2*freq^2)+(c3*freq)+c4;
k = x+y+z;
% if test==0
% for i=1:1:1000
% hl(i) = ((c*sqrt(gt*gr))/(4*pi*freq*d1(i)))*exp(-0.5*k*d1(i))
% beta(i) =sqrt( (gt*gr*lambda^2)./(4*pi*d1(i)).^eta)
% end
% 
% % beta = (gt*gr*lambda^2)./(4*pi*d1).^eta;
% figure(1)
% hold all
% plot(d1,-20*log10(hl),'b') 
% plot(d1,-20*log10((beta)),'r')
% xlabel('Distance')
% ylabel('Path Gain')
% legend ('THz','Free Space')
% end
