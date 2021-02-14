function h = get_hf_thz(m_inp, K_inp, Delta_inp, iter_inp)
m = m_inp;
K = K_inp;
Delta = Delta_inp;
Sigma=1/sqrt(2*(1+K));
iter = iter_inp;
V1=0.5*Sigma*(sqrt(2*K*(1+Delta))+sqrt(2*K*(1-Delta))); 
V2=0.5*Sigma*(sqrt(2*K*(1+Delta))-sqrt(2*K*(1-Delta))); 

zeta = random('nakagami',m,1,1,iter); %distribution for zeta
X = Sigma*randn(1,iter); %Normal
Y = Sigma*randn(1,iter); 
phi1 = 2*pi*rand(1,iter); % Uniform
phi2 = 2*pi*rand(1,iter);

% FTR Channel
h = sqrt(zeta).*V1.*exp(1i*phi1)+sqrt(zeta).*V2.*exp(1i*phi2) + X + 1i*Y;


