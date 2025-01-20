function H_Ric=Ric2deig_model(K_dB,N1,N2) %(output of func.=name)
% Rician channel model
% Input : K_dB = K factor[dB]
% Output: H = Channel vector
 K = 10^(K_dB/10);
 m = sqrt(K/(K+1));  s = sqrt(1/2*(K+1));
  H_Ric =( ( s* randn(N1,N2)+ m) + 1i *( randn(N1,N2)*s+ m));% convert from Rayleigh to Rician
 
 
