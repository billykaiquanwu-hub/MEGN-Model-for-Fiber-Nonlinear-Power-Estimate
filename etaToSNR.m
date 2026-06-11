function [SNR_opt] = etaToSNR(eta,F,G,lambda,B,Ns)
%ETATOSNR Convert nonlinear coefficient eta to optimum effective SNR.

h = 6.62607015*1e-34; % Planck's constant [J/Hz]
c = 299792.458;       % Speed of light [km/s]
v = c*1e12/lambda;    % Carrier frequency [Hz]

var_ASE = h*v*10^(F/10)*10^(G/10)*B.*Ns;
P_opt = nthroot(var_ASE/2./eta,3);
SNR_opt = 10*log10(P_opt./(3/2*var_ASE));
end
