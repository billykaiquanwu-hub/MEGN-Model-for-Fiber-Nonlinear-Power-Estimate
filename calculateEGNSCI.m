function [eta_SCI_EGN,eta_SCI_GN,G_SCI_EGN,G_SCI_GN] = calculateEGNSCI(P_chn,f0,Phi,Psi,Nsmpf1f2,T,rollOff,Rs,alpha,beta2,gamma,Ls,Ns)
%CALCULATEEGNSCI Calculate baseline EGN/GN SCI coefficients.
%
% Inputs follow the notation used in MEGN_Demo.m. Rs is expressed in TBaud.

f1Num = Nsmpf1f2;
f2Num = Nsmpf1f2;
f1Mtx = repmat(Rs.*(linspace(0,1,f1Num)')-0.5*Rs.*ones(f1Num,1),1,f2Num);
f2Mtx = repmat(Rs.*linspace(0,1,f2Num)-0.5*Rs.*ones(1,f2Num),f1Num,1);

spec1Mtx = RRC_spectra(f1Mtx,T,rollOff);
spec2Mtx = RRC_spectra(f2Mtx,T,rollOff);

[G_SCI_EGN,G_SCI_GN] = G_SCI_Cal(P_chn,f0,f1Mtx,f2Mtx,spec1Mtx,spec2Mtx,T,rollOff,Phi,Psi,Rs,alpha,beta2,gamma,Ls,Ns);

P_tot  = 2*P_chn;% total launch power 
eta_SCI_EGN = mean(G_SCI_EGN)*Rs/(P_tot^3);
eta_SCI_GN = mean(G_SCI_GN)*Rs/(P_tot^3);
end

function [pulseSpec] = RC_spectra(freq_vec,T,rollOff)
sz = size(freq_vec);
freq_vec = freq_vec(:);
pulseSpec = T.*zeros(length(freq_vec),1);

bound1 = (1-rollOff)/2/T;
bound2 = (1+rollOff)/2/T;

range1 = abs(freq_vec)<=bound1;
pulseSpec(range1) = 1;

range2 = abs(freq_vec)>bound1 & abs(freq_vec)<bound2;
pulseSpec(range2) = 1/2*(1+...
    cos(pi*T/rollOff*(abs(freq_vec(range2))-(1-rollOff)/2/T)));

pulseSpec = reshape(pulseSpec,sz);
end

function [pulseSpec] = RRC_spectra(freq_vec,T,rollOff)
[pulseSpecRC] = RC_spectra(freq_vec,T,rollOff);
pulseSpec = sqrt(pulseSpecRC);
pulseSpec = pulseSpec./max(pulseSpec(:))*T;
end

function mu = linkFunc(ftms,alpha,beta2,gamma,Ls,Ns)
sz = size(ftms);

zeta = gamma*(1-exp(-2*alpha*Ls+1j*4*pi^2*beta2*Ls.*ftms))./...
    (2*alpha-1j*4*pi^2*beta2.*ftms);

nu = sin(2*beta2*pi^2*Ls*Ns.*ftms)./...
    sin(2*beta2*pi^2*Ls.*ftms).*...
    exp(1j*2*beta2*pi^2*Ls*(Ns-1).*ftms);

nu(ftms==0) = Ns.*exp(1j*2*beta2*pi^2*Ls*(Ns-1).*ftms(ftms==0));

mu = zeta.*nu;
mu = reshape(mu,sz);
end

function [ka1,ka2,ka3] = kaCal(f0,f1Mtx,f2Mtx,spec1Mtx,spec2Mtx,T,rollOff,Rs,alpha,beta2,gamma,Ls,Ns)
Nsmpf = length(f0);
ka1 = zeros(Nsmpf,1);
ka2 = zeros(Nsmpf,1);
ka3 = zeros(Nsmpf,1);

for idx = 1:Nsmpf
    f = f0(idx);

    fcmb1 = f1Mtx + f2Mtx - f;
    spec_fcmb1 = RRC_spectra(fcmb1,T,rollOff);
    argOutA = abs(spec_fcmb1)~=0;

    fcmb2 = f1Mtx - f2Mtx + f;
    spec_fcmb2 = RRC_spectra(fcmb2,T,rollOff);
    argOutB = abs(spec_fcmb2)~=0;

    ftms = (f1Mtx-f).*(f2Mtx-f);
    mu_vec1 = linkFunc(ftms,alpha,beta2,gamma,Ls,Ns);

    ftms = (f1Mtx-f2Mtx).*(f1Mtx-2.*f2Mtx+f);
    mu_vec2 = linkFunc(ftms,alpha,beta2,gamma,Ls,Ns);

    phi1 = 64/27*Rs^3*...
        mean(mean(abs(spec1Mtx).^2.*abs(spec2Mtx).^2.*...
        abs(spec_fcmb1).^2.*argOutA.*abs(mu_vec1).^2))*Rs^2;
    ka1(idx) = phi1;

    phi2 = 64/81*Rs^2*...
        mean(abs(spec1Mtx(:,1)).^2.*...
        abs(mean(spec2Mtx.*conj(spec_fcmb1).*argOutA.*mu_vec1,2)*Rs).^2)*Rs; 
    phi3 = 64/81*Rs^2*...
        mean(abs(spec1Mtx(:,1)).^2.*...
        abs(mean(spec2Mtx.*spec_fcmb2.*argOutB.*mu_vec2,2)*Rs).^2)*Rs;
    ka2(idx) = 5*phi2 + phi3;

    phi4 = 64/81*Rs*...
        abs(mean(mean(spec1Mtx.*spec2Mtx.*conj(spec_fcmb1).*argOutA.*mu_vec1))*Rs^2)^2;
    ka3(idx) =  phi4;
end
end

function [G_SCI_EGN,G_SCI_GN] = G_SCI_Cal(P_chn,f0,f1Mtx,f2Mtx,spec1Mtx,spec2Mtx,T,rollOff,Phi,Psi,Rs,alpha,beta2,gamma,Ls,Ns)
[ka1,ka2,ka3] = kaCal(f0,f1Mtx,f2Mtx,spec1Mtx,spec2Mtx,T,rollOff,Rs,alpha,beta2,gamma,Ls,Ns);

G_SCI_EGN = 2*(ka1+Phi.*ka2+Psi.*ka3)*P_chn^3;
G_SCI_GN = 2*ka1*P_chn^3;
end
