% MEGN_Demo
% Minimal working example for the finite-memory EGN (MEGN) model.
%
% This demo compares the standard EGN estimate with the MEGN correction for
% different CCDM shaping blocklengths and mapping strategies (Fig. 9(b) of the paper).
% Kernel files are loaded from kernel_cache when available.
close all;
clearvars;

%% Link and model parameters
sys.Rs       = 32e-3;       % Symbol rate [TBaud]
sys.Ts       = 1/sys.Rs;    % Symbol period [ps]
sys.rollOff  = 0.05;        % Root-raised-cosine roll-off factor
sys.c        = 299792.458;  % Speed of light [km/s]
sys.waveLen  = 1550;        % Reference wavelength [nm]
sys.alpha_dB = 0.22;        % Fiber attenuation [dB/km]
sys.alpha    = sys.alpha_dB/20*log(10); % Fiber attenuation [1/km]
sys.Dp       = 16.7;        % Dispersion parameter [ps/nm/km]
sys.beta2    = sys.Dp*sys.waveLen.^2/(2*pi*sys.c); % Dispersion [ps^2/km]
sys.gamma    = 1.3;         % Nonlinear parameter [1/W/km]
sys.Ls       = 100;         % Span length [km]
sys.Ns       = 10;          % Number of spans used in this demo

model.Nf1f2   = 2e2*4;      % Number of f1/f2 samples per dimension
model.Nf      = 2e2/4+1;    % Number of evaluated f points; odd so f = 0 exists
model.tau_max = 50;         % MEGN memory length in symbols

%% Signal statistics
mappingGP = [1,2,4];  % 1D, 2D, and 4D mapping cases
blkVec = [10:10:100, 200:100:1000, 2000:1000:10000];
blkGP = 1:length(blkVec);

ampSet = [1,3,5,7];
Pamp = [0.4,0.3,0.2,0.1]; % Amplitude probability distribution

% Normalize amplitudes so the complex QAM average symbol energy is 1.
Ea2_unNorm = sum(Pamp .* ampSet.^2);
normFactor = 1 / sqrt(2 * Ea2_unNorm);
ampSet = ampSet * normFactor;

Ea2 = sum(Pamp .* ampSet.^2);
Ea4 = sum(Pamp .* ampSet.^4);
Ea6 = sum(Pamp .* ampSet.^6);

Ex2 = 2 * Ea2;
P_chn = Ex2;
Ex4 = 2 * Ea4 + 2 * Ea2^2;
Ex6 = 2 * Ea6 + 6 * Ea4 * Ea2;
Phi  = Ex4/Ex2^2-2;
Psi  = Ex6/Ex2^3-9*Ex4/Ex2^2+12;

%% Baseline EGN calculation
f0 = sort(sys.Rs.*linspace(0,1,401)-0.5*sys.Rs);
[eta_SCI_EGN,~,~,~] = calculateEGNSCI(f0,Phi,Psi,400,...
    sys.Ts,sys.rollOff,sys.Rs,sys.alpha,...
    sys.beta2,sys.gamma,sys.Ls,sys.Ns);

%% MEGN correction
megnModel = MEGNModel(sys,model);
eta_SCI_MEGN_Corr = zeros(length(blkGP),length(mappingGP));

for mapIdx = 1:length(mappingGP)
    for blkIdx = 1:length(blkGP)
        fprintf("\nBlocklength case %d/%d, %dD mapping\n", ...
            blkIdx,length(blkGP),mappingGP(mapIdx))
        blkLen = blkVec(blkIdx);

        fprintf("  Computing analytical symbol correlations... ")
        tic
        Cov = calculateSymbolCorrelations(Ea2,Ea4,Ea6,blkLen,mappingGP(mapIdx),model.tau_max);
        toc

        fprintf("  Applying MEGN correction... ")
        tic
        [eta_SCI,~] = megnModel.eta_SCI_MEGN_Cal(P_chn,Cov);
        toc

        eta_SCI_MEGN_Corr(blkIdx,mapIdx) = eta_SCI.corr;
    end
end

eta_SCI_MEGN = eta_SCI_EGN+eta_SCI_MEGN_Corr;

%% Convert eta to optimum effective SNR
F = 6;                   % EDFA noise figure [dB]
G = sys.Ls*sys.alpha_dB; % EDFA gain [dB]
lambda = sys.waveLen;    % Laser wavelength [nm]
B = sys.Rs*1e12;         % Receiver bandwidth [Hz]
SNR_opt_MEGN = etaToSNR(eta_SCI_MEGN,F,G,lambda,B,sys.Ns);
SNR_opt_EGN = etaToSNR(eta_SCI_EGN,F,G,lambda,B,sys.Ns);

%% Plot eta
lineStyles = {'-','-.',':'};
figure
yl = yline(10*log10(eta_SCI_EGN),'-.','EGN Model');
yl.LabelHorizontalAlignment = 'left';
yl.DisplayName = 'EGN Model';
hold on
for mapIdx = 1:length(mappingGP)
    legendText = [num2str(sys.Ns), ' x ', num2str(sys.Ls), ...
        ' km, ', num2str(mappingGP(mapIdx)), 'D mapping 64-QAM, MEGN'];
    semilogx(blkVec,10*log10(eta_SCI_MEGN(:,mapIdx)),...
        'LineStyle',lineStyles{mapIdx},'DisplayName',legendText,"LineWidth",2)
end
xlabel('CCDM shaping blocklength')
ylabel('\eta [dB (1/W^2)]')
legend
set(gca, 'XScale', 'log');
grid on

%% Plot optimum effective SNR
figure
yl = yline(SNR_opt_EGN,'-.','EGN Model');
yl.LabelHorizontalAlignment = 'left';
yl.DisplayName = 'EGN Model';
hold on
for mapIdx = 1:length(mappingGP)
    legendText = [num2str(sys.Ns), ' x ', num2str(sys.Ls), ...
        ' km, ', num2str(mappingGP(mapIdx)), 'D mapping 64-QAM, MEGN'];
    semilogx(blkVec,SNR_opt_MEGN(:,mapIdx),...
        'LineStyle',lineStyles{mapIdx},'DisplayName',legendText,"LineWidth",2)
end
xlabel('CCDM shaping blocklength')
ylabel('Effective SNR [dB]')
legend
set(gca, 'XScale', 'log');
grid on
