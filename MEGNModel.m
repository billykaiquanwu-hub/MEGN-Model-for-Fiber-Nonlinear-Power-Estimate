classdef MEGNModel
    %MEGNMODEL Finite-memory EGN model for SCI correction.
    %
    % The model loads precomputed kernels from kernel_cache when available.
    % Missing kernels are computed once and saved for later runs.

    properties
        %%% for fiber setup
        sys_Rs        % symbol rate in TBaud (1 TBd = 1000 GBd = 1e12 Bd)
        sys_Ts        % sampling interval, picosecond (1 ps = 1e-12 s)
        sys_rollOff   % roll-off factor of root raised cosine pulse
        sys_c         % Light speed [km/s]
        sys_waveLen   % Reference wavelength lambda [nm]
        sys_alpha_dB  % Fibre attenuation [dB/km]
        sys_alpha     % Fibre attenuation in [1/km]
        sys_Dp        % Dispersion parameter [ps/nm/km]
        sys_beta2     % Dispersion coefficient [ps^2/km]
        sys_gamma     % Nonlinear parameter [1/W/km]
        sys_Ls        % Span length [km]
        sys_Ns        % Number of Spans

        %%% Model parameters
        model_Nf1f2   % sample number for f1, f2 pairs
        model_Nf      % sample number for evaluated PSD f
        model_f0      % specify the Evaluated frequency points and the spectrum
        model_tau_max % sample number for evaluated PSD f

        model_kn_phi_2      % old kernels for two identical time index
        model_kn_phi_3      % old kernels for two identical time index
        model_kn_phi_4      % old kernels for one identical time index
        model_kn_chi_1      % kernels for two distinct time groups
        model_kn_chi_2      % kernels for two distinct time groups
        model_kn_chi_3      % kernels for two distinct time groups
        model_kn_xi_1      % kernels for two distinct time groups
        model_kn_xi_2      % kernels for three distinct time groups
        model_kn_psi_1      % kernels for two distinct time groups
        model_kn_psi_2      % kernels for two distinct time groups
        model_kn_psi_3      % kernels for three distinct time groups

        model_Kn            % Channel-Dependent Functions given kernels
    end

    methods
        function obj = MEGNModel(sys,model)

            %MEGNMODEL Construct an instance of this class.
            obj.sys_c         = 299792.458;
            obj.sys_Rs        = sys.Rs;
            obj.sys_Ts        = 1/obj.sys_Rs;
            obj.sys_rollOff   = sys.rollOff;
            obj.sys_waveLen   = sys.waveLen;
            obj.sys_alpha_dB  = sys.alpha_dB;
            obj.sys_alpha     = sys.alpha_dB/20*log(10);
            obj.sys_Dp        = sys.Dp;
            obj.sys_beta2     = sys.Dp*sys.waveLen.^2/(2*pi*obj.sys_c);
            obj.sys_gamma     = sys.gamma;
            obj.sys_Ls        = sys.Ls;
            obj.sys_Ns        = sys.Ns;

            %%% Model Parameters For estimating SCI correction
            obj.model_Nf1f2   = model.Nf1f2;
            obj.model_Nf      = model.Nf;
            obj.model_tau_max = model.tau_max;

            %%% Folder name for storing/getting MEGN model coefficients
            [folder_name,file_name] = getName(obj);
            list=dir(file_name);
            LL=length(list());
            if LL == 0
                fprintf('File %s not found!\nInitialize MEGN model coefficients\n',file_name);
                if exist(folder_name,'dir')~=7
                    mkdir(folder_name);
                end
                tic
                obj = obj.initialize_Coeffs(file_name);
                toc
            end 


            load(file_name,'f0','kn_phi_2','kn_phi_3','kn_phi_4','kn_chi_1','kn_chi_2',...
                'kn_chi_3','kn_xi_1','kn_xi_2','kn_psi_1','kn_psi_2','kn_psi_3')
            
            obj.model_f0           = f0;
            obj.model_kn_phi_2     = kn_phi_2;
            obj.model_kn_phi_3     = kn_phi_3;
            obj.model_kn_phi_4     = kn_phi_4;
            obj.model_kn_chi_1     = kn_chi_1;
            obj.model_kn_chi_2     = kn_chi_2;
            obj.model_kn_chi_3     = kn_chi_3;
            obj.model_kn_xi_1      = kn_xi_1;
            obj.model_kn_xi_2      = kn_xi_2;
            obj.model_kn_psi_1     = kn_psi_1;
            obj.model_kn_psi_2     = kn_psi_2;
            obj.model_kn_psi_3     = kn_psi_3;

            [obj.model_Kn] = Kn_SCI_MEGN_Cal(obj);
        end

        function [folder_name,file_name] = getName(obj)
            % get the folder name for the considered setup
            if ispc
                bs='\'; % for windows machines
            else
                bs='/'; % non-windows backslash
            end

            sysStr=[
                num2str(obj.sys_waveLen),'nm_',...
                num2str(obj.sys_Rs*1e3),'GBd_',...
                num2str(obj.sys_rollOff, '%.2f'),'rollOff_',...
                num2str(obj.sys_alpha_dB, '%.2f'),'alphdB_',...
                num2str(obj.sys_Dp, '%.2f'),'Dp_',...
                num2str(obj.sys_beta2, '%.2f'),'beta2_',...
                num2str(obj.sys_gamma, '%.2f'),'gamma_',...
                num2str(obj.sys_Ls),'Ls_'...
                num2str(obj.model_tau_max),'tau'];

            folder_name=['.',bs,'kernel_cache',bs,sysStr,bs];
            file_name=[folder_name,'Span_', num2str(obj.sys_Ns),'.mat'];
        end

        function [obj] = initialize_Coeffs(obj,file_name)
            
            %%% Monte-Carlo integration
            tau_max = obj.model_tau_max;

            f0  = sort(1*obj.sys_Rs.*linspace(0,1,obj.model_Nf)-0.5*obj.sys_Rs);

            f1Num = 1*obj.model_Nf1f2; f2Num = 1*obj.model_Nf1f2;
            f1Mtx = repmat(1*obj.sys_Rs.*(linspace(0,1,f1Num)')-0.5*obj.sys_Rs.*ones(f1Num,1),1,f2Num);
            f2Mtx = repmat(1*obj.sys_Rs.* linspace(0,1,f2Num) -0.5*obj.sys_Rs.*ones(1,f2Num),f1Num,1);
            
            spec1Mtx = obj.RRC_spectra(f1Mtx);
            spec2Mtx = obj.RRC_spectra(f2Mtx);

            %%% Evaluated coefficients at different frequency points 
            kn_phi_4 = zeros(1,obj.model_Nf);
            kn_phi_2 = zeros(1,obj.model_Nf);
            kn_phi_3 = zeros(1,obj.model_Nf);
            kn_chi_1 = zeros(obj.model_tau_max,obj.model_Nf);
            kn_chi_2 = zeros(obj.model_tau_max,obj.model_Nf);
            kn_chi_3 = zeros(obj.model_tau_max,obj.model_Nf);
            kn_xi_1 = zeros(obj.model_tau_max,obj.model_Nf);
            kn_xi_2 = zeros(obj.model_tau_max,2*obj.model_tau_max,obj.model_Nf);
            kn_psi_1 = zeros(obj.model_tau_max,obj.model_Nf);
            kn_psi_2 = zeros(obj.model_tau_max,obj.model_Nf);
            kn_psi_3 = zeros(obj.model_tau_max,2*obj.model_tau_max,obj.model_Nf);


            % simple progress bar
            fprintf('Progress:\n');
            fprintf(['\n' repmat('.',1,obj.model_Nf) '\n\n']);

            % For SCI coefficients
            %parfor idx = 1:obj.model_Nf
            for idx = 1:obj.model_Nf

                f = f0(idx);
                %%% pulse spectrum
                %%% f1+f2-f
                fcmb1 = f1Mtx + f2Mtx - f;
                spec_fcmb1 = obj.RRC_spectra(fcmb1); 

                %%% f1-f2+f
                fcmb2 = f1Mtx - f2Mtx + f;
                spec_fcmb2 = obj.RRC_spectra(fcmb2); 

                %%% Link Functions
                ftms = (f1Mtx-f).*(f2Mtx-f);
                mu_vec1 = linkFunc(obj,ftms);

                ftms = (f1Mtx-f2Mtx).*(f2Mtx-f);
                mu_vec2 = linkFunc(obj,ftms);

                specBtVec = spec1Mtx.*spec2Mtx.*conj(spec_fcmb1).*mu_vec1;

                %%% (1) EGN kernels
                kn_phi_4(idx) = 64/81*obj.sys_Rs*abs(mean(mean(specBtVec))*obj.sys_Rs^2)^2;
                kn_phi_2(idx) = 64/81*obj.sys_Rs^2*...
                    mean(abs(spec1Mtx(:,1)).^2.*...
                    abs(mean(spec2Mtx.*conj(spec_fcmb1).*mu_vec1,2)*obj.sys_Rs).^2)*obj.sys_Rs;
                kn_phi_3(idx) = 64/81*obj.sys_Rs^2*...
                    mean(abs(spec1Mtx(:,1)).^2.*...
                    abs(mean(spec2Mtx.*spec_fcmb2.*mu_vec2,2)*obj.sys_Rs).^2)*obj.sys_Rs;


                %%% (2) Kernels for double energy correlation (two time groups)
                for tau = 1:tau_max
                    [kn_chi_1(tau,idx),kn_chi_2(tau,idx),kn_chi_3(tau,idx)] = ...
                        obj.specSumCal2_1(f,specBtVec,f1Mtx,f2Mtx,tau); 

                    [kn_xi_1(tau,idx),kn_psi_1(tau,idx)] = ...
                        obj.specSumCal2_2(f,f1Mtx,f2Mtx,spec1Mtx,spec2Mtx,...
                        spec_fcmb1,mu_vec1,tau);

                    [kn_psi_2(tau,idx)] = obj.specSumCal2_3(f,f1Mtx,f2Mtx,...
                        spec1Mtx,spec2Mtx,spec_fcmb2,mu_vec2,tau);
                end
 

                %%% (3) Triple Energy Correlation (Three time groups)
                for tau1 = 1:tau_max
                    for tau2 = tau1+1:tau1+tau_max
                        [kn_xi_2(tau1,tau2,idx),...
                            kn_psi_3(tau1,tau2,idx)] = ...
                            obj.specSumCal3(f,specBtVec,f1Mtx,f2Mtx,tau1,tau2);
                    end
                end

                % progress bar
                fprintf('\b|\n');
            end 
  
            save(file_name,'f0','kn_phi_4','kn_phi_2','kn_phi_3','kn_chi_1','kn_chi_2',...
                'kn_chi_3','kn_xi_1','kn_xi_2','kn_psi_1','kn_psi_2','kn_psi_3')
        end

        %%% Calculate the Channel Dependent Functions given Kernels
        function [Kn] = Kn_SCI_MEGN_Cal(obj)

            %%% Symbol correlation statistics 
            for idx = 1:obj.model_Nf  
                

                %%% SPT: self-pol. temporal interactions
                %%%%%% Pure Double Energy Correlation 
                Kn.K22_SPT(:,idx) = 5.*obj.model_kn_xi_1(:,idx) ...
                    + 5.*obj.model_kn_psi_1(:,idx)...
                    + 2.*obj.model_kn_psi_2(:,idx) ...
                    - 22.*obj.model_kn_chi_1(:,idx)...
                    - 21.*obj.model_kn_chi_2(:,idx) ...
                    - 5.*obj.model_kn_chi_3(:,idx); 
                Kn.K24_SPT(:,idx) = 4.*obj.model_kn_chi_1(:,idx) ...
                    + 4.*obj.model_kn_chi_2(:,idx)...
                    + 1.*obj.model_kn_chi_3(:,idx); 
                %%%%%% Pure Triple Energy Correlation  
                Kn.K222_SPT(:,:,idx) =  4.*obj.model_kn_xi_2(:,:,idx)...
                    +  2.*obj.model_kn_psi_3(:,:,idx);  


                %%% XP: cross-pol. interactions
                %%%%%% Pure Single Energy Correlation 
                Kn.K24_XP(idx) = 3.*obj.model_kn_phi_4(idx);  
                Kn.K22_XP(idx) = -12.*obj.model_kn_phi_4(idx)...
                    + 5.*obj.model_kn_phi_2(idx)...
                    + 1.*obj.model_kn_phi_3(idx); 


                %%% XPT: cross-pol. temporal interactions 
                %%%%%% Pure Double Energy Correlation  
                Kn.K22_XPT(:,idx) = 4.*obj.model_kn_xi_1(:,idx) ...
                    + 1.*obj.model_kn_psi_1(:,idx)...
                    + 1.*obj.model_kn_psi_2(:,idx)...
                    - 14.*obj.model_kn_chi_1(:,idx) ...
                    - 9.*obj.model_kn_chi_2(:,idx) ...
                    - 1.*obj.model_kn_chi_3(:,idx); 
                Kn.K24_XPT(:,idx) = 2.*obj.model_kn_chi_1(:,idx) ...
                    + 1.*obj.model_kn_chi_2(:,idx); 
                Kn.K222_XPT_1tau(:,idx) = 6.*obj.model_kn_chi_1(:,idx) ...
                    + 5.*obj.model_kn_chi_2(:,idx)...
                    + 1.*obj.model_kn_chi_3(:,idx);  
                %%%%%% Pure Triple Energy Correlation 
                Kn.flagAppxOnly = 0;
                if sum(abs(obj.model_kn_xi_2(:)))==0 || sum(abs(obj.model_kn_psi_3(:)))==0
                    Kn.flagAppxOnly = 1;
                    fprintf("Kernels for Triple Energy not computed. Only support approximate calculation")
                end
                Kn.K222_XPT(:,:,idx) =  5.*obj.model_kn_xi_2(:,:,idx)...
                    +  1.*obj.model_kn_psi_3(:,:,idx);   


            end
            
            %%%%% Get all different SPT, XP, XPT contributions
                Kn.SPT_2D = sum(Kn.K22_SPT(:))...
                    + sum(Kn.K24_SPT(:));
                Kn.SPT_3D = sum(Kn.K222_SPT(:)); 

                Kn.XP_1D = sum(Kn.K22_XP(:))...
                    + sum(Kn.K24_XP(:));

                Kn.XPT_2D =  sum(Kn.K22_XPT(:))...
                    + sum(Kn.K24_XPT(:))...
                    + sum(Kn.K222_XPT_1tau(:));
                Kn.XPT_3D = sum(Kn.K222_XPT(:));
 
        end


        %%% Calculate the NLI correction terms due to symbol dependency Cov
        function [eta_SCI,G_SCI] = eta_SCI_MEGN_Cal(obj,P_tot,Cov)

            %%% Symbol correlation statistics
            tau_max = obj.model_tau_max;

            % For intra-pol correlations  
            Ex2 = Cov.Ex2;
            Cov22_X = Cov.Cov22_X;
            Cov24_X = Cov.Cov24_X;
            Cov222_X = Cov.Cov222_X;

            % Construct the combined triple-energy covariance.
            Cov222_Cmb_X = zeros(size(Cov222_X));
            for tau1 = 1:tau_max
                for tau2 = tau1+1:tau1+tau_max
                    Cov222_Cmb_X(tau1,tau2) = Cov222_X(tau1,tau2)...
                        -Ex2*Cov22_X(tau1)-Ex2*Cov22_X(tau2)...
                        -Ex2*Cov22_X(tau2-tau1);
                end
            end

            % For inter-pol correlations  
            % identical time index, no delay
            Cov24_XY_0tau = Cov.Cov24_XY_0tau;  
            Cov22_XY_0tau = Cov.Cov22_XY_0tau;  
            % two time groups, one delay 
            Cov22_XY = Cov.Cov22_XY;
            Cov24_XY = Cov.Cov24_XY;
            Cov222_XY = Cov.Cov222_XY;
            Cov222_XY_1tau = Cov.Cov222_XY_1tau; 

            % Construct the combined triple-energy covariance.
            Cov222_Cmb_XY = zeros(size(Cov222_XY));
            for tau1 = 1:tau_max
                for tau2 = tau1+1:tau1+tau_max
                    Cov222_Cmb_XY(tau1,tau2) = Cov222_XY(tau1,tau2)...
                        -Ex2*Cov22_XY(tau1)-Ex2*Cov22_XY(tau2)...
                        -Ex2*Cov22_X(tau2-tau1);
                end
            end

            %%% Monte-Carlo integration
            %%% Evaluated frequency points and the spectrum    

            G_SCI.corr = zeros(obj.model_Nf,1);

            for idx = 1:obj.model_Nf  


                %%% SPT: self-pol. temporal interactions
                G_SCI.corr_SPT_2D(idx) = 2*(...
                    sum(Ex2.*Cov22_X(1:tau_max).*obj.model_Kn.K22_SPT(:,idx))...
                    + sum(Cov24_X(1:tau_max).*obj.model_Kn.K24_SPT(:,idx)) );
                if obj.model_Kn.flagAppxOnly == 0
                    G_SCI.corr_SPT_3D(idx) = 2*(...
                        sum(sum(Cov222_Cmb_X.*obj.model_Kn.K222_SPT(:,:,idx))) );
                else
                    G_SCI.corr_SPT_3D(idx) = 0;
                end

                %%% XP: cross-pol. interactions
                G_SCI.corr_XP_1D(idx) = 2*( Ex2.*Cov22_XY_0tau.*obj.model_Kn.K22_XP(idx)...
                    + Cov24_XY_0tau.*obj.model_Kn.K24_XP(idx) );

                %%% XPT: self-pol. temporal interactions
                G_SCI.corr_XPT_2D(idx) = 2*(... 
                    + sum(Ex2.*Cov22_XY(1:tau_max).*obj.model_Kn.K22_XPT(:,idx))...
                    + sum(Cov24_XY(1:tau_max).*obj.model_Kn.K24_XPT(:,idx))...
                    + sum( ( Cov222_XY_1tau(1:tau_max) - Ex2.*Cov22_XY_0tau ) ...
                    .*obj.model_Kn.K222_XPT_1tau(:,idx)) );
                if obj.model_Kn.flagAppxOnly == 0
                    G_SCI.corr_XPT_3D(idx) = 2*(...
                        sum(sum(Cov222_Cmb_XY.*obj.model_Kn.K222_XPT(:,:,idx))) );
                else
                    G_SCI.corr_XPT_3D(idx) = 0;
                end

            end

            
            %%%  times 2 because of double pol effects see ref.[1] (41) (44) 
            eta_SCI.corr_SPT_2D = mean(G_SCI.corr_SPT_2D)*obj.sys_Rs/P_tot^3;
            eta_SCI.corr_SPT_3D = mean(G_SCI.corr_SPT_3D)*obj.sys_Rs/P_tot^3;
            eta_SCI.corr_XPT_2D = mean(G_SCI.corr_XPT_2D)*obj.sys_Rs/P_tot^3;
            eta_SCI.corr_XPT_3D = mean(G_SCI.corr_XPT_3D)*obj.sys_Rs/P_tot^3;
            eta_SCI.corr_XP_1D = mean(G_SCI.corr_XP_1D)*obj.sys_Rs/P_tot^3;
            
            %%%  approximate version without double delay summations
            G_SCI.corr_apprx = G_SCI.corr_SPT_2D + G_SCI.corr_XPT_2D+ G_SCI.corr_XP_1D; 
            eta_SCI.corr_apprx = mean(G_SCI.corr_apprx)*obj.sys_Rs/P_tot^3;
            
            %%%  exact version with double delay summations
            G_SCI.corr = G_SCI.corr_apprx +G_SCI.corr_SPT_3D +  G_SCI.corr_XPT_3D;
            eta_SCI.corr = mean(G_SCI.corr)*obj.sys_Rs/P_tot^3;


        end


        %% MEGN model correction helpers
        %%% Calculate kernels over spectrum summation; 

        function [kn_chi_1,kn_chi_2,kn_chi_3] = specSumCal2_1(obj,f,specBtVec,f1Mtx,f2Mtx,tau)

            Rs = obj.sys_Rs;
            %%% For two energy-dependent summation
            terms = zeros([size(f1Mtx),6]);
            terms(:,:,1) = cos(pi/Rs*tau.*(f-f2Mtx)).^2;
            terms(:,:,2) = sin(pi/Rs*tau.*(f-f2Mtx)).^2;
            terms(:,:,3) = cos(2*pi/Rs*tau.*f1Mtx);
            terms(:,:,4) = sin(2*pi/Rs*tau.*f1Mtx); 
            terms(:,:,5) = cos(2*pi/Rs*tau.*(f1Mtx+f2Mtx-f));
            terms(:,:,6) = sin(2*pi/Rs*tau.*(f1Mtx+f2Mtx-f)); 

            specSumVec2 = zeros(1,6);
 
            for idx_term = 1:6
                specSumVec2(idx_term) = abs(mean(mean(specBtVec.*terms(:,:,idx_term)))*Rs^2)^2;
            end

            %%% Typical group kn_chi_1
            coffs = [1,-1,0,0,0,0];
            kn_chi_1 = 128/81*Rs*sum(specSumVec2.*coffs); 
            %%% Typical group kn_chi_2
            coffs = [0,0,1,1,0,0];
            kn_chi_2 = 128/81*Rs*sum(specSumVec2.*coffs); 
            %%% Typical group kn_chi_3
            coffs = [0,0,0,0,1,1];
            kn_chi_3 = 128/81*Rs*sum(specSumVec2.*coffs); 
        end
 

        function [kn_xi_1,kn_psi_1] = specSumCal2_2(obj,f,f1Mtx,f2Mtx,spec1Mtx,spec2Mtx,spec_fcmb1,mu_vec1,tau)

            Rs = obj.sys_Rs;
            %%% For two energy-dependent summation
            terms = zeros([size(f1Mtx),4]);
            terms(:,:,1) = cos(pi/Rs*tau.*(f-f1Mtx)) ;
            terms(:,:,2) = sin(pi/Rs*tau.*(f-f1Mtx)) ;

            terms(:,:,3) = cos(2*pi/Rs*tau.*f2Mtx);
            terms(:,:,4) = sin(2*pi/Rs*tau.*f2Mtx);  

            kn_xi_1 = 128/81*Rs^2* ( mean(abs(spec1Mtx(:,1).*terms(:,1,1)).^2.* ...
               abs(mean(spec2Mtx.*conj(spec_fcmb1).*mu_vec1,2)*Rs).^2)*Rs ...
               - mean(abs(spec1Mtx(:,1).*terms(:,1,2)).^2.* ...
               abs(mean(spec2Mtx.*conj(spec_fcmb1).*mu_vec1,2)*Rs).^2)*Rs  );

            kn_psi_1 = 128/81*Rs^2* ( mean(abs(spec1Mtx(:,1)).^2.* ...
               abs(mean(spec2Mtx.*conj(spec_fcmb1).*mu_vec1.*terms(:,:,3),2)*Rs).^2)*Rs ...
               + mean(abs(spec1Mtx(:,1)).^2.* ...
               abs(mean(spec2Mtx.*conj(spec_fcmb1).*mu_vec1.*terms(:,:,4),2)*Rs).^2)*Rs  );
 
        end

        function [kn_psi_2] = specSumCal2_3(obj,f,f1Mtx,f2Mtx,spec1Mtx,spec2Mtx,spec_fcmb2,mu_vec2,tau)

            Rs = obj.sys_Rs;
            %%% For two energy-dependent summation
            terms = zeros([size(f1Mtx),2]); 
            terms(:,:,1) = cos(2*pi/Rs*tau.*(f1Mtx-f2Mtx+f));
            terms(:,:,2) = sin(2*pi/Rs*tau.*(f1Mtx-f2Mtx+f));   

            kn_psi_2 = 128/81*Rs^2* ( mean(abs(spec1Mtx(:,1)).^2.* ...
               abs(mean(spec2Mtx.*conj(spec_fcmb2).*mu_vec2.*terms(:,:,1),2)*Rs).^2)*Rs ...
               + mean(abs(spec1Mtx(:,1)).^2.* ...
               abs(mean(spec2Mtx.*conj(spec_fcmb2).*mu_vec2.*terms(:,:,2),2)*Rs).^2)*Rs  );
 
        end
       


        %%% Calculate spectrum summation in
        function [kn_xi_2,kn_psi_3] = specSumCal3(obj,f,specBtVec,f1Mtx,f2Mtx,tau1,tau2)

            Rs = obj.sys_Rs;

            %%% For three energy-dependent summation
            terms = zeros([size(specBtVec),15]);

            %%% kn_xi_2
            terms(:,:,2) = exp(-1j*2*pi/Rs*tau1.*(f-f1Mtx));
            terms(:,:,3) = exp(-1j*2*pi/Rs*tau2.*(f-f1Mtx));
            terms(:,:,1) = terms(:,:,2)+terms(:,:,3);

            terms(:,:,5) = exp(-1j*2*pi/Rs*tau1.*(f1Mtx-f));
            terms(:,:,6) = exp(-1j*2*pi/Rs*(tau2-tau1).*(f-f1Mtx));
            terms(:,:,4) = terms(:,:,5)+terms(:,:,6);
            
            terms(:,:,8) = exp(-1j*2*pi/Rs*tau2.*(f1Mtx-f));
            terms(:,:,9) = exp(-1j*2*pi/Rs*(tau1-tau2).*(f-f1Mtx));
            terms(:,:,7) = terms(:,:,8)+terms(:,:,9);

            %%% kn_psi_3
            terms(:,:,10) = exp(-1j*2*pi/Rs*(tau1.*f1Mtx+tau2.*f2Mtx));
            terms(:,:,11) = exp(-1j*2*pi/Rs*(tau2.*f1Mtx+tau1.*f2Mtx));
            terms(:,:,12) = exp(-1j*2*pi/Rs*(-tau1.*f1Mtx+(tau2-tau1).*f2Mtx));
            terms(:,:,13) = exp(-1j*2*pi/Rs*((tau2-tau1).*f1Mtx-tau1.*f2Mtx));
            terms(:,:,14) = exp(-1j*2*pi/Rs*((tau1-tau2).*f1Mtx-tau2.*f2Mtx));
            terms(:,:,15) = exp(-1j*2*pi/Rs*(-tau2.*f1Mtx+(tau1-tau2).*f2Mtx));

            specSumVec3 = zeros(1,15);
            for idx_term = 1:15
                specSumVec3(idx_term) = abs(mean(mean(specBtVec.*terms(:,:,idx_term)))*Rs^2)^2;
            end

            %%% Typical group L3a
            coffs =  [1,-1,-1,1,-1,-1,1,-1,-1,0,0,0,0,0,0];
            kn_xi_2 = 64/81*Rs*sum(specSumVec3.*coffs);  
            %%% Typical group L3b
            coffs = [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1];
            kn_psi_3 = 64/81*Rs*sum(specSumVec3.*coffs); 


        end

        function [pulseSpec] = RC_spectra(obj,freq_vec)
            T = obj.sys_Ts;
            rollOff = obj.sys_rollOff;

            sz = size(freq_vec);
            freq_vec = freq_vec(:);
            pulseSpec = T.*zeros(length(freq_vec),1);

            bound1 = (1-rollOff)/2/T;
            bound2 = (1+rollOff)/2/T;

            range1 = abs(freq_vec)<=bound1;
            pulseSpec(range1) =1;

            range2 = abs(freq_vec)>bound1 & abs(freq_vec)<bound2;
            pulseSpec(range2) = 1/2*(1+...
                cos(pi*T/rollOff*(abs(freq_vec(range2))-(1-rollOff)/2/T)));

            pulseSpec = reshape(pulseSpec,sz);
        end


        function [pulseSpec] = RRC_spectra(obj,freq_vec)
            [pulseSpecRC] = RC_spectra(obj,freq_vec);
            pulseSpec = sqrt(pulseSpecRC);
            %%% RRC spectrum top value T = 1/Rs
            pulseSpec = pulseSpec./max(pulseSpec(:))*obj.sys_Ts; 
        end

        %%% Calculate link function mu
        function mu = linkFunc(obj, ftms)
            alpha = obj.sys_alpha;
            beta2 = obj.sys_beta2;
            gamma = obj.sys_gamma;
            Ls    = obj.sys_Ls;
            Ns    = obj.sys_Ns;

            sz = size(ftms);
            %%% ftms = (f1-f).*(f2-f);
            %%% -Rs < f1-f, f2-f < Rs
            %%% -Rs^2 < f_cmb < Rs^2

            %%% the efficiency of non-degenerate FWM occurring among three spectral
            %%% components of the signal placed at frequencies f1, f2, f3=(f1+f2-f)
            zeta = gamma*(1-exp(-2*alpha*Ls+1j*4*pi^2*beta2*Ls.*ftms))./...
                (2*alpha-1j*4*pi^2*beta2.*ftms);

            %%% The factor nu relates to the coherent interference of the NLI field
            %%% contributions produced in different spans, when they are summed up at
            %%% the receiver location.
            nu = sin(2*beta2*pi^2*Ls*Ns.*ftms)./...
                sin(2*beta2*pi^2*Ls.*ftms).*...
                exp(1j*2*beta2*pi^2*Ls*(Ns-1).*ftms);

            %%% Avoid nan values
            nu(ftms==0) = Ns.* exp(1j*2*beta2*pi^2*Ls*(Ns-1).*ftms(ftms==0));

            %%% the link function weighs the generation of NLI and depends only
            %%% on the link parameters, but not on the launched signal
            mu = zeta.*nu;
            mu = reshape(mu,sz);

        end


    end
end



 



 
