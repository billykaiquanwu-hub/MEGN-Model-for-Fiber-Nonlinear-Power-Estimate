function [Cov] = calculateSymbolCorrelations(Ea2,Ea4,Ea6,n,mapping,tauMax)
%CALCULATESYMBOLCORRELATIONS Analytical CCDM symbol correlations.

corrLen = round(n/mapping);

[Rc,Ex2,Ex4] = CCDMSymCorrAna_WithinBlock(Ea2,Ea4,Ea6,n,mapping);
Cov.Ex2 = Ex2;

Cov.Cov22_X = linspace((corrLen-1)/corrLen*Rc.Rc22_X+(1/corrLen-1)*Ex2*Ex2,0,corrLen).';
if corrLen < 2*tauMax
    Cov.Cov22_X = [Cov.Cov22_X;zeros(2*tauMax-corrLen,1)];
else
    Cov.Cov22_X = Cov.Cov22_X(1:2*tauMax);
end

Cov.Cov24_X = linspace((corrLen-1)/corrLen*Rc.Rc24_X+(1/corrLen-1)*Ex2*Ex4,0,corrLen).';
if corrLen < 2*tauMax
    Cov.Cov24_X = [Cov.Cov24_X;zeros(2*tauMax-corrLen,1)];
else
    Cov.Cov24_X = Cov.Cov24_X(1:2*tauMax);
end
Cov.CovPar24_X = 0;

if mapping == 4
    Cov.Cov22_XY_0tau = Rc.Rc22_XY - Ex2*Ex2;
    Cov.Cov24_XY_0tau = Rc.Rc24_XY - Ex2*Ex4;
else
    Cov.Cov22_XY_0tau = 0;
    Cov.Cov24_XY_0tau = 0;
end

if mapping == 4
    Cov.Cov24_XY = linspace((corrLen-1)/corrLen*Rc.Rc24_XY+(1/corrLen-1)*Ex2*Ex4,0,corrLen).';
    Cov.Cov22_XY = linspace((corrLen-1)/corrLen*Rc.Rc22_XY+(1/corrLen-1)*Ex2*Ex2,0,corrLen).';

    if corrLen < 2*tauMax
        Cov.Cov24_XY = [Cov.Cov24_XY;zeros(2*tauMax-corrLen,1)];
        Cov.Cov22_XY = [Cov.Cov22_XY;zeros(2*tauMax-corrLen,1)];
    else
        Cov.Cov24_XY = Cov.Cov24_XY(1:2*tauMax);
        Cov.Cov22_XY = Cov.Cov22_XY(1:2*tauMax);
    end
else
    Cov.Cov24_XY = zeros(size(Cov.Cov24_X));
    Cov.Cov22_XY = zeros(size(Cov.Cov22_X));
end

Cov.CovPar24_XY = 0;
Cov.CovPar222_X = zeros(tauMax,2*tauMax);
Cov.CovPar222_XY = zeros(tauMax,2*tauMax);

Cov.Cov222_X = zeros(tauMax,2*tauMax);
for tau1 = 1:tauMax
    for tau2 = tau1+1:tau1+tauMax
        if tau2 < corrLen
            Cov.Cov222_X(tau1,tau2) = (corrLen-tau2)/corrLen*Rc.Rc222_X + tau2/corrLen*Rc.Rc22_X*Ex2-Ex2^3;
        elseif tau2>= corrLen && tau1 < corrLen && tau2-tau1>=corrLen
            Cov.Cov222_X(tau1,tau2) = (corrLen-tau1)/corrLen*Rc.Rc22_X*Ex2 + tau1/corrLen*Ex2^3-Ex2^3;
        elseif tau2>= corrLen && tau1 < corrLen && tau2-tau1<corrLen
            Cov.Cov222_X(tau1,tau2) = (2*corrLen-tau2)/corrLen*Rc.Rc22_X*Ex2 + (tau2-corrLen)/corrLen*Ex2^3-Ex2^3;
        elseif tau1>= corrLen && tau2 - tau1 < corrLen
            Cov.Cov222_X(tau1,tau2) = (corrLen-tau2+tau1)/corrLen*Rc.Rc22_X*Ex2 + (tau2-tau1)/corrLen*Ex2^3-Ex2^3;
        else
            Cov.Cov222_X(tau1,tau2) = 0;
        end
    end
end

Cov.Cov222_XY_1tau = linspace((corrLen-1)/corrLen*Rc.Rc222_XY+(1/corrLen-1)*Ex2^3,Rc.Rc22_XY*Ex2-Ex2^3,corrLen).';
if corrLen < 2*tauMax
    Cov.Cov222_XY_1tau = [Cov.Cov222_XY_1tau;repmat(Rc.Rc22_XY*Ex2-Ex2^3,2*tauMax-corrLen,1)];
else
    Cov.Cov222_XY_1tau = Cov.Cov222_XY_1tau(1:2*tauMax);
end

Cov.Cov222_XY = zeros(tauMax,2*tauMax);
if mapping == 4
    for tau1 = 1:tauMax
        for tau2 = tau1+1:tau1+tauMax
            if tau2 < corrLen
                Cov.Cov222_XY(tau1,tau2) = (corrLen-tau2)/corrLen*Rc.Rc222_XY + tau2/corrLen*Rc.Rc22_XY*Ex2-Ex2^3;
            elseif tau2>= corrLen && tau1 < corrLen && tau2-tau1>=corrLen
                Cov.Cov222_XY(tau1,tau2) = (corrLen-tau1)/corrLen*Rc.Rc22_XY*Ex2 + tau1/corrLen*Ex2^3-Ex2^3;
            elseif tau2>= corrLen && tau1 < corrLen && tau2-tau1<corrLen
                Cov.Cov222_XY(tau1,tau2) = (2*corrLen-tau2)/corrLen*Rc.Rc22_XY*Ex2 + (tau2-corrLen)/corrLen*Ex2^3-Ex2^3;
            elseif tau1>= corrLen && tau2 - tau1 < corrLen
                Cov.Cov222_XY(tau1,tau2) = (corrLen-tau2+tau1)/corrLen*Rc.Rc22_X*Ex2 + (tau2-tau1)/corrLen*Ex2^3-Ex2^3;
            else
                Cov.Cov222_XY(tau1,tau2) = 0;
            end
        end
    end
else
    for tau1 = 1:tauMax
        for tau2 = tau1+1:tau1+tauMax
            if tau2 - tau1 < corrLen
                Cov.Cov222_XY(tau1,tau2) = (corrLen-(tau2-tau1))/corrLen*Rc.Rc22_X*Ex2 + (tau2-tau1)/corrLen*Ex2^3 - Ex2^3;
            end
        end
    end
end
end

function [Rc,Ex2,Ex4] = CCDMSymCorrAna_WithinBlock(Ea2,Ea4,Ea6,n,mapping)
Ea2a2 = (n*Ea2^2-Ea4)/(n-1);
Ea2a4  = (n*Ea2*Ea4-Ea6)/(n-1);
Ea2a2a2  = n^2/(n-1)/(n-2)*Ea2^3+2*n/(n-1)^2/(n-2)*Ea6-...
    n*(3*n-1)/(n-1)^2/(n-2)*Ea2*Ea4;

Ex2 = 2*Ea2;

if mapping == 1
    Rc.Rc22_X = 2*Ea2a2+2*Ea2^2;
    Rc.Rc24_X = 2*Ea2a4+2*Ea2*Ea4+4*Ea2a2*Ea2;
    Rc.Rc222_X = 2*Ea2a2a2+6*Ea2a2*Ea2;
    Ex4 = 2*Ea4+2*Ea2^2;
elseif mapping == 2 || mapping == 4
    Rc.Rc22_X = 4*Ea2a2;
    Rc.Rc24_X = 4*Ea2a4+4*Ea2a2a2;
    Rc.Rc222_X = 8*Ea2a2a2;
    Ex4 = 2*Ea4+2*Ea2a2;
end

if mapping == 1
    Rc.Rc22_XY = 4*Ea2^2;
    Rc.Rc24_XY = 4*Ea2*Ea4+4*Ea2^3;
    Rc.Rc222_XY = 4*Ea2a2*Ea2+4*Ea2^3;
elseif mapping == 2
    Rc.Rc22_XY = 4*Ea2^2;
    Rc.Rc24_XY = 4*Ea2*Ea4+4*Ea2a2*Ea2;
    Rc.Rc222_XY = 8*Ea2a2*Ea2;
elseif mapping == 4
    Rc.Rc22_XY = 4*Ea2a2;
    Rc.Rc24_XY = 4*Ea2a4+4*Ea2a2a2;
    Rc.Rc222_XY = 8*Ea2a2a2;
end
end
