clear
fs = 44100;
lengthSound = 10*fs;
k=1/fs; %time step

Bow = false;
plateTension = false;

stringConn = true;
numStrings = 6;

tubeConn = true;

%% Define plate parameters

Lx = 0.5; %side length (x)
Ly = 0.7; %side length (y)
sigma0P = 0.0005; %frequency independent damping
sigma1P = 0.02; %frequency dependent damping
rhoP = 7700; %material density
HP = 7 * 10^(-3); %Thickness
EP = 2.5 * 10^(11); %Young's modulus
nuP = 0.3; %Poisson's ratio
T60 = 0;
D = EP*HP^3/(12*(1-nuP^2)); % stifness coefficient
kappaP = sqrt(D/(rhoP*HP)); % stifness paramater
if plateTension == true
    TP = 200000;
    %gammaP= sqrt(TP/(rhoP*H));
    cPSq = TP/(rhoP*HP);
    hP=4*(2*sqrt((cPSq*k^2+4*sigma1P*k+sqrt(cPSq*k^2+4*sigma1P*k)^2+4*kappaP^2*k^2)/2));
else
    hP = 2*sqrt(k*(sigma1P+sqrt(kappaP^2+sigma1P^2)));
end
Nx = floor(Lx/hP); %grid steps (x)
Ny = floor(Ly/hP); %grid steps (y)
hP = Lx /Nx; %redefine grid spacing (use the smallest dimension)
muP = kappaP*k/hP^2;
S = 2 * sigma1P * k /hP^2;
if plateTension == true
    lambdaPSq = cPSq * k * k /(hP*hP);
    dP= 1/(1+sigma0P*k);
end

%% Define string parameters

LS = 0.5;                      % Length [m]
rhoS = 7850;                    % Material density [kg/m^3]
rS = 0.001;                    % Radius [m]
AS = pi * rS^2;                 % Cross-sectional area [m^2]
ES = 2e11;                      % Young's modulus [Pa]
IS = pi * rS^4 / 4;             % Area moment of inertia [m^4]

sigma0S = 0.2;                    % Frequency-independent damping [1/s]
sigma1S = 0.005;                % Frequency-dependent damping [m^2/s]
kappaSSq = ES * IS / (rhoS * AS);   % Stiffness term
As = 1 + sigma0S * k;

TS = zeros(numStrings,1);
cSSq = zeros(numStrings,1);
%hS = zeros(numStrings,1);
NS = zeros(numStrings,1);
hS = zeros(numStrings,1);
lambdaSSq = zeros(numStrings,1);
muSSq = zeros(numStrings,1);
uS1  = zeros(numStrings,1);
uS2 = zeros(numStrings,1);
uS3 = zeros(numStrings,1);
stringConnTerm = zeros(numStrings,1);
for i = 1 : numStrings
    TS(i) = 1200*(i/numStrings);               % Tension [N]
    cSSq(i) = TS(i) / (rhoS * AS);        % Wave speed 
    hS(i) = sqrt((cSSq(i) * k^2 + 4 * sigma1S * k ...
        + sqrt((cSSq(i) * k^2 + 4 * sigma1S * k)^2 + 16 * kappaSSq * k^2))/2);
    NS(i) = floor(LS / hS(i));
    hS(i) = LS / NS(i);
    lambdaSSq(i) = cSSq(i)*k*k/hS(i)^2;
    muSSq(i) = kappaSSq*k*k/hS(i)^4;
    uS1(i)  = 2-2*lambdaSSq(i)-6*muSSq(i)-4*sigma1S*k/hS(i)^2;
    uS2(i) = lambdaSSq(i)+4*muSSq(i)+2*sigma1S*k/hS(i)^2;
    uS3(i) = -1+sigma0S*k + 4*sigma1S*k/hS(i)^2;
    stringConnTerm(i) = k * k / (rhoS * AS * hS(i) * (1.0 + sigma0S * k));
end

%% Define bar parameters
LB = 0.5;
rhoB = 7500;
rB = 3;
AB = pi * rB^2;                 % Cross-sectional area [m^2]
EB = 2e11;                      % Young's modulus [Pa]
IB = pi * rB^4 / 4;             % Area moment of inertia [m^4]
sigma0B = 0.2;                    % Frequency-independent damping [1/s]
sigma1B = 0.005;                % Frequency-dependent damping [m^2/s]

kappaBSq = EB * IB / (rhoB * AB);   % Stiffness term
hB = sqrt(2*k*(sigma1B^2+sqrt(kappaBSq+sigma1B^2)));
NB = floor(LB / hB);
hB = LB / NS;

%% Define tube parameters
cylinderLength = 1;
bellLength = 0.2;
cylinderRadius = 0.01;
bellEndRadius = 0.1;
LT = cylinderLength+bellLength;

cT = 343;
h = cT* k;

N_C = floor(cylinderLength / h);
N_B = floor(bellLength / h);
NT = N_B + N_C;
 
bellCurve = 1;
S_C = zeros(N_C+1,1);
S_C(:) = pi*cylinderRadius^2;
S_B = zeros(N_B,1);
if bellCurve == 1
    rT = cylinderRadius;
    rGrowth = (bellEndRadius-cylinderRadius)/N_B;
    for i = 1 : N_B
        rT = rT + rGrowth;
        S_B(i) = pi*rT^2;
    end
end

if bellCurve == 2;
    rT = cylinderRadius;
    rGrowth = exp(log(bellEndRadius/cylinderRadius)/N_B);
    for i = 1 : N_B
        rT = cylinderRadius*rGrowth^i;
        S_B(i) = pi*rT^2;
    end
end

if bellCurve == 3;
    rT = cylinderRadius;
    rGrowth = (-cylinderRadius+bellEndRadius)/log(N_B);
    for i = 1 : N_B
        rT = cylinderRadius+rGrowth*log(i);
        S_B(i) = pi*rT^2;
    end
end


ST = cat(1,S_C,S_B);


hT= LT/NT;
lambdaT = cT*k/hT;
rhoT = 1.225;

pTNext = zeros(NT+1,1);
pT= pTNext;
pTPrev = pTNext;

vTNext = zeros(NT+1,1);
vT = vTNext;

lcT = 1;

R_1 = rhoT*cT;
R_2 = 0.505* rhoT*cT;
L_r = 0.613*rhoT*sqrt(ST(NT+1)/pi);
C_r = 1.111* sqrt(ST(NT+1))/(rhoT*cT^2*sqrt(pi));

zeta_1=(2*R_2*k)/(2*R_1*R_2*C_r+k*(R_1+R_2));
zeta_2 = (2*R_1*R_2*C_r-k*(R_1+R_2))/(2*R_1*R_2*C_r+k*(R_1+R_2));
zeta_3 = k/(2*L_r)+zeta_1/(2*R_2)+(C_r*zeta_1)/k;
zeta_4 = (zeta_2+1)/(2*R_2)+(C_r*zeta_2-C_r)/k;

v_int=0;
p_int=0;

%% Connection parameters
K1 = 1e3;
K3 = 1e10;
R = 0.01;

%K1 = 1e8;
%K3 = 0;
%R = 0.01;

% Connection locations
connSposRatio = zeros(numStrings, 1);
connXposRatio = connSposRatio;
connYposRatio = connSposRatio;
connSpos = connSposRatio;
connXpos = connSposRatio;
connYpos = connSposRatio;
connSpos2 = connSposRatio;
connXpos2 = connSposRatio;
connYpos2 = connSposRatio;

for i = 1 : numStrings
    connSposRatio(i) = 0.1;    % along string
    connXposRatio(i) = ((i-1)/numStrings)+0.2;   % along x-direction plate
    connYposRatio(i) = 0.25;   % along y-direction plate
    connSpos(i) = connSposRatio(i) * LS;   % along string
    connXpos(i) = connXposRatio(i) * Lx;   % along x-direction plate
    connYpos(i) = connYposRatio(i) * Ly;   % along y-direction plate
    connSpos2(i) = LS - connSpos(i);
    connXpos2(i) = connXpos(i);
    connYpos2(i) = connYpos(i)+LS-connSpos(i)*2;
end
% connection force
connF = 0;

%% Position parameters
excXposRatio = 0.2;
excYposRatio = linspace(0.2,0.9,lengthSound+1);
excXpos = excXposRatio*Lx;
excYpos = excYposRatio*Ly;
outXposRatio = 0.5;
outYposRatio = 0.5;
outXpos = outXposRatio*Lx;
outYpos = outYposRatio*Ly;
hx = Lx/Nx;
hy = Ly/Ny;

li = floor(excXpos/hx);
mi = floor(excYpos(1)/hy);
alphaX = excXpos/hx-li;
alphaY = excYpos(1)/hy-mi;

%% Bow parameters
attack = linspace(0,1,0.1*lengthSound);
sustain= linspace(1,1,0.3*lengthSound);
release = linspace(1, 0, 0.1*lengthSound);
env = cat(2, attack,sustain,release, zeros(1,0.5*lengthSound),[0]);
vB = env.*0.1;
FB = env.*10;
a = 0.000000000001;
vRel = 0;
vRelPrev = 0;
%phiVRel = 0;

tol = 0.0001;

%% Raised cosine excitation parameters
maxForce = 10; %maximum force
excTime = 0.001;
t0 = 0;
t = 0;
%malletForce = maxForce/2*(1-cos((2*pi*(t-t0))/excTime)); 

%%
J=0;
uP = zeros(Nx+1,Ny+1);
uPNext = uP;
uPPrev = uP;
out = zeros(lengthSound+1,1);
uPStar = 0;

%%
NSMax = max(NS);
for nS = 1 : numStrings
    uS = zeros(NSMax+1,numStrings);
    uSNext = uS;
    uSPrev = uS;
end
%uSStar = 0;
%%  discrete connection position and alpha for interpolation and spreading
if tubeConn == true
    numStrings = 1;
    stringConn = false;
end
lcS = zeros(numStrings, 1);
alphaConnS = lcS;
lcP = lcS;
mcP = lcS;
alphaConnX = lcS;
alphaConnY = lcS;
lcS2 = lcS;
alphaConnS2 = lcS;
lcP2 = lcS;
mcP2 = lcS;
alphaConnX2 = lcS;
alphaConnY2 = lcS;
for i = 1 : numStrings
    lcS(i) = floor(connSpos(i)/hS(i));
    alphaConnS(i) = connSpos(i)/hS(i)-lcS(i);
    lcP(i) = floor(connXpos(i)/hx); 
    mcP(i) = floor(connYpos(i)/hy);
    alphaConnX(i) = connXpos(i)/hx-lcP(i);
    alphaConnY(i) = connYpos(i)/hy-mcP(i);
    lcS2(i) = floor(connSpos2(i)/hS(i));
    alphaConnS2(i) = connSpos2(i)/hS(i)-lcS2(i);
    lcP2(i) = floor(connXpos2(i)/hx); 
    mcP2(i) = floor(connYpos2(i)/hy);
    alphaConnX2(i) = connXpos2(i)/hx-lcP2(i);
    alphaConnY2(i) = connYpos2(i)/hy-mcP2(i);
end
eta = zeros(numStrings,1);
etaPrev = zeros(numStrings,1);
etaNext = zeros(numStrings,1);
rPlus = zeros(numStrings,1);
rMinus = zeros(numStrings,1);
connF = zeros(numStrings,1);
connFTot = 0;

plateConnTerm = k^2 / (rhoP * HP * hP^2 * (1.0 + sigma0P * k));

tubeConnTerm = k^2 / (rhoT * ST(1)*hT);
lcT = 1;
%%
for n = 1 : lengthSound+1
    if Bow == true
        li = floor(excXpos/hx);
        mi = floor(excYpos(n)/hy);
        alphaX = excXpos/hx-li;
        alphaY = excYpos(n)/hy-mi;
         %b = ((2/k) + 2*sigma0)*vB - ((2/k^2)*u[li, mi]-uPrev[li, mi])...
         %+ (kappa^2/h^4)*(u(li+2,mi)+u(li-2,mi)+u(li,mi+2)+u(li,mi-2)+2*(u[li+1, mi+1]+u[li+1, mi-1]+u[li-1, mi+1]+u[li-1, mi-1])-8*(u[li+1, mi]+u[li-1, mi]+u[li, mi+1]+u[li, mi-1])+20*u[li, mi])...
         %- 2*sigma1/(k*h^2)*(u[li+1, mi]+u[li-1, mi]+u[li, mi+1]+u[li, mi-1]-u[li+1, mi]-u[li-1, mi]-u[li, mi+1]-u[li, mi-1]-4*(u(li,mi)-uPrev(li,mi));
        b = ((2/k) + 2*sigma0P)*vB(n) - ((2/k^2)*interpolation2D(uP,li,mi,alphaX, alphaY)-interpolation2D(uPPrev,li,mi,alphaX, alphaY))...
         + (kappaP^2/hP^4)*(interpolation2D(uP,li+2,mi,alphaX, alphaY)+interpolation2D(uP,li-2,mi,alphaX, alphaY)+interpolation2D(uP,li,mi+2,alphaX, alphaY)+interpolation2D(uP,li,mi-2,alphaX, alphaY))...
         + 2 * (interpolation2D(uP,li+1,mi+1,alphaX, alphaY)+interpolation2D(uP,li+1,mi-1,alphaX, alphaY)+interpolation2D(uP,li-1,mi+1,alphaX, alphaY)+interpolation2D(uP,li-1,mi-1,alphaX, alphaY)-8*(interpolation2D(uP,li+1,mi,alphaX, alphaY)+interpolation2D(uP,li-1,mi,alphaX, alphaY)+interpolation2D(uP,li,mi+1,alphaX, alphaY)+interpolation2D(uP,li,mi-1,alphaX, alphaY))+20*interpolation2D(uP,li,mi,alphaX, alphaY))...
         - 2*sigma1P/(k*hP^2)*(interpolation2D(uP,li+1,mi,alphaX, alphaY)+interpolation2D(uP,li-1,mi,alphaX, alphaY)+interpolation2D(uP,li,mi+1,alphaX, alphaY)+interpolation2D(uP,li,mi-1,alphaX, alphaY)-interpolation2D(uP,li+1,mi,alphaX, alphaY)-interpolation2D(uP,li-1,mi,alphaX, alphaY)-interpolation2D(uP,li,mi+1,alphaX, alphaY)-interpolation2D(uP,li,mi-1,alphaX, alphaY)-4*(interpolation2D(uP,li,mi,alphaX, alphaY)-interpolation2D(uPPrev,li,mi,alphaX, alphaY)));
        eps = 1;
        i = 0 ;
        while eps > tol
           %phiVRel = sqrt(2*a)*vRel*exp(-a*vRel^2+0.5);
           %phiVRelalt = sqrt(2*a)*(1-2*a(vRel)^2)*exp(-a*(vRel^2+1/2));
           vRel = vRelPrev - (((2/k+2*sigma0P)*vRelPrev+FB(n)*sqrt(2*a)*vRelPrev*exp(-a*vRelPrev^2+0.5)+b)/(2/k + 2*sigma0P+FB(n)*sqrt(2*a)*(1-2*a*(vRel)^2)*exp(-a*(vRel^2+0.5))));
           eps=  abs(vRel-vRelPrev);
           vRelPrev = vRel;
           i = i+1;
           if 1000 < i
              break
           end
        end
        exc = sqrt(2*a)*vRel*exp(-a*vRel^2+0.5)*FB(n);
    else
        if n-1 < floor(excTime*fs)
        malletForce = maxForce/2*(1-cos((2*pi*(t-t0))/(excTime)));
        t = t + k;
        else
        malletForce = 0;
        end
        exc = -malletForce;
    end
    for l = 3 : Nx-1
        for m = 3 : Ny-1
            if l == li && m == mi
                J = ((1-alphaX)*(1-alphaY))/(hx*hy);
            elseif l == li && m == mi+1
                J = ((1-alphaX)*alphaY)/(hx*hy);
            elseif l == li+1 && m == mi
                J = (alphaX*(1-alphaY))/(hx*hy);
            elseif l == li+1 && m == mi+1
                J = (alphaX*alphaY)/(hx*hy);
            else
                J = 0;
            end
            if plateTension == true
                uPNext(l,m)=(2-20*muP^2-4*S-4*lambdaPSq)*dP*uP(l,m)...
                +(8*muP^2+S+lambdaPSq)*dP*(uP(l+1,m)+uP(l-1,m)+uP(l,m+1)+uP(l,m-1))...
                -2*muP^2*dP*(uP(l+1,m+1)+uP(l-1,m+1)+uP(l+1,m-1)+uP(l-1,m-1))...
                -muP^2*dP*(uP(l+2,m)+uP(l-2,m)+uP(l,m+2)+uP(l,m-2))...
                +(sigma0P*k-1+4*S)*dP*uPPrev(l,m)...
                -S*dP*(uPPrev(l+1,m)+uPPrev(l-1,m)+uPPrev(l,m+1)+uPPrev(l,m-1))...
                -J*exc;
            else
                uPNext(l,m)=(2-20*muP^2-4*S)*uP(l,m)...
                    +(8*muP^2+S)*(uP(l+1,m)+uP(l-1,m)+uP(l,m+1)+uP(l,m-1))...
                    -2*muP^2*(uP(l+1,m+1)+uP(l-1,m+1)+uP(l+1,m-1)+uP(l-1,m-1))...
                    -muP^2*(uP(l+2,m)+uP(l-2,m)+uP(l,m+2)+uP(l,m-2))...
                    +(sigma0P*k-1+4*S)*uPPrev(l,m)...
                    -S*(uPPrev(l+1,m)+uPPrev(l-1,m)+uPPrev(l,m+1)+uPPrev(l,m-1))...
                    -J*exc;
            end
                %+J*FB(n)*phiVRel;
        end
    end

    if tubeConn == true
        % Iterate over interleaved spacial grid
        for l = 1:NT
            % Update equation for velocity
            vTNext(l)= vT(l)-((lambdaT/(rhoT*cT))*(pT(l+1)-pT(l)));
        end

        % Iterate over spacial grid
        for l = 2:NT
            % Update interleaved cross section
            S_minus = 0.5*(ST(l)+ST(l-1));
            S_plus = 0.5*(ST(l)+ST(l+1));
            % Update equation for pressure
            pTNext(l) = pT(l)-((rhoT*cT*lambdaT)/((S_plus+S_minus)/2))*(vTNext(l)*S_plus-vTNext(l-1)*S_minus);
        end

        % Update right boundary
        pTNext(NT+1) = (1-rhoT*cT*lambdaT*zeta_3)/(1+rhoT*cT*lambdaT*zeta_3)*pT(NT+1)-((2*rhoT*cT*lambdaT)/(1+rhoT*cT*lambdaT*zeta_3))*(v_int+zeta_4*p_int-(0.5*(ST(NT+1)+ST(NT))*vTNext(NT)/ST(NT+1)));
    
        % Update the internal states
        v_int=v_int+k/L_r*0.5*(pTNext(NT+1)+pT(NT+1));
        p_int= zeta_1 * 0.5*(pTNext(NT+1)+pT(NT+1)) + zeta_2 * p_int;

        etaNext = pTNext(lcT) - uPNext(lcP, mcP);
        eta = pT(lcT) - uP(lcP, mcP);
        etaPrev = pTPrev(lcT) - uPPrev(lcP, mcP);
        rPlus = 0.5 * K1 + 0.5 * K3 * eta^2 + 0.5 * fs * R;
        rMinus = 0.5 * K1 + 0.5 * K3 * eta^2 - 0.5 * fs * R;

        connF = (etaNext  + rMinus / rPlus * etaPrev) / (1/rPlus + plateConnTerm + tubeConnTerm);

        pTNext(lcT) = pTNext(lcT) - connF*tubeConnTerm;
        uPNext(lcP,mcP) = uPNext(lcP,mcP)+connF*plateConnTerm;

        uPPrev = uP;
        uP = uPNext;

        pTPrev = pT;
        pT = pTNext;
        vT = vTNext;
        out(n)= pT(NT);
    end
    if stringConn == true
        for nS = 1 : numStrings
            for l = 3 : NS(nS)-2
                uSNext(l,nS) = (uS1(nS) * uS(l, nS) + uS2(nS) * (uS(l+1, nS)+uS(l-1, nS)) - muSSq(nS) * (uS(l+2, nS)+uS(l-2, nS)) + uS3(nS) * uSPrev(l, nS)-2*sigma1S*k/hS(nS)^2*(uSPrev(l+1, nS)+uSPrev(l-1, nS)))/As;
            end
        end
        for nS = 1 : numStrings
            etaNext(nS) = uSNext(lcS(nS), nS) - uPNext(lcP(nS), mcP(nS));
            eta(nS) = uS(lcS(nS), nS) - uP(lcP(nS), mcP(nS));
            etaPrev(nS) = uSPrev(lcS(nS), nS) - uPPrev(lcP(nS), mcP(nS));
            rPlus(nS) = 0.5 * K1 + 0.5 * K3 * eta(nS)^2 + 0.5 * fs * R;
            rMinus(nS) = 0.5 * K1 + 0.5 * K3 * eta(nS)^2 - 0.5 * fs * R;
            connF(nS) = (etaNext(nS)  + rMinus(nS) / rPlus(nS) * etaPrev(nS)) / (1/rPlus(nS) + plateConnTerm + stringConnTerm(nS));
            connFTot = connFTot + connF(nS);
            uSNext(lcS(nS), nS) = uSNext(lcS(nS), nS) - connF(nS)*stringConnTerm(nS);
            uPNext(lcP(nS),mcP(nS)) = uPNext(lcP(nS),mcP(nS))+connF(nS)*plateConnTerm;
    
            etaNext(nS) = uSNext(lcS2(nS), nS) - uPNext(lcP2(nS), mcP2(nS));
            eta(nS) = uS(lcS2(nS), nS) - uP(lcP2(nS), mcP2(nS));
            etaPrev(nS) = uSPrev(lcS2(nS), nS) - uPPrev(lcP2(nS), mcP2(nS));
            rPlus(nS) = 0.5 * K1 + 0.5 * K3 * eta(nS)^2 + 0.5 * fs * R;
            rMinus(nS) = 0.5 * K1 + 0.5 * K3 * eta(nS)^2 - 0.5 * fs * R;
            connF(nS) = (etaNext(nS)  + rMinus(nS) / rPlus(nS) * etaPrev(nS)) / (1/rPlus(nS) + plateConnTerm + stringConnTerm(nS));
            connFTot = connFTot + connF(nS);
            uSNext(lcS2(nS), nS) = uSNext(lcS2(nS), nS) - connF(nS)*stringConnTerm(nS);
            uPNext(lcP2(nS),mcP2(nS)) = uPNext(lcP2(nS),mcP2(nS))+connF(nS)*plateConnTerm;
        end
        uSOut = 0;
        for nS = 1:numStrings
            uSOut = uSOut + uS(12,nS);
        end
        out(n)= uP(12,12)+uSOut;
        %out(n)= uS(12);
        uPPrev = uP;
        uP = uPNext;

        uSPrev = uS;
        uS = uSNext;    
    end

    
    %etaNext = interpolation1D(uSNext, lcS, alphaConnS) - interpolation2D(uPNext, lcP, mcP, alphaConnX, alphaConnY);
    %eta = interpolation1D(uS, lcS, alphaConnS) - interpolation2D(uP, lcP, mcP, alphaConnX, alphaConnY);
    %etaPrev = interpolation1D(uSPrev, lcS, alphaConnS) - interpolation2D(uPPrev, lcP, mcP, alphaConnX, alphaConnY);


    %% Calculate force for conection

    %rPlus = K1 / 4 + K3 * eta^2 / 2 + R / (2*k);
    %rMinus = K1 / 4 + K3 * eta^2 / 2 - R / (2*k);


    %uSStar = interpolation1D(uSNext, lcS, alphaConnS);
    %uPStar = interpolation2D(uPNext, lcP, mcP, alphaConnX, alphaConnY);
    %connF = (uSStar - uPStar + K1 / (2 * rPlus) * eta + rMinus / rPlus * etaPrev) ...
    %        / (1/rPlus + Iu * Ju * k^2 / (rhoS * A * (1+sig0S * k)) + Iw * Jw * k^2 / (rhoP * HP * (1+sig0P * k)));

    %connF = (uSStar - uPStar + K1 / (2 * rPlus) * eta + rMinus / rPlus * etaPrev) ...
    %        / (1/rPlus + interpolation1D(k^2 , lcS, alphaConnS) * Ju * k^2 / (rhoS * A * (1+sig0S * k)) + interpolation2D(k^2, lcP, mcP, alphaConnX, alphaConnY)) * Jw * k^2 / (rhoP * HP * (1+sig0P * k)));
    %vRel=(1-alphaX)*(1-alphaY)*(1/(2*k))*(uNext(li, mi)-uPrev(li, mi)) ...
    %    + ((1-alphaX)*alphaY)*(1/(2*k))*(uNext(li, mi+1)-uPrev(li, mi+1)) ...
    %    + (alphaX*(1-alphaY))*(1/(2*k))*(uNext(li+1, mi)-uPrev(li+1, mi)) ...
    %    + (alphaX*alphaY)*(1/(2*k))*(uNext(li+1, mi+1)-uPrev(li+1, mi+1)) ...
    %    - vB(n);
    %vRel = (1/(2*k))*(uNext(li, mi)-uPrev(li, mi)) - vB(n);
    %phiVRel = sqrt(2*a)*vRel*exp(-a*vRel^2+0.5); 

    %mesh(uP);
    %drawnow
end
soundsc(out,fs)

function interp = interpolation1D(u, l, alpha)
    interp= (1-alpha)*u(l)+alpha*u(l+1);
end

function interp = interpolation2D(u, li, mi, alphaX, alphaY)
    interp=(1-alphaX)*(1-alphaY)*u(li,mi)+((1-alphaX)*alphaY)*u(li,mi+1)...
        +(alphaX*(1-alphaY))*u(li+1,mi)+(alphaX*alphaY)*u(li+1,mi+1);
end



%uNext(l+m*Nx)=(2-20*mu^2-4*S)*u(l+m*Nx)...
%                +(8*mu^2+S)*(u(l+1+m*Nx)+u(l-1+m*Nx)+u(l+(m+1)*Nx)+u(l+(m-1)*Nx))...
%                -2*mu^2*(u(l+1+(m+1)*Nx)+u(l-1+(m+1)*Nx)+u(l+1+(m-1)*Nx)+u(l-1+(m-1)*Nx))...
%                -mu^2*(u(l+2+m*Nx)+u(l-2+m*Nx)+u(l+(m+2)*Nx)+u(l+(m-2)*Nx))...
%                +(sigma0*k-1+4*S)*uPrev(l+m*Nx)...
%                -S*(uPrev(l+1+m*Nx)+uPrev(l-1+m*Nx)+uPrev(l+(m+1)*Nx)+uPrev(l+(m-1)*Nx))...
%                +J*malletForce;
