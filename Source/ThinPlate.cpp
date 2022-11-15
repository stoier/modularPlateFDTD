/*
  ==============================================================================

    ThinPlate.cpp
    Created: 21 Sep 2022 3:33:46pm
    Author:  Benjamin Støier

  ==============================================================================
*/

#include "ThinPlate.h"


//==============================================================================
ThinPlate::ThinPlate (double kIn) : k (kIn) // <- This is an initialiser list. It initialises the member variable 'k' (in the "private" section in OneDWave.h), using the argument of the constructor 'kIn'.
{
    Lx= 0.5;
    Ly = 0.5; //side length (y)
    c = 343; //speed of sound
    sigma0 = 1; //frequency independent damping
    sigma1 = 0.005; //frequency dependent damping
    rho = 7700; //material density
    H = 5 * pow(10,-3); //Thickness
    E = 2 * pow(10,11); //Young's modulus
    nu = 0.3; //Poisson's ratio
    
    D = E*pow(H,3)/(12*(1-pow(nu,2))); // stifness coefficient
    kappa = sqrt(D/(rho*H)); // stifness paramater
    
    h = 2*sqrt(k*(sigma1+sqrt(pow(kappa,2)+pow(sigma1,2))));
    Nx = floor(Lx/h); //grid steps (x)
    Ny = floor(Ly/h); //grid steps (y)
    h = std::min (Lx / Nx, Ly / Ny); //redefine grid spacing (using the smallest dimension)

    N = (Nx+1)*(Ny+1);
    mu = kappa*k/pow(h,2);
    muSq = mu * mu;
    S = 2 * sigma1 * k /pow(h,2);
    maxForce = 10; //maximum force
    excTime = 0.001;
    t0 = 0;
    t = 0;
    excXpos = 0.3;
    excYpos = 0.3;
    hx = Lx/Nx;
    hy = Ly/Ny;
    excXidx = floor(excXpos/hx);
    excYidx = floor(excYpos/hy);
    alphaX = excXpos/hx-excXidx;
    alphaY = excYpos/hy-excYidx;
    isBowing = false;
    bowEnd = true;
    lisXpos = 0.5;
    lisYpos = 0.5;
    n=0;
    
    J=0;
    
    
    uStates = std::vector<std::vector<std::vector<double>>> (3, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0)));
    uNext = &uStates[0][0];
    u = &uStates[1][0];
    uPrev = &uStates[2][0];
    
    numStrings = 1;
    uStringStates = std::vector<std::vector<double>> (3,
                                           std::vector<double>(NS+1, 0));
    uString.resize (3, nullptr);
    
    for (int i = 0; i < 3; ++i)
            uString[i] = &uStringStates[i][0];
    
    eta.resize(numStrings);
    etaPrev.resize(numStrings);
    etaNext.resize(numStrings);
    rPlus.resize(numStrings);
    rMinus.resize(numStrings);
    connF.resize(numStrings);
    connSPos.resize(numStrings);
    connXPos.resize(numStrings);
    connYPos.resize(numStrings);
    connSPos2.resize(numStrings);
    connXPos2.resize(numStrings);
    connYPos2.resize(numStrings);
    lcS.resize(numStrings);
    lcP.resize(numStrings);
    mcP.resize(numStrings);
    lcS2.resize(numStrings);
    lcP2.resize(numStrings);
    mcP2.resize(numStrings);
    alphaConnS.resize(numStrings);
    alphaConnX.resize(numStrings);
    alphaConnY.resize(numStrings);
    alphaConnS2.resize(numStrings);
    alphaConnX2.resize(numStrings);
    alphaConnY2.resize(numStrings);
    stringConnTerm.resize(numStrings);
}

ThinPlate::~ThinPlate()
{
}

void ThinPlate::getSampleRate(double fsToSet)
{
    fs = fsToSet;
}



void ThinPlate::initParameters()
{
    D = E*pow(H,3)/(12*(1-pow(nu,2))); // stifness coefficient
    kappa = sqrt(D/(rho*H)); // stifness paramater
    h = 2*sqrt(k*(sigma1+sqrt(pow(kappa,2)+pow(sigma1,2))));
    Nx = floor(Lx/h); //grid steps (x)
    Ny = floor(Ly/h); //grid steps (y)
    h = std::min (Lx / Nx, Ly / Ny); //redefine grid spacing (using the smallest dimension)
    N = (Nx+1)*(Ny+1);
    
    mu = kappa*k/pow(h,2);
    muSq = mu * mu;
    S = 2 * sigma1 * k /pow(h,2);
    
    hx = Lx/Nx;
    hy = Ly/Ny;
    
    excXpos = excXposRatio*Lx;
    excYpos = excYposRatio*Ly;
    excXidx = floor(excXpos/hx);
    excYidx = floor(excYpos/hy);
    alphaX = excXpos/hx-excXidx;
    alphaY = excYpos/hy-excYidx;
    
    J=0;
    
    uStates = std::vector<std::vector<std::vector<double>>> (3, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0)));
    uNext = &uStates[0][0]; //Initialise time step u^n+1
    u = &uStates[1][0]; //Initialise time step u^n
    uPrev = &uStates[2][0]; //Initialise time step u^n-1
    
    uStringStates = std::vector<std::vector<double>> (3,
                                           std::vector<double>(NS+1, 0));
    uString.resize (3, nullptr);
    
    for (int i = 0; i < 3; ++i)
            uString[i] = &uStringStates[i][0];
    
    eta.resize(numStrings);
    etaPrev.resize(numStrings);
    etaNext.resize(numStrings);
    rPlus.resize(numStrings);
    rMinus.resize(numStrings);
    connF.resize(numStrings);
    connSPos.resize(numStrings);
    connXPos.resize(numStrings);
    connYPos.resize(numStrings);
    connSPos2.resize(numStrings);
    connXPos2.resize(numStrings);
    connYPos2.resize(numStrings);
    lcS.resize(numStrings);
    lcP.resize(numStrings);
    mcP.resize(numStrings);
    lcS2.resize(numStrings);
    lcP2.resize(numStrings);
    mcP2.resize(numStrings);
    alphaConnS.resize(numStrings);
    alphaConnX.resize(numStrings);
    alphaConnY.resize(numStrings);
    alphaConnS2.resize(numStrings);
    alphaConnX2.resize(numStrings);
    alphaConnY2.resize(numStrings);
    stringConnTerm.resize(numStrings);
    hS.resize(numStrings);
    
    for (int nS = 0; nS < numStrings; ++nS)
    {
        lcS[nS]  = floor(connSPos[nS]/hS[nS]);
        //alphaConnS[nS] = connSPos[nS]/hS[nS] - lcS(i)
        lcP[nS] = floor(connXPos[nS]/hx);
        mcP[nS] = floor(connYPos[nS]/hy);
        lcS2[nS]  = floor(connSPos2[nS]/hS[nS]);
        lcP2[nS] = floor(connXPos2[nS]/hx);
        mcP2[nS] = floor(connYPos2[nS]/hy);
    }
    
    plateConnTerm = (k*k)/(rho*H*h*h*(1+sigma0*k));
}

void ThinPlate::updateParameters(const double sig0ToSet, const double sig1ToSet, const double LxToSet, const double LyToSet, const double excXToSet, const double excYToSet, const double lisXToSet, const double lisYToSet, const double thicknessToSet, const double excFToSet, const double excTToSet, const double vBToSet, const double FBToSet, const double aToSet, const int excTypeId, const double  bAtt1ToSet, const double bDec1ToSet, const double  bSus1ToSet, const double bRel1ToSet, const double FBEnv1ToSet, const double vBEnv1ToSet, const double lfoRateToSet, const double xPosModToSet, const double yPosModToSet)
{
    sigma0 = sig0ToSet;
    sigma1 = sig1ToSet;
    Lx = LxToSet;
    Ly = LyToSet;
    excXposRatio = excXToSet;
    excYposRatio = excYToSet;
    lisXpos = lisXToSet;
    lisYpos = lisYToSet;
    lfoRate= lfoRateToSet;
    xPosMod = xPosModToSet/100;
    yPosMod = yPosModToSet/100;
    
    H = thicknessToSet*0.001f;
    
    if (excTypeId == 1)
    {
        vB = vBToSet;
        FB = FBToSet;
        /*
        D = E*pow(H,3)/(12*(1-pow(nu,2))); // stifness coefficient
        kappa = sqrt(D/(rho*H)); // stifness paramater
        h = 2*sqrt(k*(sigma1+sqrt(pow(kappa,2)+pow(sigma1,2))));
        Nx = floor(Lx/h); //grid steps (x)
        Ny = floor(Ly/h); //grid steps (y)
        h = std::min (Lx / Nx, Ly / Ny); //redefine grid spacing (using the smallest dimension)
        N = (Nx+1)*(Ny+1);
        mu = kappa*k/pow(h,2);
        muSq = mu * mu;
        S = 2 * sigma1 * k /pow(h,2);
        uStates = std::vector<std::vector<std::vector<double>>> (3, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0)));
        uNext = &uStates[0][0]; //Initialise time step u^n+1
        u = &uStates[1][0]; //Initialise time step u^n
        uPrev = &uStates[2][0]; //Initialise time step u^n-1
        hx = Lx/Nx;
        hy = Ly/Ny;
         */
        
        a = aToSet * 0.00000000001;
        //a =   0.00000000001;
        FBEnv1 = FBEnv1ToSet;
        vBEnv1 = vBEnv1ToSet;
        bAtt1 = bAtt1ToSet;
        bDec1 = bDec1ToSet;
        bSus1 = bSus1ToSet;
        bRel1 = bRel1ToSet;
        setADSR(fs);
        excType = Bow;
    }
    else
    {
        maxForce = excFToSet;
        excTime = excTToSet*0.001f;
        excType = Mallet;
    }
}

void ThinPlate::updatePlateMaterial(int plateMaterialToSet)
{
    plateMaterial=plateMaterialToSet;
    if (plateMaterial == 1) //brass
    {
        rho=8530;
        E = 110 * pow(10, 9);
        nu = 0.34;
    }
    else if (plateMaterial == 2) //bronze
    {
        rho = 8770;
        E = 103 * pow(10, 9);
        nu = 0.34;
    }
    else if (plateMaterial == 3) //iron
    {
        rho = 7874;
        E = 211 * pow(10, 9);
        nu = 0.26;
    }
    else if (plateMaterial == 4) // aluminium
    {
        rho = 2700;
        E = 70 * pow(10, 9);
        nu = 0.33;
    }
    else if (plateMaterial == 5) //gold
    {
        rho = 19300;
        E = 79 * pow(10, 9);
        nu = 0.42;
    }
    else if (plateMaterial == 6) //silver
    {
        rho = 10490;
        E = 83 * pow(10, 9);
        nu = 0.37;
    }
    else if (plateMaterial == 7) //copper
    {
        rho = 8920;
        E = 120 * pow(10, 9);
        nu = 0.36;
    }
}

void ThinPlate::plateHit()
{
    n = 0;
    t = 0;
}

void ThinPlate::startBow()
{
    if (bowEnd == true)
    {
        adsr1.noteOn();
        isBowing = true;
        bowEnd = false;
    }
}

void ThinPlate::endBow()
{
    adsr1.noteOff();
    isBowing = false;
    bowEnd = true;
}

void ThinPlate::calculateScheme()
{
    //Get number of cycles pr. sample the number of cycles pr. sample for the lfo
    auto cyclesPerSampleLFO= lfoRate / fs;
    
    //Update the amount that the phase angle of the vibrato needs to increment for each sample
    angleDeltaLFO = cyclesPerSampleLFO * 2.0 * juce::MathConstants<double>::pi;
    
    currentXMod = std::sin(currentAngleLFO)*xPosMod;
    currentYMod = std::sin(currentAngleLFO)*yPosMod;
    currentAngleLFO += angleDeltaLFO;
    
    switch (excType) {
        case Mallet:
            if (n < floor(excTime*fs))
            {
                malletForce= maxForce/2*(1-std::cos((2*juce::MathConstants<double>::pi*(t-t0))/(excTime)));
                t = t+k;
            }
            else
            {
                malletForce = 0;
            }
            excitation = malletForce;
            n++;
            break;
            
        case Bow:
           if (isBowing == true)
           {
               excXpos = excXposRatio*Lx;
               excYpos = excYposRatio*Ly;
               if (0 < xPosMod)
               {
                   excXpos = excXpos + (currentXMod *  excXpos);
               }
               if (0 < yPosMod)
               {
                   excYpos = excYpos + (currentYMod *  excYpos);
               }
               excXidx = floor(excXpos/hx);
               excYidx = floor(excYpos/hy);
               alphaX = excXpos/hx-excXidx;
               alphaY = excYpos/hy-excYidx;
                
               nextAdsr1 = adsr1.getNextSample();
               b = ((2/k) + 2*sigma0)*(vB*nextAdsr1) - ((2/(k*k)*interpolation(u,excXidx,excYidx,alphaX, alphaY)-interpolation(uPrev,excXidx,excYidx,alphaX, alphaY)) + ((kappa*kappa)/(h*h*h*h))*(interpolation(u,excXidx+2,excYidx,alphaX, alphaY)+interpolation(u,excXidx-2,excYidx,alphaX, alphaY)+interpolation(u,excXidx,excYidx+2,alphaX, alphaY)+interpolation(u,excXidx,excYidx-2,alphaX, alphaY))
                + 2 * (interpolation(u,excXidx+1,excYidx+1,alphaX, alphaY)+interpolation(u,excXidx+1,excYidx-1,alphaX, alphaY)+interpolation(u,excXidx-1,excYidx+1,alphaX, alphaY)+interpolation(u,excXidx-1,excYidx-1,alphaX, alphaY)-8*(interpolation(u,excXidx+1,excYidx,alphaX, alphaY)+interpolation(u,excXidx-1,excYidx,alphaX, alphaY)+interpolation(u,excXidx,excYidx+1,alphaX, alphaY)+interpolation(u,excXidx,excYidx-1,alphaX, alphaY))+20*interpolation(u,excXidx,excYidx,alphaX, alphaY))
                - 2*sigma1/(k*h*h)*(interpolation(u,excXidx+1,excYidx,alphaX, alphaY)+interpolation(u,excXidx-1,excYidx,alphaX, alphaY)+interpolation(u,excXidx,excYidx+1,alphaX, alphaY)+interpolation(u,excXidx,excYidx-1,alphaX, alphaY)-interpolation(u,excXidx+1,excYidx,alphaX, alphaY)-interpolation(u,excXidx-1,excYidx,alphaX, alphaY)-interpolation(u,excXidx,excYidx+1,alphaX, alphaY)-interpolation(u,excXidx,excYidx-1,alphaX, alphaY)-4*(interpolation(u,excXidx,excYidx,alphaX, alphaY)-interpolation(uPrev,excXidx,excYidx,alphaX, alphaY))));
                eps = 1;
                int i = 0;
                while (eps > tol && i < 200)
                {
                    vRel =  vRelPrev - (((2/k+2*sigma0)*vRelPrev+(FB*nextAdsr1)*sqrt(2*a)*vRelPrev*exp(-a*vRelPrev*vRelPrev+0.5)+b)/(2/k + 2*sigma0+(FB*nextAdsr1)*sqrt(2*a)*(1-2*a*vRel*vRel)*exp(-a*(vRel*vRel+0.5))));
                    eps = std::abs(vRel-vRelPrev);
                    vRelPrev = vRel;
                    ++i;
                }
                excitation = sqrt(2*a)*vRel*exp(-a*vRel*vRel+0.5)*(FB*nextAdsr1);
           }
            break;
    }
    
    for (int l = 2; l < Nx-2; ++l) // clamped boundaries
    {
        for (int m = 2; m < Ny-2; ++m) // clamped boundaries
        {
            if (l == excXidx && m == excYidx)
                J = ((1-alphaX)*(1-alphaY))/(hx*hy);
            else if (l == excXidx && m == excYidx+1)
                J = ((1-alphaX)*alphaY)/(hx*hy);
            else if (l == excXidx+1 && m == excYidx)
                J = (alphaX*(1-alphaY))/(hx*hy);
            else if (l == excXidx+1 && m == excYidx+1)
                J = (alphaX*alphaY)/(hx*hy);
            else
                J = 0;
            
            uNext[l][m] =
            (2-20*muSq-4*S)*u[l][m]
            + (8*muSq+S) * (u[l+1][m] + u[l-1][m] + u[l][m+1] + u[l][m-1])
            - 2*muSq * (u[l+1][m+1] + u[l-1][m+1] + u[l+1][m-1] + u[l-1][m-1])
            - muSq * (u[l+2][m] + u[l-2][m] + u[l][m+2] + u[l][m-2])
            + (sigma0*k-1+4*S) * uPrev[l][m]
            - S * (uPrev[l+1][m] + uPrev[l-1][m] + uPrev[l][m+1] + u[l][m-1])
            + J * excitation;
        }
    }
    if (stringConn == true)
    {
        for (int nS = 0; nS < numStrings; ++nS)
        {
            for (int l = 2; l < NS-1; l++)
            {
                uString[0][l] =
                uS1 * uString[1][l]
                + uS2 * (uString[1][l+1]+uString[1][l-1]) - muSSq * (uString[1][l+2]+uString[1][l-2])
                + uS3 * uString[2][l]-2*sigma1S*k/(hS[nS]*hS[nS])*(uString[2][l+1]+uString[2][l-1])/As;
            }
            etaNext[nS] = uString[0][lcS[nS]]-uNext[lcP[nS]][mcP[nS]];
            eta[nS] = uString[1][lcS[nS]]-u[lcP[nS]][mcP[nS]];
            etaPrev[nS] = uString[2][lcS[nS]]-uPrev[lcP[nS]][mcP[nS]];
            rPlus[nS] = 0.5 * K1 + 0.5 * K3 * eta[nS] * eta[nS] + 0.5 * fs * R;
            rMinus[nS] = 0.5 * K1 + 0.5 * K3 * eta[nS] * eta[nS] + 0.5 * fs * R;
            connF[nS]  = (etaNext[nS] + rMinus[nS] / rPlus[nS] * etaPrev[nS]) / (1/rPlus[nS] + plateConnTerm + stringConnTerm[nS]);
        }
    }
    updateStates();
}

void ThinPlate::updateStates()
{
        std::vector<double>* uTmp = uPrev;
        uPrev = u;
        u= uNext;
        uNext = uTmp;
    
        double* uStringTmp = uString[2];
        uString[2] = uString[1];
        uString[1] = uString[0];
        uString[0] = uStringTmp;
}

void ThinPlate::setADSR(double sampleRate)
{
    adsr1.setSampleRate(sampleRate);
    adsr1Params.attack = bAtt1;
    adsr1Params.decay = bDec1;
    adsr1Params.sustain = bSus1;
    adsr1Params.release = bRel1;
    adsr1.setParameters(adsr1Params);
    
}

