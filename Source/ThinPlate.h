/*
  ==============================================================================

    ThinPlate.h
    Created: 21 Sep 2022 3:33:46pm
    Author:  Benjamin Støier

  ==============================================================================
*/


#pragma once

#include <JuceHeader.h>



class ThinPlate  : public juce::Component
{
public:
    ThinPlate (double k); // initialise the model with the time step
    ~ThinPlate() override;
    
void initParameters();

void updateParameters(const double sig0ToSet, const double sig1ToSet, const double LxToSet, const double LyToSet, const double excXToSet, const double excYToSet, const double lisXToSet, const double lisYToSet, const double thicknessToSet, const double excFToSet, const double excTToSet, const double vBToSet, const double fBToSet, const double aToSet, const int excTypeId, const double  bAtt1ToSet, const double bDec1ToSet, const double  bSus1ToSet, const double bRel1ToSet, const double FBEnv1ToSet, const double vBEnv1ToSet, const double lfoRateToSet, const double xPosModToSet, const double yPosModToSet);
    
void updatePlateMaterial(int plateMaterialToSet);
    
void getSampleRate(double fsToSet);

void calculateScheme();
    
void plateHit();
    
void startBow();
    
void endBow();

void updateStates();

float getOutput()
{
    switch (excType) {
        case Mallet:
            return u[static_cast <int> (floor(lisXpos*Nx))][static_cast <int> (floor(lisYpos*Ny))]*0.000001;
            break;
        case Bow:
            return u[static_cast <int> (floor(lisXpos*Nx))][static_cast <int> (floor(lisYpos*Ny))]*0.00001;
            break;
    }
    
};
    
float interpolation(std::vector<double>* u, int excXidx, int excYidx, double alphaX, double alphaY)
{
    if (excXidx < 0)
    {
        excXidx = 0;
    }
    if (excYidx < 0)
    {
        excYidx = 0;
    }
    
    if (excXidx > Nx)
    {
        excXidx = Nx;
    }
    if (excYidx > Ny)
    {
        excYidx = Ny;
    }
    return (1.0 - alphaX) * (1.0 - alphaY) * u[excXidx][excYidx]+(1.0-alphaX)*alphaY*u[excXidx][excYidx+1]+alphaX*(1-alphaY)*u[excXidx+1][excYidx]+alphaX*alphaY*u[excXidx+1][excYidx+1];
}
  
void setADSR(double sampleRate);
    
private:
    enum ExcitationType
    {
        Mallet,
        Bow
    };
    
    juce::ADSR adsr1;
    juce::ADSR::Parameters adsr1Params;

    // Variables
    double k; // Time step (in s)
    double c; // Wave speed (in m/s)
    double h; // Grid spacing (in m)
    double Lx; // Length of x dimension(in m)
    double Ly; // Length of y dimension(in m)
    
    double T60;
    double sigma0; // frequency independent damping
    double sigma1; //frequency dependent damping
    double rho; //material density
    double H; //plate thickness
    double E; //Young's modulus
    double nu; //Poisson's ratio
    
    double D; // stifness coefficient
    double kappa; //stifness parameter
    double mu;
    double muSq;
    double S;
    double maxForce;
    double excTime;
    double excXpos;
    double excYpos;
    double excXposRatio;
    double excYposRatio;
    double lisXpos;
    double lisYpos;
    double hx;
    double hy;
    double lambda;
    int excXidx;
    int excYidx;
    double b;
    double a;
    double vB, FB, vRel, vRelPrev, eps, tol;
    
    double alphaX;
    double alphaY;
    double J;
    double malletForce;
    
    double angleDeltaLFO;
    double currentXMod;
    double currentYMod;
    double currentAngleLFO;
    double lfoRate;
    
    double xPosMod;
    double yPosMod;
    
    int n;
    int plateMaterial;
    
    int numStrings;
    
    double K1, K3, R; // Connection parameters
    std::vector<double> connSPos, connXPos, connYPos, connSPos2, connXPos2, connYPos2; //Connection posisistion (contenios domain)
    std::vector<int> lcS, lcP, mcP, lcS2, lcP2, mcP2, lcT; //Connection posistions (discrete domain)
    std::vector<double> alphaConnS, alphaConnX, alphaConnY, alphaConnS2, alphaConnX2, alphaConnY2;
    std::vector<double> connF; // connection force
    std::vector<double> eta, etaPrev, etaNext; //Connection distance
    std::vector<double> rPlus, rMinus;
    std::vector<double> stringConnTerm;
    double plateConnTerm, tubeConnTerm;
    
    bool stringConn, tubeConn;
    
    int Nx; // grid steps (x)
    int Ny; // grid steps (y)
    int N; // grid steps (total)
    
    ExcitationType excType;
    
    double fs; //sample rate
    double t;
    double t0;
    
    double FBEnv1, vBEnv1;
    double bAtt1, bDec1, bSus1, bRel1;
    double nextAdsr1;
    bool isBowing, bowEnd;
    
    double excitation;
    std::vector<double>* uNext;
    std::vector<double>* u; // state pointers
    std::vector<double>* uPrev;
    std::vector<std::vector<std::vector<double>>> uStates;
    
    //String parameters
    double LS; //Length
    double rhoS; //Density
    double rS; //Radius
    double AS; // Cross-sectional area
    double ES; // Young's modulus
    double IS; //Area moment of inertia
    double sigma0S, sigma1S; // frequency-independent damping and frequency-dependent damping
    double kappaSSq; // Stiffness term
    double As;
    double TS;
    double cSSq;
    double NS;
    std::vector<double> hS;
    double lambdaSSq;
    double muSSq;
    double uS1, uS2, uS3;
    
    std::vector<double*> uString;
    std::vector<std::vector<double>> uStringStates;
    
    
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ThinPlate)
};


