#include <stdlib.h>
#include <TF1.h>
#include <math.h>
#include <TMath.h>
#include "CFGlobal.h"
#include "ReadPar.h"
#include "SourceModelCMSHyperBt.h"
#define XSWITCH -0.53

// #define PAROUT2 0.7
// #define PAROUT4 1.0
// #define PAROUT5 0.0
// #define PAROUT6 2.0
// #define PAROUT7 0.69

// #define PARSIDE2 2.35
// #define PARLONG2 0.67
// #define PARTIME2 2.1

#define PAROUT2 45.0
#define PAROUT4 1.0
#define PAROUT5 0.0
#define PAROUT6 0.4
#define PAROUT7 1.73

#define PARSIDE2 5.0
#define PARLONG2 0.6
#define PARTIME2 2.0

TRandom *SourceModelCMSHyperBt::mRandom = new TRandom();

SourceModelCMSHyperBt::SourceModelCMSHyperBt()
{
  mNPar = 1;
  for (int iter=0; iter<mNPar; iter++) {
    mParameters[iter] = 0.0;
  }
  
  mNRandVar = 4;
  mRandom->SetSeed(23133);
}

SourceModelCMSHyperBt::~SourceModelCMSHyperBt()
{
}

void       
SourceModelCMSHyperBt::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++) 
    mParameters[ti] = aParTable[ti];
  // Parameters for RQMD 0,1
//   mSigmaOut  = 7.3+mParameters[0]*(-1.2);
//   mMeanOut   = -2.32;
//   mSigmaSide = 5.3+mParameters[0]*(-2.2);
//   mSigmaLong = 21.3+mParameters[0]*(-14.5);
//   mSigmaTime = 9.92+mParameters[0]*(-2.55);
//   mMeanTime  = 4.1;
// Parameters for RQMD 1,1
//   mSigmaOut  = 7.67+mParameters[0]*(-1.47);
//   mMeanOut   = -2.05;
//   mSigmaSide = 6.3+mParameters[0]*(-2.63);
//   mSigmaLong = 22.5+mParameters[0]*(-15.8);
//   mSigmaTime = 10.1+mParameters[0]*(2.9)+mParameters[0]*mParameters[0]*(-5.3);
//   mMeanTime  = 5.2;
// Parameters for BFPW Bt6
  mSigmaOut  = 2.72*mParameters[0];
  mMeanOut   = 0.0;
  mSigmaSide = 2.77*mParameters[0];
  mSigmaLong = 7.00*mParameters[0];
  mSigmaTime = 21.77*mParameters[0];
  mMeanTime  = 0.0;
  PRINT_DEBUG_3("Pars:" << mSigmaOut << " " << mMeanOut << "  " << mSigmaSide << " " << mSigmaLong << " " << mSigmaTime << " " << mMeanTime);
}

double*
SourceModelCMSHyperBt::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelCMSHyperBt::InitRandVar(double** aRandVar)
{
  Double_t tTemp;
  
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti=0; ti<mNRandVar; ti++)
    (*aRandVar)[ti] = -100.0;
  (*aRandVar)[0] = GetHyperDouble(0.0, PAROUT2, 1.0, 0.0, PAROUT6*PAROUT2, PAROUT7);
  (*aRandVar)[1] = GetHyper(0.0, PARSIDE2, 1.0);
  (*aRandVar)[2] = GetHyper(0.0, PARLONG2, 1.0);
  (*aRandVar)[3] = GetHyper(0.0, PARTIME2, 1.0);
}

void
SourceModelCMSHyperBt::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
{
  aPos->part2.x=0.;
  aPos->part2.y=0.;
  aPos->part2.z=0.;
  aPos->part2.t=0.;
  
  aRandVar[0] = GetHyperDouble(0.0, PAROUT2, 1.0, 0.0, PAROUT6*PAROUT2, PAROUT7);
  aRandVar[1] = GetHyper(0.0, PARSIDE2, 1.0);
  aRandVar[2] = GetHyper(0.0, PARLONG2, 1.0);
  aRandVar[3] = GetHyper(0.0, PARTIME2, 1.0);
  
  double tPx = aMom->part1.x+aMom->part2.x;
  double tPy = aMom->part1.y+aMom->part2.y;
  double tPz = aMom->part1.z+aMom->part2.z;
  double tE  = aMom->part1.t+aMom->part2.t;
  double tPt = tPx*tPx + tPy*tPy;
  //mCVK = tPz*tPz;
  double tMt = tE*tE - tPz*tPz;//mCVK;
  //mCVK += tPt;
  //mCVK = sqrt(mCVK);
  double tM =   sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  
//   double tROut  = GetHyperDouble(mParameters[1], PAROUT2, mParameters[0], mParameters[1], PAROUT6*PAROUT2, mParameters[0]*PAROUT7);
//   double tRSide = GetHyper(0.0, PARSIDE2, mParameters[2]);
//   aPos->part1.z = GetHyper(0.0, PARLONG2, mParameters[3]);
//   aPos->part1.t = GetHyper(0.0, PARTIME2, mParameters[4]);

//   double tROut  = TMath::Power(aRandVar[0],mParInverse[0])+mParameters[1];
//   double tRSide = TMath::Power(aRandVar[1],mParInverse[2]);
//   aPos->part1.z = TMath::Power(aRandVar[2],mParInverse[3]);
//   aPos->part1.t = TMath::Power(aRandVar[3],mParInverse[4]);

  double tROut  = ScaleHyper(mSigmaOut, aRandVar[0]) + mMeanOut;
  double tRSide = ScaleHyper(mSigmaSide, aRandVar[1]);
  aPos->part1.z = ScaleHyper(mSigmaLong, aRandVar[2]);
  aPos->part1.t = ScaleHyper(mSigmaTime, aRandVar[3]) + mMeanTime;
  
  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part1.x = tROut*tPx-tRSide*tPy;
  aPos->part1.y = tROut*tPy+tRSide*tPx;
}

const char*       
SourceModelCMSHyperBt::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "Bt range";
    break;
  }
}

const char* 
SourceModelCMSHyperBt::GetModelIdentifier()
{
  return "SMExpHyperbolicCMSBt";
}

double 
SourceModelCMSHyperBt::GetHyper(Double_t aMu, Double_t aAlfa, Double_t aSigma)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;
  Double_t tSigSq = aSigma*aSigma;
  Double_t tAlSq  = aAlfa*aAlfa;
  Double_t tAlExp = TMath::Exp(aAlfa/aSigma);
  Double_t tHalfWidth = TMath::Sqrt(TMath::Power((aAlfa-TMath::Log(0.5)),2)-tAlSq);
  
  do {
    tTemp = (mRandom->Rndm()*16.0-8.0)*tHalfWidth;
    tArgSq = (tTemp-aMu)*(tTemp-aMu);
    tVal  = tAlExp*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq)+tAlSq));
  }
  while (tVal < mRandom->Rndm());

  return tTemp;
}

double 
SourceModelCMSHyperBt::GetHyperDouble(Double_t aMu1, Double_t aAlfa1, Double_t aSigma1, Double_t aMu2, Double_t aAlfa2, Double_t aSigma2)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;

  Double_t tSigSq1 = aSigma1*aSigma1;
  Double_t tAlSq1  = aAlfa1*aAlfa1;
  Double_t tAlExp1 = TMath::Exp(aAlfa1);

  Double_t tSigSq2 = aSigma2*aSigma2;
  Double_t tAlSq2  = aAlfa2*aAlfa2;
  Double_t tAlExp2 = TMath::Exp(aAlfa2);

  Double_t tHalfWidth = (aAlfa1>aAlfa2) ? TMath::Sqrt(TMath::Power((aAlfa1-TMath::Log(0.5)),2)-tAlSq1) : TMath::Sqrt(TMath::Power((aAlfa2-TMath::Log(0.5)),2)-tAlSq2);

  do {
    tTemp = (mRandom->Rndm()*16.0-8.0)*tHalfWidth;
    if (tTemp > XSWITCH) {
      tArgSq = (tTemp-aMu2)*(tTemp-aMu2);
      tVal  = tAlExp2*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq2)+tAlSq2));
    }
    else {
      tArgSq = (tTemp-aMu1)*(tTemp-aMu1);
      tVal  = tAlExp1*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq1)+tAlSq1));
    }
  }
  while (tVal < mRandom->Rndm());

  return tTemp;
}

inline double 
SourceModelCMSHyperBt::ScaleHyper(Double_t aScale, Double_t aX)
{
  return aX*aScale;
}


