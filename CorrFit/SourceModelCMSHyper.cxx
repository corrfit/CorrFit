#include <stdlib.h>
#include <TF1.h>
#include "CFGlobal.h"
#include "ReadPar.h"
#include "SourceModelCMSHyper.h"
#include <math.h>
#include <TMath.h>
#define XSWITCH -0.0

#define PAROUT2 66.2
#define PAROUT4 1.0
#define PAROUT5 0.0
#define PAROUT6 1.0
#define PAROUT7 1.0

#define PARSIDE2 67.0
#define PARLONG2 0.48
#define PARTIME2 1.8

TRandom *SourceModelCMSHyper::mRandom = new TRandom();

SourceModelCMSHyper::SourceModelCMSHyper()
{
  mNPar = 5;
  for (int iter=0; iter<mNPar; iter++) {
    mParameters[iter] = 0.0;
  }
  
  mNRandVar = 4;
  mRandom->SetSeed(23133);
}

SourceModelCMSHyper::~SourceModelCMSHyper()
{
}

void       
SourceModelCMSHyper::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++) 
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelCMSHyper::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelCMSHyper::InitRandVar(double** aRandVar)
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
SourceModelCMSHyper::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
{
  aPos->part2.x=0.;
  aPos->part2.y=0.;
  aPos->part2.z=0.;
  aPos->part2.t=0.;
  
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

  double tROut  = ScaleHyper(mParameters[0], aRandVar[0]) + mParameters[1];
  double tRSide = ScaleHyper(mParameters[2], aRandVar[1]);
  aPos->part1.z = ScaleHyper(mParameters[3], aRandVar[2]);
  aPos->part1.t = ScaleHyper(mParameters[4], aRandVar[3]);
  
  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part1.x = tROut*tPx-tRSide*tPy;
  aPos->part1.y = tROut*tPy+tRSide*tPx;
}

const char*       
SourceModelCMSHyper::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "ROut distribution sigma";
    break;
  case 1:
    return "ROut distribution mean";
    break;
  case 2:
    return "RSide distribution sigma";
    break;
  case 3:
    return "RLong distribution sigma";
    break;
  case 4:
    return "RTime distribution sigma";
    break;
  }
}

const char* 
SourceModelCMSHyper::GetModelIdentifier()
{
  return "SMExpHyperbolicCMS";
}

double 
SourceModelCMSHyper::GetHyper(Double_t aMu, Double_t aAlfa, Double_t aSigma)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;
  Double_t tSigSq = aSigma*aSigma;
  Double_t tAlSq  = aAlfa*aAlfa;
  Double_t tAlExp = TMath::Exp(aAlfa/aSigma);

  do {
    tTemp = (mRandom->Rndm()*16.0-8.0)*aSigma;
    tArgSq = (tTemp-aMu)*(tTemp-aMu);
    tVal  = tAlExp*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq)+tAlSq));
  }
  while (tVal < mRandom->Rndm());

  return tTemp;
}

double 
SourceModelCMSHyper::GetHyperDouble(Double_t aMu1, Double_t aAlfa1, Double_t aSigma1, Double_t aMu2, Double_t aAlfa2, Double_t aSigma2)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;

  Double_t tSigSq1 = aSigma1*aSigma1;
  Double_t tAlSq1  = aAlfa1*aAlfa1;
  Double_t tAlExp1 = TMath::Exp(aAlfa1);

  Double_t tSigSq2 = aSigma2*aSigma2;
  Double_t tAlSq2  = aAlfa2*aAlfa2;
  Double_t tAlExp2 = TMath::Exp(aAlfa2);

  Double_t tSigma = (aSigma1>aSigma2) ? aSigma1 : aSigma2;
  do {
    tTemp = (mRandom->Rndm()*16.0-8.0)*tSigma;
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
SourceModelCMSHyper::ScaleHyper(Double_t aScale, Double_t aX)
{
  return aX*aScale;
}


