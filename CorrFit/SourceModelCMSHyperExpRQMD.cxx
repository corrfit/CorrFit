#include <stdlib.h>
#include <TF1.h>
#include "CFGlobal.h"
#include "ReadPar.h"
#include "SourceModelCMSHyperExpRQMD.h"
#define XSWITCH -0.13
#define TSWITCH 0.13

#define PAROUT2 2.0
#define PAROUT4 1.0
#define PAROUT5 0.0
#define PAROUT6 1.2
#define PAROUT7 0.8
#define PAROUT9 0.031
#define PAROUT10 0.025*5.09

#define PARTIME2 2.0
#define PARTIME4 1.0
#define PARTIME5 0.0
#define PARTIME6 1.2
#define PARTIME7 0.8

#define PARSIDE2 2.0
#define PARSIDE4 0.031
#define PARSIDE5 0.025*4.59

#define PARLONG2 2.0
#define PARLONG4 0.21
#define PARLONG5 0.025*10.6

TRandom *SourceModelCMSHyperExpRQMD::mRandom = new TRandom();

SourceModelCMSHyperExpRQMD::SourceModelCMSHyperExpRQMD()
{
  Double_t tLogHalf = TMath::Log(0.5);

  mNPar = 5;
  for (int iter=0; iter<mNPar; iter++) {
    mParameters[iter] = 0.0;
  }
  
  mNRandVar = 4;
  mRandom->SetSeed(23133);
  mHalfWides[0] = 1.0*TMath::Sqrt(TMath::Power(PAROUT2-tLogHalf,2)-PAROUT2*PAROUT2);
  mHalfWides[1] = 1.0*TMath::Sqrt(TMath::Power(PARSIDE2-tLogHalf,2)-PARSIDE2*PARSIDE2);
  mHalfWides[2] = 1.0*TMath::Sqrt(TMath::Power(PARLONG2-tLogHalf,2)-PARLONG2*PARLONG2);
  mHalfWides[3] = 1.0*TMath::Sqrt(TMath::Power(PARTIME2-tLogHalf,2)-PARTIME2*PARTIME2);
}

SourceModelCMSHyperExpRQMD::~SourceModelCMSHyperExpRQMD()
{
}

void       
SourceModelCMSHyperExpRQMD::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++) 
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelCMSHyperExpRQMD::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelCMSHyperExpRQMD::InitRandVar(double** aRandVar)
{
  Double_t tTemp;
  
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti=0; ti<mNRandVar; ti++)
    (*aRandVar)[ti] = -100.0;
  (*aRandVar)[0] = GetHyperDoubleExp(0.0, PAROUT2, 1.0, 0.0, PAROUT6*PAROUT2, PAROUT7, mHalfWides[0], XSWITCH, PAROUT9, PAROUT10);
  (*aRandVar)[1] = GetHyperExp(0.0, PARSIDE2, 1.0, mHalfWides[1], PARSIDE4, PARSIDE5);
  (*aRandVar)[2] = GetHyperExp(0.0, PARLONG2, 1.0, mHalfWides[2], PARLONG4, PARLONG5);
  (*aRandVar)[3] = GetHyperDouble(0.0, PARTIME6*PARTIME2, PARTIME7, 0.0, PARTIME2, 1.0, mHalfWides[3], TSWITCH);
  //(*aRandVar)[3] = GetHyper(0.0, PARTIME2, 1.0, mHalfWides[3]);

}

void
SourceModelCMSHyperExpRQMD::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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
  double tRSide = ScaleHyper(mParameters[0]*mParameters[2], aRandVar[1]);
  aPos->part1.z = ScaleHyper(mParameters[0]*mParameters[3], aRandVar[2]);
  aPos->part1.t = ScaleHyper(mParameters[0]*mParameters[4], aRandVar[3]) - mParameters[1];
  PRINT_DEBUG_3("Scaled " << aRandVar[3] << " to " << aPos->part1.t);
  
  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part1.x = tROut*tPx-tRSide*tPy;
  aPos->part1.y = tROut*tPy+tRSide*tPx;
}

const char*       
SourceModelCMSHyperExpRQMD::GetParameterName(int aPar)
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
SourceModelCMSHyperExpRQMD::GetModelIdentifier()
{
  return "SMExpHyperbolicExpCMSRQMD";
}

double 
SourceModelCMSHyperExpRQMD::GetHyper(Double_t aMu, Double_t aAlfa, Double_t aSigma, Double_t aHalfWide)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;
  Double_t tSigSq = aSigma*aSigma;
  Double_t tAlSq  = aAlfa*aAlfa;
  Double_t tAlExp = TMath::Exp(aAlfa/aSigma);

  do {
    tTemp = (mRandom->Rndm()*50.0-25.0)*aHalfWide;
    tArgSq = (tTemp-aMu)*(tTemp-aMu);
    tVal  = tAlExp*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq)+tAlSq));
  }
  while (tVal < mRandom->Rndm());

  return tTemp;
}

double 
SourceModelCMSHyperExpRQMD::GetHyperDouble(Double_t aMu1, Double_t aAlfa1, Double_t aSigma1, Double_t aMu2, Double_t aAlfa2, Double_t aSigma2, Double_t aHalfWide, Double_t aSwitch)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;

  Double_t tSigSq1 = aSigma1*aSigma1;
  Double_t tAlSq1  = aAlfa1*aAlfa1;
  Double_t tAlExp1 = TMath::Exp(aAlfa1);

  Double_t tSigSq2 = aSigma2*aSigma2;
  Double_t tAlSq2  = aAlfa2*aAlfa2;
  Double_t tAlExp2 = TMath::Exp(aAlfa2);

  do {
    tTemp = (mRandom->Rndm()*50.0-25.0)*aHalfWide;
    if (tTemp > aSwitch) {
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

double 
SourceModelCMSHyperExpRQMD::GetHyperExp(Double_t aMu, Double_t aAlfa, Double_t aSigma, Double_t aHalfWide, Double_t aExpScale, Double_t aExpSlope)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;
  Double_t tSigSq = aSigma*aSigma;
  Double_t tAlSq  = aAlfa*aAlfa;
  Double_t tAlExp = TMath::Exp(aAlfa/aSigma);

  do {
    tTemp = (mRandom->Rndm()*50.0-25.0)*aHalfWide;
    tArgSq = (tTemp-aMu)*(tTemp-aMu);
    tVal  = tAlExp*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq)+tAlSq));
    if ((tTemp-aMu) > 0.0)
      tVal += aExpScale*TMath::Exp(-aExpSlope*(tTemp-aMu));
    else
      tVal += aExpScale*TMath::Exp(+aExpSlope*(tTemp-aMu));
  }
  while (tVal < mRandom->Rndm()*(1.0+aExpScale));

  return tTemp;
}

double 
SourceModelCMSHyperExpRQMD::GetHyperDoubleExp(Double_t aMu1, Double_t aAlfa1, Double_t aSigma1, Double_t aMu2, Double_t aAlfa2, Double_t aSigma2, Double_t aHalfWide, Double_t aSwitch, Double_t aExpScale, Double_t aExpSlope)
{
  Double_t tTemp;
  Double_t tVal, tArgSq;

  Double_t tSigSq1 = aSigma1*aSigma1;
  Double_t tAlSq1  = aAlfa1*aAlfa1;
  Double_t tAlExp1 = TMath::Exp(aAlfa1);

  Double_t tSigSq2 = aSigma2*aSigma2;
  Double_t tAlSq2  = aAlfa2*aAlfa2;
  Double_t tAlExp2 = TMath::Exp(aAlfa2);

  do {
    tTemp = (mRandom->Rndm()*50.0-25.0)*aHalfWide;
    if (tTemp > aSwitch) {
      tArgSq = (tTemp-aMu2)*(tTemp-aMu2);
      tVal  = tAlExp2*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq2)+tAlSq2));
      if ((tTemp-aMu2) > 0.0)
	tVal += aExpScale*TMath::Exp(-aExpSlope*(tTemp-aMu2));
      else
	tVal += aExpScale*TMath::Exp(+aExpSlope*(tTemp-aMu2));
    }
    else {
      tArgSq = (tTemp-aMu1)*(tTemp-aMu1);
      tVal  = tAlExp1*TMath::Exp(-TMath::Sqrt((tArgSq/tSigSq1)+tAlSq1));
      if ((tTemp-aMu1) > 0.0)
	tVal += aExpScale*TMath::Exp(-aExpSlope*(tTemp-aMu1));
      else
	tVal += aExpScale*TMath::Exp(+aExpSlope*(tTemp-aMu1));
    }
  }
  while (tVal < mRandom->Rndm()*(1.0+aExpScale));

  return tTemp;
}

inline double 
SourceModelCMSHyperExpRQMD::ScaleHyper(Double_t aScale, Double_t aX)
{
  return aX*aScale;
}


