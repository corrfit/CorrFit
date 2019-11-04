#include <stdlib.h>
#include <TF1.h>
#include "SourceModelEHRInvShift.h"

TRandom *SourceModelEHRInvShift::mRandom = new TRandom();

SourceModelEHRInvShift::SourceModelEHRInvShift()
{
  mNPar = 3;
  for (int iter=0; iter<mNPar; iter++)
    mParameters[iter] = 0.0;
  mNRandVar = 4;
}

SourceModelEHRInvShift::~SourceModelEHRInvShift()
{
}

void       
SourceModelEHRInvShift::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelEHRInvShift::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelEHRInvShift::InitRandVar(double** aRandVar)
{
  double temp;
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  do {
    (*aRandVar)[0] = GetHyper(0.0,1.0,1.0);
  }
  while (mRandom->Rndm() > ((*aRandVar)[0]*(*aRandVar)[0]*TMath::Abs((*aRandVar)[0])) / 256.0 );
  do {
    (*aRandVar)[1] = GetHyper(0.0,1.0,1.0);
  }
  while (TMath::Abs((*aRandVar)[1]) > TMath::Abs((*aRandVar)[0]));
  for (int ti = 2; ti<mNRandVar; ti++) {
    (*aRandVar)[ti] = mRandom->Rndm();
  }
}

void
SourceModelEHRInvShift::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
{
  aPos->part1.x=0.;
  aPos->part1.y=0.;
  aPos->part1.z=0.;
  aPos->part1.t=0.;
  
  double tPx = aMom->part1.x+aMom->part2.x;
  double tPy = aMom->part1.y+aMom->part2.y;
  double tPz = aMom->part1.z+aMom->part2.z;
  double tE  = aMom->part1.t+aMom->part2.t;
  double tPt = tPx*tPx + tPy*tPy;
  double tMt = tE*tE - tPz*tPz;
  double tM  = sqrt(tMt - tPt);

  double gammat = 1.0/TMath::Sqrt(1.0-tPt/tMt);

  tMt = sqrt(tMt);
  tPt = sqrt(tPt);

  double betat  = tPt/tMt;

  double betaz  = tPz/tE;
  double gammaz = 1.0/TMath::Sqrt(1.0-betaz*betaz);

  double tC1, tC2, tC3;
  double tROutS, tRSideS, tRLongS, tRTimeS;
  double tROut,  tRSide,  tRLong,  tRStar;
  double phi, theta;
  
  tRStar = mParameters[0] * aRandVar[0];
  tRTimeS = mParameters[0] * aRandVar[1];

  tC1 = TMath::Sqrt(tRStar*tRStar - tRTimeS*tRTimeS);
  phi = aRandVar[2] * TMath::Pi() * 2.0;
  if (aRandVar[3] < 0.5)
    theta = TMath::ACos(aRandVar[3] * 4.0 - 1.0);
  else
    theta = TMath::ACos((aRandVar[3]-0.5) * 4.0 - 1.0) + TMath::Pi();
  tROutS  = tC1 * TMath::Sin(theta) * TMath::Cos(phi);
  tRSideS = tC1 * TMath::Sin(theta) * TMath::Sin(phi);
  tRLongS = tC1 * TMath::Cos(theta);

  tROutS += mParameters[1];

  tROut = gammat * (tROutS + betat * tRTimeS);
  double tDtL  = gammat * (tRTimeS + betat * tROutS);

  tRLong = gammaz * (tRLongS + betaz * tDtL);
  double tDt    = gammaz * (tDtL + betaz * tRLongS);

  tPx /= tPt;
  tPy /= tPt;

  aPos->part2.x = tROut*tPx - tRSideS*tPy;
  aPos->part2.y = tROut*tPy + tRSideS*tPx;
  aPos->part2.z = tRLong;
  aPos->part2.t = tDt;
}

const char*       
SourceModelEHRInvShift::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "RInv (two-particle) distribution EH sigma";
    break;
  case 1:
    return "ROutStar (two-particle) distribution mean";
    break;
  case 2:
    return "RInv (two-particle) distribution EH alpha";
    break;
  }
}

const char* SourceModelEHRInvShift::GetModelIdentifier()
{
  return "SMEHRInvShift";
}

double 
SourceModelEHRInvShift::GetHyper(Double_t aMu, Double_t aAlfa, Double_t aSigma)
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

