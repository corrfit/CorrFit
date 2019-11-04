#include <stdlib.h>
#include <TF1.h>
#include "SourceModelExpGausRInv.h"

TRandom *SourceModelExpGausRInv::mRandom = new TRandom();

SourceModelExpGausRInv::SourceModelExpGausRInv()
{
  for (int iter=0; iter<3; iter++)
    mParameters[iter] = 0.0;
  mNPar = 3;
  mNRandVar = 7;
  mExpProbability = 0.0;
}

SourceModelExpGausRInv::~SourceModelExpGausRInv()
{
}

void       
SourceModelExpGausRInv::SetModelParameters(double *aParTable)
{
  Double_t tEta = 0;
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
//   tEta = mParameters[1]/(1.0-mParameters[1]);
//   mExpProbability = tEta*8*TMath::Pi()/TMath::Power(mParameters[2],3) / 
//     (tEta*8*TMath::Pi()/TMath::Power(mParameters[2],3) +
//      8*TMath::Power(TMath::Pi(), 1.5)*TMath::Power(mParameters[0],3));
  mExpProbability = mParameters[1];
}

double*
SourceModelExpGausRInv::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelExpGausRInv::InitRandVar(double** aRandVar)
{
  double temp;
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  (*aRandVar)[0] = mRandom->Rndm();
  (*aRandVar)[1] = mRandom->Gaus(0.0,1.0);
  (*aRandVar)[2] = mRandom->Gaus(0.0,1.0);
  (*aRandVar)[3] = mRandom->Gaus(0.0,1.0);
  (*aRandVar)[4] = mRandom->Rndm()*2.0;
  (*aRandVar)[5] = mRandom->Rndm()*2.0 - 1.0;
  (*aRandVar)[6] = mRandom->Rndm()*TMath::Pi()*2.0;
}

void
SourceModelExpGausRInv::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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

  if (aRandVar[0] > mExpProbability) {
    // We generate from Gaus
    tROutS = mParameters[0] * aRandVar[1];
    tRSideS = mParameters[0] * aRandVar[2];
    tRLongS = mParameters[0] * aRandVar[3];
  }
  else {
    tRStar  = GetXFromY(aRandVar[4], mParameters[2]);
    tROutS  = tRStar*TMath::Sqrt(1-aRandVar[5]*aRandVar[5]);
    tRSideS = tROutS*TMath::Sin(aRandVar[6]);
    tROutS  = tROutS*TMath::Cos(aRandVar[6]);
    tRLongS = tRStar*aRandVar[5];
  }
  tRTimeS = 0.0;

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
SourceModelExpGausRInv::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "RInv (two-particle) Gaus distribution sigma";
    break;
  case 1:
    return "Zeta parameter - Exp to Gaus ratio";
    break;
  case 2:
    return "RInv (two-particle) Exp distribution slope";
    break;
  }
}

const char* SourceModelExpGausRInv::GetModelIdentifier()
{
  return "SMGausRInvExpGaus";
}


double 
SourceModelExpGausRInv::GetXFromY(double aY, double aPar)
{
  double tPrecision = 1.0e-4;

  double tX = 2.5/aPar;
  double tStep = tX/2.0;
  double tY = FunExp(tX, aPar);
  int tWasUp = -1;

  while (fabs(tY - aY) > tPrecision) {
    if (tY > aY) {
      tX += tStep;
      if (tWasUp == 0)
	tStep *= 0.5;
      tWasUp = 1;
      if (tX < 0) tX=0;
    }
    else {
      tX -= tStep;
      if (tWasUp == 1)
	tStep *= 0.5;
      tWasUp = 0;
    }
    tY = FunExp(tX, aPar);
  }
  
  return tX;
}

double 
SourceModelExpGausRInv::FunExp(double aArg, double aPar)
{
  return exp(-aPar*aArg)*(aPar*aPar*aArg*aArg + 2.0*aArg*aPar + 2);
}
