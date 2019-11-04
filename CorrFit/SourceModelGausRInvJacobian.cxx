#include <stdlib.h>
#include <TF1.h>
#include "SourceModelGausRInvJacobian.h"

TRandom *SourceModelGausRInvJacobian::mRandom = new TRandom();

SourceModelGausRInvJacobian::SourceModelGausRInvJacobian()
{
  for (int iter=0; iter<1; iter++)
    mParameters[iter] = 0.0;
  mNPar = 1;
  mNRandVar = 7;
}

SourceModelGausRInvJacobian::~SourceModelGausRInvJacobian()
{
}

void       
SourceModelGausRInvJacobian::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausRInvJacobian::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausRInvJacobian::InitRandVar(double** aRandVar)
{
  double temp;
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  do {
    (*aRandVar[0]) = mRandom->Gaus(0.0,1.0);
  }
  while (mRandom->Rndm() > ((*aRandVar[0])*(*aRandVar[0])) / 9.0 );
  for (int ti = 1; ti<mNRandVar; ti++) {
    (*aRandVar)[ti] = mRandom->Rndm();
  }
}

void
SourceModelGausRInvJacobian::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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

  tRStar = mParameters[0] * aRandVar[0];

  tC1 = aRandVar[1] * tRStar;
  tC2 = aRandVar[2] * (TMath::Sqrt(tRStar*tRStar - tC1 * tC1));
  tC3 = TMath::Sqrt(tRStar*tRStar - tC1 * tC1 - tC2 * tC2);

  if (aRandVar[4] > 0.5) tC1 *= -1.0;
  if (aRandVar[5] > 0.5) tC2 *= -1.0;
  if (aRandVar[6] > 0.5) tC3 *= -1.0;

  if (aRandVar[3] < (1.0/6.0))
    {
      tROutS  = tC1;
      tRSideS = tC2;
      tRLongS = tC3;
    }
  else if (aRandVar[3] < (2.0/6.0))
    {
      tROutS  = tC1;
      tRSideS = tC3;
      tRLongS = tC2;
    }
  else if (aRandVar[3] < (3.0/6.0))
    {
      tROutS  = tC2;
      tRSideS = tC1;
      tRLongS = tC3;
    }
  else if (aRandVar[3] < (4.0/6.0))
    {
      tROutS  = tC2;
      tRSideS = tC3;
      tRLongS = tC1;
    }
  else if (aRandVar[3] < (5.0/6.0))
    {
      tROutS  = tC3;
      tRSideS = tC1;
      tRLongS = tC2;
    }
  else 
    {
      tROutS  = tC3;
      tRSideS = tC2;
      tRLongS = tC1;
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
SourceModelGausRInvJacobian::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "RInv (two-particle) distribution sigma";
    break;
  }
}

const char* SourceModelGausRInvJacobian::GetModelIdentifier()
{
  return "SMGausRInvGausJacob";
}


