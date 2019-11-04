#include <stdlib.h>
#include <TF1.h>
#include "SourceModelGausCMS.h"

#define PAROUT4 1.0
#define PAROUT5 0.0
#define PAROUT7 1.0
#define XSWITCH 0.0

TRandom *SourceModelGausCMS::mRandom = new TRandom();

SourceModelGausCMS::SourceModelGausCMS()
{
  for (int iter=0; iter<5; iter++)
    mParameters[iter] = 0.0;
  mNPar = 5;
  mNRandVar = 4;

  TF1 *fun1 = new TF1("fun1","(1.0/[2])*[0]*exp(-((x-[1])*(x-[1]))/([2]*[2]))");
  fun1->SetParameter(0,1.0);
  fun1->SetParameter(1,0.0);
  fun1->SetParameter(2,1.0);
  Double_t aint = fun1->Integral(-10.0,XSWITCH);
  fun1->SetParameter(0,1.0);
  fun1->SetParameter(1,0.0);
  fun1->SetParameter(2,1.0*PAROUT7);
  Double_t bint = fun1->Integral(XSWITCH,10.0);
  mNegProbability = aint/(aint+bint);
  delete fun1;
}

SourceModelGausCMS::~SourceModelGausCMS()
{
}

void       
SourceModelGausCMS::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausCMS::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausCMS::InitRandVar(double** aRandVar)
{
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti = 1; ti<mNRandVar; ti++) {
    (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
  }
  if (mRandom->Rndm() > mNegProbability) 
    {
      (*aRandVar)[0] = TMath::Abs(mRandom->Gaus()*PAROUT7);
    }
  else
    {
      (*aRandVar)[0] = -TMath::Abs(mRandom->Gaus());
    }
  
}

void
SourceModelGausCMS::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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
  double tMt = tE*tE - tPz*tPz;
  double tM  = sqrt(tMt - tPt);

  tMt = sqrt(tMt);
  tPt = sqrt(tPt);

  double tROut;
  if (aRandVar[0] > 0.0)
    tROut  = mParameters[0] * aRandVar[0] + mParameters[1];
  else
    tROut  = mParameters[0] * aRandVar[0] * PAROUT7 + mParameters[1];
  
  double tRSide = mParameters[0] * mParameters[2] * aRandVar[1];
  aPos->part1.z = mParameters[0] * mParameters[3] * aRandVar[2];
  aPos->part1.t = mParameters[0] * mParameters[4] * aRandVar[3];
  
  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part1.x = tROut*tPx-tRSide*tPy;
  aPos->part1.y = tROut*tPy+tRSide*tPx;
}

const char*       
SourceModelGausCMS::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "ROut distribution sigma";
    break;
  case 1:
    return "ROut distribution mean";
    break;
  case 2:
    return "RSide distribution sigma multiplier";
    break;
  case 3:
    return "RLong distribution sigma multiplier";
    break;
  case 4:
    return "RTime distribution sigma multiplier";
    break;
  }
}

const char* SourceModelGausCMS::GetModelIdentifier()
{
  return "SMGausROutCMS";
}


