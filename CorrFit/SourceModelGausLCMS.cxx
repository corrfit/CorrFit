#include <stdlib.h>
#include <TF1.h>
#include "SourceModelGausLCMS.h"

#define PAROUT4 1.0
#define PAROUT5 0.0
#define PAROUT7 1.0
#define XSWITCH 0.0

TRandom *SourceModelGausLCMS::mRandom = new TRandom();

SourceModelGausLCMS::SourceModelGausLCMS()
{
  for (int iter=0; iter<5; iter++)
    mParameters[iter] = 0.0;
  mNPar = 5;
  mNRandVar = 4;
}

SourceModelGausLCMS::~SourceModelGausLCMS()
{
}

void       
SourceModelGausLCMS::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausLCMS::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausLCMS::InitRandVar(double** aRandVar)
{
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti = 0; ti<mNRandVar; ti++) {
    (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
  }  
}

void
SourceModelGausLCMS::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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

  double betaz  = tPz/tE;
  double gammaz = 1.0/TMath::Sqrt(1.0-betaz*betaz);

  double tROut;
  tROut  = mParameters[0] * aRandVar[0] + mParameters[1];
  
  double tRLongS;
  double tDtL;

  double tRSide = mParameters[0] * mParameters[2] * aRandVar[1];
  tRLongS = mParameters[0] * mParameters[3] * aRandVar[2];
  tDtL = 0.0;
  
  double tRLong, tDt;
  tRLong = gammaz * (tRLongS + betaz * tDtL);
  tDt    = gammaz * (tDtL + betaz * tRLongS);

  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part1.x = tROut*tPx-tRSide*tPy;
  aPos->part1.y = tROut*tPy+tRSide*tPx;
  aPos->part1.z = tRLong;
  aPos->part1.t = tDt;
}

const char*       
SourceModelGausLCMS::GetParameterName(int aPar)
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

const char* SourceModelGausLCMS::GetModelIdentifier()
{
  return "SMGausROutLCMS";
}


