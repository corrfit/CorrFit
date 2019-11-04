#include <stdlib.h>
#include "SourceModelGausROut.h"
#include <math.h>

TRandom *SourceModelGausROut::mRandom = new TRandom();

SourceModelGausROut::SourceModelGausROut()
{
  for (int iter=0; iter<2; iter++)
    mParameters[iter] = 0.0;
  mNPar = 2;
  mNRandVar = 3;
  mRandom->SetSeed(235431);
}

SourceModelGausROut::~SourceModelGausROut()
{
}

void       
SourceModelGausROut::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausROut::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausROut::InitRandVar(double** aRandVar)
{
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti = 0; ti<mNRandVar; ti++) {
    (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
  }
}

void
SourceModelGausROut::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
{
  aPos->part2.x=0.;
  aPos->part2.y=0.;
  aPos->part2.z=0.;
  aPos->part2.t=0.;
  
  double tPx = aMom->part2.x+aMom->part1.x;
  double tPy = aMom->part2.y+aMom->part1.y;
  double tPz = aMom->part2.z+aMom->part1.z;
  double tE  = aMom->part2.t+aMom->part1.t;
  double tPt = tPx*tPx + tPy*tPy;
  //mCVK = tPz*tPz;
  double tMt = tE*tE - tPz*tPz;//mCVK;
  //mCVK += tPt;
  //mCVK = sqrt(mCVK);
  double tM =   sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  
  double tROut  = aRandVar[0]*mParameters[0]+mParameters[1];
  double tRSide = aRandVar[1]*mParameters[0];
  aPos->part1.z = aRandVar[2]*mParameters[0];
  aPos->part1.t = tROut; // =0 | Just a computing trick for the boost
  
  tROut  *= (tMt/tM); // Rout*gammaT
  aPos->part1.t *= (tPt/tM); // Rout*betaT*gammaT
  double ttDTime = aPos->part1.t; 
  aPos->part1.t += (tPz/tE*aPos->part1.z);
  aPos->part1.t *= (tE/tMt);
  aPos->part1.z += (tPz/tE*ttDTime); 
  aPos->part1.z *= (tE/tMt);
  
  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part1.x = tROut*tPx-tRSide*tPy;
  aPos->part1.y = tROut*tPy+tRSide*tPx;
}

const char*       
SourceModelGausROut::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "ROut distribution sigma";
    break;
  case 1:
    return "ROut distribution mean";
    break;
  }
}

const char* 
SourceModelGausROut::GetModelIdentifier()
{
  return "SMGausROut";
}

