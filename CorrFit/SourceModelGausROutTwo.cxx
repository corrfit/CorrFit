#include <stdlib.h>
#include "SourceModelGausROutTwo.h"
#include <math.h>
#define GAUSROUTTWONPARS 3

TRandom *SourceModelGausROutTwo::mRandom = new TRandom();

SourceModelGausROutTwo::SourceModelGausROutTwo()
{
  for (int iter=0; iter<GAUSROUTTWONPARS; iter++)
    mParameters[iter] = 0.0;
  mNPar = GAUSROUTTWONPARS;
  mNRandVar = 3;
}

SourceModelGausROutTwo::~SourceModelGausROutTwo()
{
}

void       
SourceModelGausROutTwo::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausROutTwo::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausROutTwo::InitRandVar(double** aRandVar)
{
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti = 0; ti<mNRandVar; ti++) {
    (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
  }
}

void
SourceModelGausROutTwo::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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
  //mCVK = tPz*tPz;
  double tMt = tE*tE - tPz*tPz;//mCVK;
  //mCVK += tPt;
  //mCVK = sqrt(mCVK);
  double tM =   sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  
  double tROut  = aRandVar[0]>0 ? aRandVar[0]*mParameters[0]+mParameters[2] : aRandVar[0]*mParameters[1]+mParameters[2];
  double tRSide = aRandVar[1]*mParameters[0];
  aPos->part2.z = aRandVar[2]*mParameters[0];
  aPos->part2.t = tROut; // =0 | Just a computing trick for the boost
  
  tROut  *= (tMt/tM); // Rout*gammaT
  aPos->part2.t *= (tPt/tM); // Rout*betaT*gammaT
  double ttDTime = aPos->part2.t; 
  aPos->part2.t += (tPz/tE*aPos->part2.z);
  aPos->part2.t *= (tE/tMt);
  aPos->part2.z += (tPz/tE*ttDTime); 
  aPos->part2.z *= (tE/tMt);
  
  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part2.x = tROut*tPx-tRSide*tPy;
  aPos->part2.y = tROut*tPy+tRSide*tPx;
}

const char*       
SourceModelGausROutTwo::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "ROut distribution sigma for r*Out>mean";
    break;
  case 1:
    return "ROut distribution sigma for r*Out<mean";
    break;
  case 2:
    return "ROut distribution mean";
    break;
  }
}

const char* 
SourceModelGausROutTwo::GetModelIdentifier()
{
  return "SMGausROutTwo";
}

