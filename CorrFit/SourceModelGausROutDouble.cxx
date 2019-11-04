#include <stdlib.h>
#include <math.h>
#include "SourceModelGausROutDouble.h"
#define GAUSROUTTWONPARS 2
#define SIGMAMULT 0.8
#define XSWITCH 0.13
 
TRandom *SourceModelGausROutDouble::mRandom = new TRandom();

SourceModelGausROutDouble::SourceModelGausROutDouble()
{
  for (int iter=0; iter<GAUSROUTTWONPARS; iter++)
    mParameters[iter] = 0.0;
  mNPar = GAUSROUTTWONPARS;
  mNRandVar = 3;

  TF1 *fun1 = new TF1("fun1","[0]*exp(-((x-[1])*(x-[1]))/([2]*[2]))");
  fun1->SetParameter(0,1.0);
  fun1->SetParameter(1,0.0);
  fun1->SetParameter(2,1.0);
  Double_t aint = fun1->Integral(-10.0,XSWITCH);
  fun1->SetParameter(0,1.0);
  fun1->SetParameter(1,0.0);
  fun1->SetParameter(2,1.0/SIGMAMULT);
  Double_t bint = fun1->Integral(XSWITCH,10.0);
  mWideProbability = bint/(aint+bint);
  delete fun1;
}

SourceModelGausROutDouble::~SourceModelGausROutDouble()
{
}

void       
SourceModelGausROutDouble::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausROutDouble::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausROutDouble::InitRandVar(double** aRandVar)
{
  Double_t tTemp;
  
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti = 0; ti<mNRandVar; ti++) {
    tTemp = mRandom->Rndm();
    if (tTemp < mWideProbability) {
      if (ti == 0)
	do {
	  (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
	}
	while ((*aRandVar)[ti] > XSWITCH);
      else
	(*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
    }
    else {
      if (ti == 0)
	do {
	  (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0/SIGMAMULT);
	}
	while ((*aRandVar)[ti] < XSWITCH);
      else
	(*aRandVar)[ti] = mRandom->Gaus(0.0,1.0/SIGMAMULT);
    }
  }
}

void
SourceModelGausROutDouble::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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

  double tROut  = aRandVar[0]*mParameters[0]+mParameters[1];
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
SourceModelGausROutDouble::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "ROut distribution sigma for r*Out>mean";
    break;
  case 1:
    return "ROut distribution mean";
    break;
  }
}

const char* 
SourceModelGausROutDouble::GetModelIdentifier()
{
  return "SMGausROutDouble";
}

