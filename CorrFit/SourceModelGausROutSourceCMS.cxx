#include <stdlib.h>
#include <TF1.h>
#include "SourceModelGausROutSourceCMS.h"
#include <math.h>
#define SIGMAMULT 0.8
#define XSWITCH 0.5

TRandom *SourceModelGausROutSourceCMS::mRandom = new TRandom();

SourceModelGausROutSourceCMS::SourceModelGausROutSourceCMS()
{
  for (int iter=0; iter<2; iter++)
    mParameters[iter] = 0.0;
  mNPar = 5;
  mNRandVar = 4;
  mRandom->SetSeed(235431);

  TF1 *fun1 = new TF1("fun1","([0]/(sqrt(2*3.14159)*[2]))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))");
  fun1->SetParameter(0,1.0);
  fun1->SetParameter(1,0.0);
  fun1->SetParameter(2,1.0);
  Double_t bscale = fun1->Eval(XSWITCH);
  fun1->SetParameter(0,1.0/SIGMAMULT);
  fun1->SetParameter(1,0.0);
  fun1->SetParameter(2,1.0/SIGMAMULT);
  Double_t ascale = fun1->Eval(XSWITCH);
  fun1->SetParameter(0,1.0);
  fun1->SetParameter(2,1.0);
  Double_t aint = fun1->Integral(-3.0,XSWITCH);
  fun1->SetParameter(0,1.0/SIGMAMULT/(ascale/bscale));
  fun1->SetParameter(2,1.0/SIGMAMULT);
  Double_t bint = fun1->Integral(XSWITCH,3.0);
  mWideProbability = aint/(aint+bint);
  delete fun1;
}

SourceModelGausROutSourceCMS::~SourceModelGausROutSourceCMS()
{
}

void       
SourceModelGausROutSourceCMS::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausROutSourceCMS::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausROutSourceCMS::InitRandVar(double** aRandVar)
{
  Double_t tTemp;
  
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
//   for (int ti = 0; ti<mNRandVar; ti++) {
//     (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
//   }
  for (int ti = 0; ti<mNRandVar; ti++) {
    tTemp = mRandom->Rndm();
    if (tTemp > mWideProbability) {
      if (ti == 0)
	do {
	  (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0/SIGMAMULT);
	}
	while ((*aRandVar)[ti] < XSWITCH);
      else
	(*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
    }
    else {
      if (ti == 0)
	do {
	  (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
	}
	while ((*aRandVar)[ti] > XSWITCH);
      else
	(*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
    }
  }

}

void
SourceModelGausROutSourceCMS::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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
  double tRSide = aRandVar[1]*mParameters[0]/mParameters[2];
  aPos->part2.z = aRandVar[2]*mParameters[0]/mParameters[3];
  aPos->part2.t = aRandVar[3]*mParameters[0]/mParameters[4];   // =0 | Just a computing trick for the boost
  
//   tROut  *= (tMt/tM); // Rout*gammaT
//   aPos->part2.t *= (tPt/tM); // Rout*betaT*gammaT
//   double ttDTime = aPos->part2.t; 
//   aPos->part2.t += (tPz/tE*aPos->part2.z);
//   aPos->part2.t *= (tE/tMt);
//   aPos->part2.z += (tPz/tE*ttDTime); 
//   aPos->part2.z *= (tE/tMt);
  
  tPx /= tPt;
  tPy /= tPt;
  
  aPos->part2.x = tROut*tPx-tRSide*tPy;
  aPos->part2.y = tROut*tPy+tRSide*tPx;
}

const char*       
SourceModelGausROutSourceCMS::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "ROut distribution sigma";
    break;
  case 1:
    return "ROut distribution mean";
    break;
  case 2:
    return "ROut/RSide distribution sigma";
    break;
  case 3:
    return "ROut/RLong distribution sigma";
    break;
  case 4:
    return "ROut/RTime distribution sigma";
    break;
  }
}

const char* 
SourceModelGausROutSourceCMS::GetModelIdentifier()
{
  return "SMGausROutSourceCMS";
}

