#include <stdlib.h>
#include <iostream>
#include <TMath.h>
#include "CFGlobal.h"
#include "SourceModelGausRStar.h"

TRandom *SourceModelGausRStar::mRandom = new TRandom();

SourceModelGausRStar::SourceModelGausRStar()
{
  for (int iter=0; iter<1; iter++)
    mParameters[iter] = 0.0;
  mNPar = 1;
  mNRandVar = 3;
  mRandom->SetSeed(235431);
}

SourceModelGausRStar::~SourceModelGausRStar()
{
}

void       
SourceModelGausRStar::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGausRStar::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGausRStar::InitRandVar(double** aRandVar)
{
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
//   for (int ti = 0; ti<mNRandVar; ti++) {
//     (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
//   }
  (*aRandVar)[0] = mRandom->Gaus(0.0,1.0);             // r* value
  (*aRandVar)[1] = mRandom->Rndm() * 2.0 * TMath::Pi();  // azimuthal angle
  (*aRandVar)[2] = mRandom->Rndm() * 2.0 - 1.0;          // cosine of the polar angle
}

void
SourceModelGausRStar::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
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
  double gammat = 1.0/TMath::Sqrt(1.0-tPt/tMt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);


//   double tROut  = aRandVar[0]*mParameters[0]+mParameters[1];
//   double tRSide = aRandVar[1]*mParameters[0];
//   aPos->part2.z = aRandVar[2]*mParameters[0];
//   aPos->part2.t = tROut; // =0 | Just a computing trick for the boost
  
//   tROut  *= (tMt/tM); // Rout*gammaT
//   aPos->part2.t *= (tPt/tM); // Rout*betaT*gammaT
//   double ttDTime = aPos->part2.t; 
//   aPos->part2.t += (tPz/tE*aPos->part2.z);
//   aPos->part2.t *= (tE/tMt);
//   aPos->part2.z += (tPz/tE*ttDTime); 
//   aPos->part2.z *= (tE/tMt);
  
//   tPx /= tPt;
//   tPy /= tPt;
  
//   aPos->part2.x = tROut*tPx-tRSide*tPy;
//   aPos->part2.y = tROut*tPy+tRSide*tPx;

  double betat  = tPt/tMt;
  double tRStar = aRandVar[0]*mParameters[0];

  double sintheta = tRStar*TMath::Sqrt(1-aRandVar[2]*aRandVar[2]);
  double tROutS  = tRStar*aRandVar[2]; // rout = r* * cos(theta)
  double tRSideS = sintheta*TMath::Sin(aRandVar[1]);

  double tRLongS = sintheta*TMath::Cos(aRandVar[1]);
  double tRTimeS = 0; // =0 | Just a computing trick for the boost
  
//   tROut         *= gammat; // Rout*gammaT
//   aPos->part2.t  = tROut*betat; // Rout*betaT*gammaT
//   double ttDTime = aPos->part2.t; 
  double betaz  = tPz/tE;
  double gammaz = 1.0/TMath::Sqrt(1.0-betaz*betaz);
//   aPos->part2.t += (betaz*aPos->part2.z);
//   aPos->part2.t *= gammaz;
//   aPos->part2.t = 0;

//   aPos->part2.z += (betaz*ttDTime); 
//   aPos->part2.z *= gammaz;

  
//   tPx /= tPt;
//   tPy /= tPt;
  
//   aPos->part2.x = tROut*tPx-tRSide*tPy;
//   aPos->part2.y = tROut*tPy+tRSide*tPx;

  double tROut = gammat * (tROutS + betat * tRTimeS);
  double tDtL  = gammat * (tRTimeS + betat * tROutS);

  double tRLong = gammaz * (tRLongS + betaz * tDtL);
  double tDt    = gammaz * (tDtL + betaz * tRLongS);

  tPx /= tPt;
  tPy /= tPt;

  aPos->part2.x = tROut*tPx - tRSideS*tPy;
  aPos->part2.y = tROut*tPy + tRSideS*tPx;
  aPos->part2.z = tRLong;
  aPos->part2.t = tDt;

//   if ((aMom->part1.x > 0.118557) && (aMom->part1.x < 0.118559) &&
//       (aMom->part2.x > 0.160845) && (aMom->part2.x < 0.160847)) {
//     PRINT_MESSAGE("Got Px1: " << aMom->part1.x);
//     PRINT_MESSAGE("Got Px2: " << aMom->part2.x);
//     PRINT_MESSAGE("Got RStar " << tRStar);
//     PRINT_MESSAGE("Got ROutS " << tROutS);
//     PRINT_MESSAGE("Got RSideS " << tRSideS);
//     PRINT_MESSAGE("Got RLongS " << tRLongS);
//     PRINT_MESSAGE("Got RTimeS " << tRTimeS);
//     PRINT_MESSAGE("Got ROut  " << aPos->part2.x);
//     PRINT_MESSAGE("Got RSide " << aPos->part2.y);
//     PRINT_MESSAGE("Got RLong " << aPos->part2.z);
//     PRINT_MESSAGE("Got RTime " << aPos->part2.t);
//   }
}

const char*       
SourceModelGausRStar::GetParameterName(int aPar)
{
  switch (aPar) {
  case 0:
    return "RStar distribution sigma";
    break;
  }
}

const char* 
SourceModelGausRStar::GetModelIdentifier()
{
  return "SMGausRStar";
}

