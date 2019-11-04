/***************************************************************************
 *
 * $Id: 
 *
 * Author: Laurent Conin, Fabrice Retiere, Subatech, France
 ***************************************************************************
 * 
 * Description : Create pair from HiddenInfo
 *
 ***************************************************************************
 *
 * $Log: 
 *
 ***************************************************************************/
#include <math.h>
#include "Pair.h"
#include <TRandom2.h>
#include <iostream>
#include <strstream>
#include <stdlib.h>

#include <TMath.h>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include "CFGlobal.h"

MomResParameters *Pair::mMomRes1 = 0;
MomResParameters *Pair::mMomRes2 = 0;

double Pair::mMomResScale = 0.0;

SourceModel *Pair::mSourceModel = 0;

TRandom2 *Pair::mRand = new TRandom2;

TF1 *fRndLapGaus = new TF1("fRndLapGaus", "[0]*TMath::LaplaceDist(x,0.,[1])+(1.-[0])*TMath::Gaus(x,0.,[1])", 0, 10);

double *Pair::mom1 = (double *)malloc(sizeof(double) * 4);
double *Pair::mom2 = (double *)malloc(sizeof(double) * 4);
double *Pair::pos1 = (double *)malloc(sizeof(double) * 4);
double *Pair::pos2 = (double *)malloc(sizeof(double) * 4);

int Pair::mRelocated = 0;
int Pair::mRelocateTable[8] = {0, 1, 2, 3, 4, 5, 6, 7};

// int Pair::mRelocateTable[1] = 1;
// int Pair::mRelocateTable[2] = 2;
// int Pair::mRelocateTable[3] = 3;
// int Pair::mRelocateTable[4] = 4;
// int Pair::mRelocateTable[5] = 5;
// int Pair::mRelocateTable[6] = 6;
// int Pair::mRelocateTable[7] = 7;

Pair::Pair() : mBad(0)
{
  mRandVar = 0;
  //  mRelocated=0;
  //   for (int ti=0; ti<8; ti++)
  //     mRelocateTable[ti] = ti;
};

Pair &Pair::operator=(Pair &aPair)
{
  // Only momentum is copied. Not position.
  for (int ti = 0; ti < 4; ti++)
  {
    mMom = aPair.mMom;
    //mX1[ti] = aPair.mX1[ti];
    //mX2[ti] = aPair.mX2[ti];
    mBad = aPair.mBad;
  }
  // mRandVar=0;
  //  mRelocated = aPair.mRelocated;
  //   for (int ti=0; ti<8; ti++)
  //     mRelocateTable[ti] = aPair.mRelocateTable[ti];

  mKStarSigned = aPair.mKStarSigned;

  mKStarOut = aPair.mKStarOut;
  mKStarSide = aPair.mKStarSide;
  mKStarLong = aPair.mKStarLong;

  mBin = aPair.mBin;
  return *this;
}

Pair::~Pair()
{
  //if(mMomRes1) delete mMomRes1;
  //if(mMomRes2)delete mMomRes2;
  if (mRandVar)
    delete[] mRandVar;
}

void Pair::SetMomentum(const float *aVect)
{

  // ___________________________________________________
  // --- Calculate signed kStar and bin it belongs to
  VectorPair tMom;
  tMom.part1.x = aVect[mRelocateTable[0]];
  tMom.part1.y = aVect[mRelocateTable[1]];
  tMom.part1.z = aVect[mRelocateTable[2]];
  tMom.part1.t = aVect[mRelocateTable[3]];
  mMom.part1.x = aVect[mRelocateTable[0]];
  mMom.part1.y = aVect[mRelocateTable[1]];
  mMom.part1.z = aVect[mRelocateTable[2]];
  mMom.part1.t = aVect[mRelocateTable[3]];
  tMom.part2.x = aVect[mRelocateTable[4]];
  tMom.part2.y = aVect[mRelocateTable[5]];
  tMom.part2.z = aVect[mRelocateTable[6]];
  tMom.part2.t = aVect[mRelocateTable[7]];
  mMom.part2.x = aVect[mRelocateTable[4]];
  mMom.part2.y = aVect[mRelocateTable[5]];
  mMom.part2.z = aVect[mRelocateTable[6]];
  mMom.part2.t = aVect[mRelocateTable[7]];

  //  PRINT_DEBUG("px1 " << tMom.part1.x << " py1 " << tMom.part1.y << " pz1 " << tMom.part1.z << " e1 " << tMom.part1.t);
  //  PRINT_DEBUG("px2 " << tMom.part2.x << " py2 " << tMom.part2.y << " pz2 " << tMom.part2.z << " e2 " << tMom.part2.t);

  if (mMomRes1 && mMomRes2)
  {

    if (0)
    {
      double tPhi = atan2(tMom.part1.y, tMom.part1.x);
      double tP = tMom.part1.x * tMom.part1.x + tMom.part1.y * tMom.part1.y + tMom.part1.z * tMom.part1.z;
      //    double tPt = TMath::Hypot(tMom.part1.x, tMom.part1.y);
      //    double tMassSquare  = tMom.part1.t*tMom.part1.t - tP;
      tP = sqrt(tP);
      //    double tTheta = acos(tMom.part1.z/tP);

      double per = mMomRes1->PResA + mMomRes1->PResB * pow(tP, mMomRes1->PResAlpha) + mMomRes1->PResC * tP;
      double thetaan = TMath::ATan2(hypot(tMom.part1.x, tMom.part1.y), tMom.part1.z);
      double phier = mMomRes1->PhiA + mMomRes1->PhiB * pow(tP, mMomRes1->PhiAlpha);
      double thetaer = mMomRes1->ThetaA + mMomRes1->ThetaB * pow(tP, mMomRes1->ThetaAlpha);
      double pshift = (mMomRes1->PMeanA + mMomRes1->PMeanB * pow(tP, mMomRes1->PMeanAlpha));
      double tmass = TMath::Sqrt(tMom.part1.t * tMom.part1.t - tP * tP);

      double rescale = (tP - pshift * mMomResScale) / tP;

      double Deltapx = pow(TMath::Abs(tMom.part1.x) / tP * per, 2) + pow(TMath::Abs(tMom.part1.y) * phier, 2) + pow(TMath::Abs(tMom.part1.x * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      double Deltapy = pow(TMath::Abs(tMom.part1.y) / tP * per, 2) + pow(TMath::Abs(tMom.part1.x) * phier, 2) + pow(TMath::Abs(tMom.part1.y * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      double Deltapz = pow(TMath::Abs(tMom.part1.z) / tP * per, 2) + pow(TMath::Abs(tMom.part1.z * TMath::Tan(thetaan)) * thetaer, 2);

      Deltapx = sqrt(Deltapx);
      Deltapy = sqrt(Deltapy);
      Deltapz = sqrt(Deltapz);

      mMom.part1.x = ((tMom.part1.x + mRand->Gaus(0, fabs(Deltapx) * mMomResScale)) * rescale);
      mMom.part1.y = ((tMom.part1.y + mRand->Gaus(0, fabs(Deltapy) * mMomResScale)) * rescale);
      mMom.part1.z = ((tMom.part1.z + mRand->Gaus(0, fabs(Deltapz) * mMomResScale)) * rescale);
      mMom.part1.t = TMath::Sqrt(mMom.part1.x * mMom.part1.x +
                                 mMom.part1.y * mMom.part1.y +
                                 mMom.part1.z * mMom.part1.z +
                                 tmass * tmass);

      //    double tD=0;
      //     tD = mMomRes1->PhiA+mMomRes1->PhiB*pow(tP,mMomRes1->PhiAlpha);
      //     tPhi += mRand->Gaus(0.,tD);
      //     tD = mMomRes1->ThetaA+mMomRes1->ThetaB*pow(tP,mMomRes1->ThetaAlpha);
      //     tTheta += mRand->Gaus(0.,tD);
      //     tD = tP*(mMomRes1->PResA+mMomRes1->PResB*pow(tP,mMomRes1->PResAlpha)+
      // 	     mMomRes1->PResC*tP);
      //     tP += (mMomRes1->PMeanA+mMomRes1->PMeanB*pow(tP,mMomRes1->PMeanAlpha));
      //     tP += mRand->Gaus(0.,tD);

      //     tMom.part1.x = tP*sin(tTheta)*cos(tPhi);
      //     tMom.part1.y = tP*sin(tTheta)*sin(tPhi);
      //     tMom.part1.z = tP*cos(tTheta);
      //     tMom.part1.t = sqrt(tMassSquare+tP*tP);

      //     tPhi = atan2(tMom.part2.y,tMom.part2.x);
      //     tP = tMom.part2.x*tMom.part2.x+tMom.part2.y*tMom.part2.y+tMom.part2.z*tMom.part2.z;
      //     tMassSquare  = tMom.part2.t*tMom.part2.t - tP;
      //     tP = sqrt(tP);
      //     tTheta = acos(tMom.part2.z/tP);

      //     tD=0;
      //     tD = mMomRes2->PhiA+mMomRes2->PhiB*pow(tP,mMomRes2->PhiAlpha);
      //     tPhi += mRand->Gaus(0.,tD);
      //     tD = mMomRes2->ThetaA+mMomRes2->ThetaB*pow(tP,mMomRes2->ThetaAlpha);
      //     tTheta += mRand->Gaus(0.,tD);
      //     tD = tP*(mMomRes2->PResA+mMomRes2->PResB*pow(tP,mMomRes2->PResAlpha)+
      // 	     mMomRes2->PResC*tP);
      //     tP += (mMomRes2->PMeanA+mMomRes2->PMeanB*pow(tP,mMomRes2->PMeanAlpha));
      //     tP += mRand->Gaus(0.,tD);

      //     tMom.part2.x = tP*sin(tTheta)*cos(tPhi);
      //     tMom.part2.y = tP*sin(tTheta)*sin(tPhi);
      //     tMom.part2.z = tP*cos(tTheta);
      //     tMom.part2.t = sqrt(tMassSquare+tP*tP);

      tPhi = atan2(tMom.part2.y, tMom.part2.x);
      tP = tMom.part2.x * tMom.part2.x + tMom.part2.y * tMom.part2.y + tMom.part2.z * tMom.part2.z;
      tP = sqrt(tP);

      per = mMomRes2->PResA + mMomRes2->PResB * pow(tP, mMomRes2->PResAlpha) + mMomRes2->PResC * tP;
      per *= tP;
      thetaan = TMath::ATan2(hypot(tMom.part2.x, tMom.part2.y), tMom.part2.z);
      phier = mMomRes2->PhiA + mMomRes2->PhiB * pow(tP, mMomRes2->PhiAlpha);
      thetaer = mMomRes2->ThetaA + mMomRes2->ThetaB * pow(tP, mMomRes2->ThetaAlpha);
      pshift = (mMomRes2->PMeanA + mMomRes2->PMeanB * pow(tP, mMomRes2->PMeanAlpha));
      tmass = TMath::Sqrt(tMom.part2.t * tMom.part2.t - tP * tP);

      rescale = (tP - pshift * mMomResScale) / tP;

      Deltapx = TMath::Abs(tMom.part2.x) * per + TMath::Abs(tMom.part2.y) * phier + TMath::Abs(tMom.part2.x * (1 / TMath::Tan(thetaan))) * thetaer;
      Deltapy = TMath::Abs(tMom.part2.y) * per + TMath::Abs(tMom.part2.x) * phier + TMath::Abs(tMom.part2.y * (1 / TMath::Tan(thetaan))) * thetaer;
      Deltapz = TMath::Abs(tMom.part2.z) * per + TMath::Abs(tMom.part2.z * TMath::Tan(thetaan)) * thetaer;

      mMom.part2.x = ((tMom.part2.x + mRand->Gaus(0, fabs(Deltapx) * mMomResScale)) * rescale);
      mMom.part2.y = ((tMom.part2.y + mRand->Gaus(0, fabs(Deltapy) * mMomResScale)) * rescale);
      mMom.part2.z = ((tMom.part2.z + mRand->Gaus(0, fabs(Deltapz) * mMomResScale)) * rescale);
      mMom.part2.t = TMath::Sqrt(mMom.part2.x * mMom.part2.x +
                                 mMom.part2.y * mMom.part2.y +
                                 mMom.part2.z * mMom.part2.z +
                                 tmass * tmass);
    }
    if (1)
    {
      // Reworked momentum resolution calculation for ALICE
      // with the correct parameterization

      //      double phian = atan2(tMom.part1.y,tMom.part1.x);
      double tP = tMom.part1.x * tMom.part1.x + tMom.part1.y * tMom.part1.y + tMom.part1.z * tMom.part1.z;
      tP = sqrt(tP);
      double per = mMomRes1->PResA + mMomRes1->PResB * pow(tP, mMomRes1->PResAlpha) + mMomRes1->PResC * tP;
      double thetaan = TMath::ATan2(hypot(tMom.part1.x, tMom.part1.y), tMom.part1.z);
      double phier = 0.0;
      if (fabs(mMomRes1->PhiMu) < 0.00000001)
        phier = mMomRes1->PhiA + mMomRes1->PhiB * pow(tP, mMomRes1->PhiAlpha);
      else
        phier = mMomRes1->PhiA + mMomRes1->PhiB * exp(mMomRes1->PhiAlpha * (tP - mMomRes1->PhiMu));
      double thetaer = mMomRes1->ThetaA + mMomRes1->ThetaB * pow(tP, mMomRes1->ThetaAlpha);

      double tmass = TMath::Sqrt(tMom.part1.t * tMom.part1.t - tP * tP);

      //      double pshift = (mMomRes1->PMeanA+mMomRes1->PMeanB*pow(tP,mMomRes1->PMeanAlpha));
      // New Energy loss from 0808.2041
      double pshift = (mMomRes1->PMeanA + mMomRes1->PMeanB * pow(1.0 + tmass * tmass / (tP * tP), mMomRes1->PMeanAlpha));
      double rescale = (tP - pshift * mMomResScale) / tP;

      double Deltapx = pow(TMath::Abs(tMom.part1.x) / tP * per, 2) + pow(TMath::Abs(tMom.part1.y) * phier, 2) + pow(TMath::Abs(tMom.part1.x * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      double Deltapy = pow(TMath::Abs(tMom.part1.y) / tP * per, 2) + pow(TMath::Abs(tMom.part1.x) * phier, 2) + pow(TMath::Abs(tMom.part1.y * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      double Deltapz = pow(TMath::Abs(tMom.part1.z) / tP * per, 2) + pow(TMath::Abs(tMom.part1.z * TMath::Tan(thetaan)) * thetaer, 2);

      Deltapx = sqrt(Deltapx);
      Deltapy = sqrt(Deltapy);
      Deltapz = sqrt(Deltapz);

      if (mRand->Rndm() > 0.4)
      {
        mMom.part1.x = ((tMom.part1.x + mRand->Gaus(0, fabs(Deltapx) * mMomResScale)) * rescale);
        mMom.part1.y = ((tMom.part1.y + mRand->Gaus(0, fabs(Deltapy) * mMomResScale)) * rescale);
        mMom.part1.z = ((tMom.part1.z + mRand->Gaus(0, fabs(Deltapz) * mMomResScale)) * rescale);
      }
      else
      {
        double sigx = mRand->Rndm() > 0.5 ? 1. : -1.;
        double sigy = mRand->Rndm() > 0.5 ? 1. : -1.;
        double sigz = mRand->Rndm() > 0.5 ? 1. : -1.;
        mMom.part1.x = ((tMom.part1.x + sigx * mRand->Exp((Deltapx)*mMomResScale)) * rescale);
        mMom.part1.y = ((tMom.part1.y + sigy * mRand->Exp((Deltapy)*mMomResScale)) * rescale);
        mMom.part1.z = ((tMom.part1.z + sigz * mRand->Exp((Deltapz)*mMomResScale)) * rescale);
      }

      mMom.part1.t = TMath::Sqrt(mMom.part1.x * mMom.part1.x +
                                 mMom.part1.y * mMom.part1.y +
                                 mMom.part1.z * mMom.part1.z +
                                 tmass * tmass);

      //      phian = atan2(tMom.part2.y,tMom.part2.x);
      tP = tMom.part2.x * tMom.part2.x + tMom.part2.y * tMom.part2.y + tMom.part2.z * tMom.part2.z;
      tP = sqrt(tP);

      per = mMomRes2->PResA + mMomRes2->PResB * pow(tP, mMomRes2->PResAlpha) + mMomRes2->PResC * tP;
      thetaan = TMath::ATan2(hypot(tMom.part2.x, tMom.part2.y), tMom.part2.z);
      if (fabs(mMomRes2->PhiMu) < 0.00000001)
        phier = mMomRes2->PhiA + mMomRes2->PhiB * pow(tP, mMomRes2->PhiAlpha);
      else
        phier = mMomRes2->PhiA + mMomRes2->PhiB * exp(mMomRes2->PhiAlpha * (tP - mMomRes2->PhiMu));
      thetaer = mMomRes2->ThetaA + mMomRes2->ThetaB * pow(tP, mMomRes2->ThetaAlpha);
      tmass = TMath::Sqrt(tMom.part2.t * tMom.part2.t - tP * tP);

      //      pshift = (mMomRes2->PMeanA+mMomRes2->PMeanB*pow(tP,mMomRes2->PMeanAlpha));
      // New Energy loss from 0808.2041
      pshift = (mMomRes2->PMeanA + mMomRes2->PMeanB * pow(1.0 + tmass * tmass / (tP * tP), mMomRes2->PMeanAlpha));
      rescale = (tP - pshift * mMomResScale) / tP;

      Deltapx = pow(TMath::Abs(tMom.part2.x) / tP * per, 2) + pow(TMath::Abs(tMom.part2.y) * phier, 2) + pow(TMath::Abs(tMom.part2.x * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      Deltapy = pow(TMath::Abs(tMom.part2.y) / tP * per, 2) + pow(TMath::Abs(tMom.part2.x) * phier, 2) + pow(TMath::Abs(tMom.part2.y * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      Deltapz = pow(TMath::Abs(tMom.part2.z) / tP * per, 2) + pow(TMath::Abs(tMom.part2.z * TMath::Tan(thetaan)) * thetaer, 2);

      Deltapx = sqrt(Deltapx);
      Deltapy = sqrt(Deltapy);
      Deltapz = sqrt(Deltapz);

      // fraction of gaus and laplace:
      // Momentum: 0.4,
      // Phi = 0.4
      // Theta = 0.25

      if (mRand->Rndm() > 0.4)
      {
        mMom.part2.x = ((tMom.part2.x + mRand->Gaus(0, fabs(Deltapx) * mMomResScale)) * rescale);
        mMom.part2.y = ((tMom.part2.y + mRand->Gaus(0, fabs(Deltapy) * mMomResScale)) * rescale);
        mMom.part2.z = ((tMom.part2.z + mRand->Gaus(0, fabs(Deltapz) * mMomResScale)) * rescale);
      }
      else
      {
        double sigx = mRand->Rndm() > 0.5 ? 1. : -1.;
        double sigy = mRand->Rndm() > 0.5 ? 1. : -1.;
        double sigz = mRand->Rndm() > 0.5 ? 1. : -1.;
        mMom.part2.x = ((tMom.part2.x + sigx * mRand->Exp((Deltapx)*mMomResScale)) * rescale);
        mMom.part2.y = ((tMom.part2.y + sigy * mRand->Exp((Deltapy)*mMomResScale)) * rescale);
        mMom.part2.z = ((tMom.part2.z + sigz * mRand->Exp((Deltapz)*mMomResScale)) * rescale);
      }

      mMom.part2.t = TMath::Sqrt(mMom.part2.x * mMom.part2.x +
                                 mMom.part2.y * mMom.part2.y +
                                 mMom.part2.z * mMom.part2.z +
                                 tmass * tmass);
    }
  }

  double tPx = tMom.part1.x + tMom.part2.x;
  double tPy = tMom.part1.y + tMom.part2.y;
  double tPz = tMom.part1.z + tMom.part2.z;
  double tE = tMom.part1.t + tMom.part2.t;
  double tPt = tPx * tPx + tPy * tPy;
  //mCVK = tPz*tPz;
  double tMt = tE * tE - tPz * tPz; //mCVK;
  //mCVK += tPt;
  //mCVK = sqrt(mCVK);
  double tM = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);

  // Boost to LCMS
  double tBeta = tPz / tE;
  double tGamma = tE / tMt;
  mKStarLong = tGamma * (tMom.part1.z - tBeta * tMom.part1.t);
  double tE1L = tGamma * (tMom.part1.t - tBeta * tMom.part1.z);

  // Rotate in transverse plane
  mKStarOut = (tMom.part1.x * tPx + tMom.part1.y * tPy) / tPt;
  mKStarSide = (-tMom.part1.x * tPy + tMom.part1.y * tPx) / tPt;

  // Boost to pair cms
  mKStarOut = tMt / tM * (mKStarOut - tPt / tMt * tE1L);

  mKStarSigned = mKStarOut > 0. ? 1. : -1.;
  mKStarSigned *= sqrt(mKStarSide * mKStarSide + mKStarOut * mKStarOut + mKStarLong * mKStarLong);

  //mCVK = (mDKOut*tPt + mDKLong*tPz)/kStarCalc/mCVK;
  //mDKLong = (tPz>=0)* mDKLong - (tPz<0)*mDKLong;
}

void Pair::Set(const float *aVect)
{
  SetMomentum(aVect);
}

void Pair::SetPosition()
{
  if (!mRandVar)
    InitRandVar();
  mSourceModel->GeneratePos(&mMom, &mPos, mRandVar);
}

void Pair::InitRandVar()
{
  mSourceModel->InitRandVar(&mRandVar);
  //  PRINT_DEBUG_3("Initialized RandVar to: " << mRandVar[0] << " " << mRandVar[1] << " " << mRandVar[2] << " " << mRandVar[3]);
}

double Pair::GetBetat()
{
  double tPx = mMom.part1.x + mMom.part2.x;
  double tPy = mMom.part1.y + mMom.part2.y;
  double tPz = mMom.part1.z + mMom.part2.z;
  double tE = mMom.part1.t + mMom.part2.t;
  return sqrt((tPx * tPx + tPy * tPy) / (tE * tE - tPz * tPz));
}
double Pair::GetKt()
{
  double tPx = mMom.part1.x + mMom.part2.x;
  double tPy = mMom.part1.y + mMom.part2.y;
  double tPz = mMom.part1.z + mMom.part2.z;
  return sqrt(tPx * tPx + tPy * tPy);
}

void Pair::SetMomRes(double aMomResScale, int aPairSystem)
{
  mMomResScale = aMomResScale;

  switch (aPairSystem)
  {
  case 1:
  case 2:
  case 3:
  case 4:
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetKaonMomResParameters();
    break;
  case 11:
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
    break;
  case 12:
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetAntiProtonMomResParameters();
    break;
  case 13:
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
    break;
  case 14:
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetAntiProtonMomResParameters();
    break;
  case 21:
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
    break;
  case 22:
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetAntiProtonMomResParameters();
    break;
  case 23:
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
    break;
  case 24:
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetAntiProtonMomResParameters();
    break;
  case 31:
  case 32:
  case 33:
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetPionMomResParameters();
    break;
  case 41:
  case 42:
  case 43:
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetKaonMomResParameters();
    break;
  case 51:
    mMomRes1 = GetProtonMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
    break;
  }
}

double Pair::GetMomRes()
{
  return mMomResScale;
}

MomResParameters *Pair::GetPionMomResParameters()
{
  MomResParameters *t = new MomResParameters;
  //   t->PhiA = 0.00204489;
  //   t->PhiB = 0.00100642;
  //   t->PhiAlpha = -1.28315;
  //   t->ThetaA = 0.000930104;
  //   t->ThetaB = 0.00123727;
  //   t->ThetaAlpha = -1.14998;
  //   t->PResA = 0.0113174;
  //   t->PResB = 0.00168459;
  //   t->PResAlpha = -1.04594;
  //   t->PResC = 0.00948437;
  //   t->PMeanA = 0.; // no energy loss (Kalman filter)
  //   t->PMeanB = 0.;
  //   t->PMeanAlpha = 0.;
  /* 
  t->PhiA = 0.00201;
  t->PhiB = 0.001018;
  t->PhiAlpha = -1.274;
  t->PhiMu = 0.0;
  t->ThetaA = 0.000908;
  t->ThetaB = 0.001255;
  t->ThetaAlpha = -1.141;
  t->PResA = 0.01074;
  t->PResB = 0.001918;
  t->PResAlpha = -0.9895;
  t->PResC = 0.009706;
  t->PMeanA = 0.; // no energy loss (Kalman filter)
  t->PMeanB = 0.;
  t->PMeanAlpha = 0.;
  */

  // New parameters for ALICE PbPb 2.7 TeV

  t->PResA = -0.0235497;
  t->PResB = 0.0175338;
  t->PResC = 0.0145192;
  t->PResAlpha = -0.171772;
  t->ThetaA = 6.56387e-05;
  t->ThetaB = 0.00382418;
  t->ThetaAlpha = -0.975946;
  t->PhiA = 0.000912177;
  t->PhiB = 0.00253254;
  t->PhiAlpha = -1.1757;
  t->PhiMu = 0.;
  t->PMeanA = 0.;
  t->PMeanB = 0.;
  t->PMeanAlpha = 0.;

  return t;
}

void Pair::SetRelocate(int number, const char *name)
{
  if (!strcmp(name, "px1"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as px1 0\n");
    mRelocateTable[0] = number;
    mRelocated &= (1 >> 0);
  }
  else if (!strcmp(name, "py1"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as py1 1\n");
    mRelocateTable[1] = number;
    mRelocated &= (1 >> 1);
  }
  else if (!strcmp(name, "pz1"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as pz1 2\n");
    mRelocateTable[2] = number;
    mRelocated &= (1 >> 2);
  }
  else if (!strcmp(name, "e1"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as e1 3\n");
    mRelocateTable[3] = number;
    mRelocated &= (1 >> 3);
  }
  else if (!strcmp(name, "px2"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as px2 4\n");
    mRelocateTable[4] = number;
    mRelocated &= (1 >> 4);
  }
  else if (!strcmp(name, "py2"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as py2 5\n");
    mRelocateTable[5] = number;
    mRelocated &= (1 >> 5);
  }
  else if (!strcmp(name, "pz2"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as pz2 6\n");
    mRelocateTable[6] = number;
    mRelocated &= (1 >> 6);
  }
  else if (!strcmp(name, "e2"))
  {
    PRINT_MESSAGE("Relocating " << name << " " << number << " number as e2 7\n");
    mRelocateTable[7] = number;
    mRelocated &= (1 >> 7);
  }
  else
  {
    PRINT_MESSAGE("Unknown branch");
  }
}

MomResParameters *Pair::GetKaonMomResParameters()
{
  MomResParameters *t = new MomResParameters;
  //   t->PResA = 0.8015;
  //   t->PResB = -2.362;
  //   t->PResC = 1.846;
  //   t->PResAlpha = 0.3445;
  //   t->PhiA = 0.0006575;
  //   t->PhiB = 0.002813;
  //   t->PhiAlpha = -1.583;
  //   t->ThetaA = 0.0002846;
  //   t->ThetaB = 0.002458;
  //   t->ThetaAlpha = -1.475;
  //   t->PMeanA = -0.006509;
  //   t->PMeanB = 0.008757;
  //   t->PMeanAlpha = -1.373;
  /*  t->PResA = 0.01981;
  t->PResB = 0.001371;
  t->PResC = 0.0;
  t->PResAlpha = -2.112;

  t->PhiA = 0.001791;
  t->PhiB = 0.001319;
  t->PhiAlpha = -1.686;
  t->PhiMu = 0.0;

  t->ThetaA = 0.0005202;
  t->ThetaB = 0.001752;
  t->ThetaAlpha = -1.352;

  //   t->PMeanA = -0.004136;
  //   t->PMeanB = 0.003511;
  //   t->PMeanAlpha = -1.192;

  // New parametrization of energy loss from STAR 0808.2041
  t->PMeanA = 0.006;
  t->PMeanB = -0.0038;
  t->PMeanAlpha = 1.10;
*/
  // New parameters for ALICE PbPb 2.7 TeV
  t->PResA = -0.00827883;
  t->PResB = 0.00158054;
  t->PResC = 0.0150525;
  t->PResAlpha = -1.58543;
  t->ThetaA = 0.000979426;
  t->ThetaB = 0.0032204;
  t->ThetaAlpha = -1.49797;
  t->PhiA = 0.00174143;
  t->PhiB = 0.00184737;
  t->PhiAlpha = -1.90668;
  t->PhiMu = 0.;
  t->PMeanA = 0.;
  t->PMeanB = 0.;
  t->PMeanAlpha = 0.;

  return t;
}

MomResParameters *Pair::GetProtonMomResParameters()
{
  MomResParameters *t = new MomResParameters;
  //   t->PhiA = 0.00179782;
  //   t->PhiB = 0.00133481;
  //   t->PhiAlpha = -1.68111;
  //   t->ThetaA = 0.000495737;
  //   t->ThetaB = 0.0017755;
  //   t->ThetaAlpha = -1.34328;
  //   t->PResA = 0.0198917;
  //   t->PResB = 0.00145757;
  //   t->PResAlpha = -2.07794;
  //   t->PResC = 0;
  //   t->PMeanA = -0.00419502;
  //   t->PMeanB = 0.00353296;
  //   t->PMeanAlpha = -1.19465;

  //   t->PResA = 0.7989;
  //   t->PResB = -2.128;
  //   t->PResAlpha = 0.5282;
  //   t->PResC = 1.41;

  // Year 1 parameterization

  //   t->PResA = 0.01708;
  //   t->PResB = 0.006794;
  //   t->PResC = 0.0;
  //   t->PResAlpha = -1.78;

  //   t->PhiA = 0.0006575;
  //   t->PhiB = 0.002813;
  //   t->PhiAlpha = -1.583;

  //   t->ThetaA = 0.0002846;
  //   t->ThetaB = 0.002458;
  //   t->ThetaAlpha = -1.475;

  //   t->PMeanA = -0.006509;
  //   t->PMeanB = 0.008757;
  //   t->PMeanAlpha = -1.373;

  // Year 2 parameterization

  //   t->PMeanA = -0.00905;
  //   t->PMeanB = 0.007839;
  //   t->PMeanAlpha = -1.69;

  // New parametrization of energy loss from STAR 0808.2041
  //t->PMeanA = 0.013;
  //t->PMeanB = -0.0081;
  //t->PMeanAlpha = 1.03;

  // Alice 2.7 TeV resolution with laplace
  t->PResA = -0.00538505;
  t->PResB = 0.000567673;
  t->PResC = 0.0139363;
  t->PResAlpha = -3.62941;
  t->ThetaA = 0.00113763;
  t->ThetaB = 0.0041515;
  t->ThetaAlpha = -1.79528;
  t->PhiA = 0.00208896;
  t->PhiB = 0.00292743;
  t->PhiAlpha = -2.30148;
  t->PhiMu = 0.;
  t->PMeanA = 0.;
  t->PMeanB = 0.;
  t->PMeanAlpha = 0.;

  return t;
}

MomResParameters *Pair::GetAntiProtonMomResParameters()
{
  MomResParameters *t = new MomResParameters;
  // Year 2 parameterization
  t->PResA = -0.01857;
  t->PResB = 0.002276;
  t->PResC = 0.05261;
  t->PResAlpha = -2.545;

  t->PhiA = 0.008614;
  t->PhiB = 0.0174;
  t->PhiAlpha = -5.717;
  t->PhiMu = 0.3319;

  t->ThetaA = 0.002326;
  t->ThetaB = 0.00596;
  t->ThetaAlpha = -1.052;

  t->PMeanA = -0.006365;
  t->PMeanB = 0.006887;
  t->PMeanAlpha = -1.747;

  return t;
}

void Pair::SetSourceModel(SourceModel *aSourceModel) { mSourceModel = aSourceModel; }

SourceModel *Pair::GetSourceModel() { return mSourceModel; }

void Pair::SetPosMom()
{
  mom1[0] = mMom.part1.x;
  mom1[1] = mMom.part1.y;
  mom1[2] = mMom.part1.z;
  mom1[3] = mMom.part1.t;

  mom2[0] = mMom.part2.x;
  mom2[1] = mMom.part2.y;
  mom2[2] = mMom.part2.z;
  mom2[3] = mMom.part2.t;

  pos1[0] = mPos.part1.x;
  pos1[1] = mPos.part1.y;
  pos1[2] = mPos.part1.z;
  pos1[3] = mPos.part1.t;

  pos2[0] = mPos.part2.x;
  pos2[1] = mPos.part2.y;
  pos2[2] = mPos.part2.z;
  pos2[3] = mPos.part2.t;
}

void Pair::SetPosMomEFirstNotation()
{
  mom1[1] = mMom.part1.x * 1000.0;
  mom1[2] = mMom.part1.y * 1000.0;
  mom1[3] = mMom.part1.z * 1000.0;
  mom1[0] = mMom.part1.t * 1000.0;

  mom2[1] = mMom.part2.x * 1000.0;
  mom2[2] = mMom.part2.y * 1000.0;
  mom2[3] = mMom.part2.z * 1000.0;
  mom2[0] = mMom.part2.t * 1000.0;

  pos1[1] = mPos.part1.x;
  pos1[2] = mPos.part1.y;
  pos1[3] = mPos.part1.z;
  pos1[0] = mPos.part1.t;

  pos2[1] = mPos.part2.x;
  pos2[2] = mPos.part2.y;
  pos2[3] = mPos.part2.z;
  pos2[0] = mPos.part2.t;
}

void Pair::SetPairNum(int aNum)
{
  mPairNum = aNum;
}

int Pair::GetPairNum()
{
  return mPairNum;
}

void Pair::GetSeed(UInt_t &s1, UInt_t &s2)
{
  s1 = mRand->GetSeed();
}

void Pair::SetSeed(UInt_t s1, UInt_t s2)
{
  mRand->SetSeed(s1);
}
