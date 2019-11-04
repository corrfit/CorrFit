#include "CalcCF1DHBT.h"
#include <math.h>
#include <TMath.h>

CalcCF1DHBT::CalcCF1DHBT() : CalcCF()
{
  mNorm = 0.0;
  
  mFitLowKStar = 0.0;
  mFitHighKStar = 0.0;
  mNFitBins = 0;
  
  mNormLowKStar = 0.0;
  mNormHighKStar = 0.0;
}

CalcCF1DHBT::~CalcCF1DHBT()
{
  if (mContent) delete [] mContent;
  if (mErr2)    delete [] mErr2;
}

CalcCF1DHBT::CalcCF1DHBT(Double_t aFitLow, Double_t aFitHigh, Int_t aFitBin,
			 Double_t aNormLow, Double_t aNormHigh) : CalcCF(NULL, aFitBin)
{
  mNorm = 0.0;
  
  mFitLowKStar = aFitLow;
  mFitHighKStar = aFitHigh;
  mNFitBins = aFitBin;
  
  mNormLowKStar = aNormLow;
  mNormHighKStar = aNormHigh;

  mContent = (double *) malloc(sizeof(double) * mNFitBins);
  mErr2    = (double *) malloc(sizeof(double) * mNFitBins);
}

CalcCF1DHBT& CalcCF1DHBT::operator=(CalcCF1DHBT& aCalcCF1DHBT) 
{
  mPairManager = aCalcCF1DHBT.mPairManager;

  mNorm = aCalcCF1DHBT.mNorm;
  
  mFitLowKStar = aCalcCF1DHBT.mFitLowKStar;
  mFitHighKStar = aCalcCF1DHBT.mFitHighKStar;
  mNFitBins = aCalcCF1DHBT.mNFitBins;
  
  mNormLowKStar = aCalcCF1DHBT.mNormLowKStar;
  mNormHighKStar = aCalcCF1DHBT.mNormHighKStar;
}

short 
CalcCF1DHBT::InFitRange(Pair *aPair)
{
  return ((aPair->GetKStar() > mFitLowKStar) && (aPair->GetKStar() < mFitHighKStar));
}

short 
CalcCF1DHBT::InNormRange(Pair *aPair)
{
  return ((aPair->GetKStar() > mNormLowKStar) && (aPair->GetKStar() < mNormHighKStar));
}

int 
CalcCF1DHBT::FitRangeBin(Pair *aPair)
{
  return floor((aPair->GetKStar()-mFitLowKStar) * (mNFitBins) / (mFitHighKStar - mFitLowKStar));
}

int 
CalcCF1DHBT::NormRangeBin(Pair *aPair)
{
  return 0;
}

void CalcCF1DHBT::Write()
{
}

void CalcCF1DHBT::Normalize()
{
}

void CalcCF1DHBT::WriteParameters()
{
}

TGraphErrors *CalcCF1DHBT::GetGraph()
{
  STR tGName;

  Double_t *tX = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);
  Double_t tBinWidth = (mFitHighKStar - mFitLowKStar) / (mNFitBins);
  for (int ti=0; ti<mNFitBins; ti++) {
    tX[ti] = mFitLowKStar + ti * tBinWidth + tBinWidth / 2;
  }
  TGraphErrors *tg = new TGraphErrors(mNFitBins, tX, mContent, 0, mErr2);
  tGName = "CalcGraph";

  for (int ti=0; ti<tg->GetN(); ti++)
    tg->GetY()[ti] = TMath::Sqrt(tg->GetY()[ti]);

  tg->SetName((tGName+GetPairSystemName(mPairSystem)).Data());
  tg->SetLineColor(2);

  return tg;
}

TGraphErrors *CalcCF1DHBT::GetGraph(double aPurity)
{
  STR tGName;
  Double_t *tX    = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);
  Double_t *tCPur = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);
  Double_t *tEPur = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);

  Double_t tBinWidth = (mFitHighKStar - mFitLowKStar) / (mNFitBins);
  for (int ti=0; ti<mNFitBins; ti++) {
    tX[ti] = mFitLowKStar + ti * tBinWidth + tBinWidth / 2;
  }
  for (int ti=0; ti<mNFitBins; ti++) {
    tCPur[ti] = (mContent[ti] - 1.0)*aPurity + 1.0;
    tEPur[ti] = mErr2[ti]*aPurity;
  }

  for (int ti=0; ti<mNFitBins; ti++)
    tEPur[ti] = TMath::Sqrt(tEPur[ti]);

  TGraphErrors *tg = new TGraphErrors(mNFitBins, tX, tCPur, 0, tEPur);
  tGName = "CalcGraph";

  tg->SetName((tGName+GetPairSystemName(mPairSystem)).Data());
  tg->SetLineColor(2);
  return tg;
}

void CalcCF1DHBT::Calculate()
{
  //  if(mAlreadyCalcCF) return;

  if(mWeightCalc->getPairType()!=mPairType){
    mWeightCalc->setPairType(mPairType);
  }
  
  // Prevent rewriting on top of stored CF.
  //  mContent = mCalcContent;
  //  mErr2 = mCalcErr2;
  
  PRINT_MESSAGE("Calculating CalcCF with " << mNFitBins << "fit bins");
  
  for(int ti=0; ti<mNFitBins; ti++){
    mContent[ti]=0.;
    mErr2[ti]=0.;
  }
  mNorm = 0;
  
  Pair *tPair;
  double tNormNP=0.;
  double tWeight;
  for(int ti=0; ti<mPairManager->GetFitPairCount(); ti++){
    tPair = mPairManager->GetFitPair(ti);

    // Treat the K0s - K0s correlation specially
    // it is a mixture of K0-K0 and K0-K0b interaction.
    // So for half of the cases we need K0 - K0 weight
    // and for the other - K0 - K0bar
    if (mPairSystem == KzeroSKzeroS) {
      if (ti % 2) {
	mPairType = KaonZeroAntiKaonZero;
	mWeightCalc->setPairType(KaonZeroAntiKaonZero);
      }
      else {
	mPairType = KaonZeroKaonZero;
	mWeightCalc->setPairType(KaonZeroKaonZero);
      }
    }
    
    tPair->SetPosition();
    tWeight = mWeightCalc->getWeight(*tPair);
    if (isnan(tWeight)) {
      PRINT_MESSAGE("Fit bin weight is NaN!!! This should not happen!");
      PRINT_MESSAGE("Pair " << ti << " " << tPair);
      PRINT_MESSAGE("beta kt kstar out side long " << 
		    tPair->GetBetat() << " " << 
		    tPair->GetKt() << " " << 
		    tPair->GetKStar() << " " <<
		    tPair->GetKStarOut() << " " << 
		    tPair->GetKStarSide() << " " <<
		    tPair->GetKStarLong());
      ((StandAloneFsiKisiel *) mWeightCalc)->PrintNext();
      tWeight = mWeightCalc->getWeight(*tPair);
    }
    else {
      mContent[tPair->GetBin()] += tWeight;
      mErr2[tPair->GetBin()] += (tWeight*tWeight);

//       if (tPair->GetBin() == 9) {
// 	PRINT_MESSAGE("Pair " << tPair << " " << tWeight);
//       }
    }
  }
  for(int ti=0; ti<mPairManager->GetNormPairCount(); ti++){
    tPair = mPairManager->GetNormPair(ti);
    tPair->SetPosition();
    tWeight=(mWeightCalc->getWeight(*tPair));
    if (isnan(tWeight)) {
      PRINT_MESSAGE("Normalization weight is NaN!!! This should not happen!");
      PRINT_MESSAGE("Pair " << ti << " " << tPair);
      mNNormPair -= 1;
    }
    else {
      mNorm += tWeight;
    }
  }

  mNorm /= mNNormPair;

  double tMean;
  int tj;
  for(int ti=0; ti<mNFitBins; ti++){
    tMean = mContent[ti]/mNPairInFitBin[ti];
    mErr2[ti] = (mErr2[ti]/mNPairInFitBin[ti]-tMean*tMean)/mNPairInFitBin[ti];
    if(mErr2[ti]<0){
    }
    mContent[ti] = tMean/mNorm;
  }
}

void CalcCF1DHBT::SetNormBins(Int_t aBins) 
{
  mNNormPair = aBins;
}

int  
CalcCF1DHBT::ReadFromStorage()
{
  TVectorD *tC;
  TVectorD *tE;

  try {
#ifdef MYSQLSTORAGE
  Pair *tPair;
  SourceModel *tSourceModel;
  double *tParameters;

  tSourceModel = tPair->GetSourceModel();
  tParameters = tSourceModel->GetModelParameters();

  mStorage->SetStoreParameters(mPairNumber,
			       tSourceModel->GetModelIdentifier(),
			       tSourceModel->GetNParameters(),
			       mNFitBins,
			       mFitLowKStar,
			       mFitHighKStar,
			       mNormLowKStar,
			       mNormHighKStar,
			       Pair::GetMomRes(),
			       mPairManager->GetFirstPairHash(),
			       mPairManager->GetFitPairCount());

    mStorage->GetStoredFunction(mPairNumber, tParameters, &tC, &tE);
#else
    mStorage->GetStoredFunction(mPairNumber, GetCFIdentifier(), &tC, &tE);
#endif
    //    mStorage->GetStoredFunction(mPairNumber, GetCFIdentifier(), &tC, &tE);
    PRINT_MESSAGE("" << GetCFIdentifier() << " [Read from storage]");
  }
  catch (CF_Storage_Exception ce)
    {
      PRINT_MESSAGE("" << GetCFIdentifier() << " [Calculating]");
      return 0;
    }
  for (int ti=0; ti<tC->GetNoElements(); ti++) {
    mContent[ti] = (*tC)(ti);
    mErr2[ti]    = (*tE)(ti);
  }

  return 1;
} 
 
void 
CalcCF1DHBT::WriteToStorage()
{
  TVectorD *tC;
  TVectorD *tE;
  
  tC = new TVectorD(mNFitBins, mContent);
  tE = new TVectorD(mNFitBins, mErr2);
  
#ifdef MYSQLSTORAGE
  Pair *tPair;
  SourceModel *tSourceModel;
  double *tParameters;

  tSourceModel = tPair->GetSourceModel();
  tParameters = tSourceModel->GetModelParameters();

  mStorage->SetStoreParameters(mPairNumber,
			       tSourceModel->GetModelIdentifier(),
			       tSourceModel->GetNParameters(),
			       mNFitBins,
			       mFitLowKStar,
			       mFitHighKStar,
			       mNormLowKStar,
			       mNormHighKStar,
			       Pair::GetMomRes(),
			       mPairManager->GetFirstPairHash(),
			       mPairManager->GetFitPairCount());

  mStorage->StoreFunction(mPairNumber,tParameters, tC, tE);
#else
  mStorage->StoreFunction(mPairNumber, GetCFIdentifier(), tC, tE);
#endif
  //  mStorage->StoreFunction(mPairNumber, GetCFIdentifier(), tC, tE);

  delete tC;
  delete tE;
}

const char *
CalcCF1DHBT::GetCFIdentifier()
{
  TString mName = "";
  char   *tName;
  char   *tBuf;
  char    tFloatBuf[20];

  Pair *tPair;
  SourceModel *tSourceModel;
  double *tParameters;
  
  tName = (char *) malloc(sizeof(char) * 100);
  
  tPair = mPairManager->GetFitPair(0);
  tSourceModel = tPair->GetSourceModel();
  tParameters = tSourceModel->GetModelParameters();
  
  mName += tSourceModel->GetModelIdentifier();
  mName += " ";
  //  D_(mName.Data());
  for (int ti=0; ti<tSourceModel->GetNParameters(); ti++) {
    sprintf(tFloatBuf, "%.3lf", tParameters[ti]);
    mName += tFloatBuf;
    mName += " ";
  }
  //  D_(mName.Data());
  sprintf(tFloatBuf, "%d", mNFitBins);
  mName += tFloatBuf;
  mName += " ";
  sprintf(tFloatBuf, "%.3lf", mFitLowKStar);
  mName += tFloatBuf;
  mName += " ";
  sprintf(tFloatBuf, "%.3lf", mFitHighKStar);
  mName += tFloatBuf;
  mName += " ";
  sprintf(tFloatBuf, "%.3lf", mNormLowKStar);
  mName += tFloatBuf;
  mName += " ";
  sprintf(tFloatBuf, "%.3lf", mNormHighKStar);
  mName += tFloatBuf;
  mName += " ";
  tBuf = mPairManager->GetFirstPairHash();
  mName += tBuf;
  //  D_(tBuf);
  
  free(tBuf);
  
  //  D_(mName.Data());
  
  strcpy(tName, mName.Data());
  
  return tName;
}

