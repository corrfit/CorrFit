#include "CalcCFNonId.h"

CalcCFNonId::CalcCFNonId() : CalcCF()
{
  mNormP = 0.0;
  mNormN = 0.0;
  
  mFitLowKStar = 0.0;
  mFitHighKStar = 0.0;
  mNFitBins = 0;
  
  mNormLowKStar = 0.0;
  mNormHighKStar = 0.0;

  mBinWidth = 0.0;
}

CalcCFNonId::~CalcCFNonId()
{
  if (mContent) delete [] mContent;
  if (mErr2)    delete [] mErr2;
}

CalcCFNonId::CalcCFNonId(Double_t aFitLow, Double_t aFitHigh, Int_t aFitBin,
			 Double_t aNormLow, Double_t aNormHigh) : CalcCF(NULL, aFitBin)
{
  mNormP = 0.0;
  mNormN = 0.0;
  
  mFitLowKStar = aFitLow;
  mFitHighKStar = aFitHigh;
  mNFitBins = aFitBin;
  
  mNormLowKStar = aNormLow;
  mNormHighKStar = aNormHigh;

  mBinWidth = 2.0 * (aFitHigh - aFitLow) / mNFitBins;

  mContent = (double *) malloc(sizeof(double) * mNFitBins);
  mErr2    = (double *) malloc(sizeof(double) * mNFitBins);

  mWeightSumP = new TH2D("WeightSumP", "WeightSumP", mNFitBins/2, aFitLow, aFitHigh, 100, 0.0, 2.0);
  mWeightSumN = new TH2D("WeightSumN", "WeightSumN", mNFitBins/2, aFitLow, aFitHigh, 100, 0.0, 2.0);
}

CalcCFNonId& CalcCFNonId::operator=(CalcCFNonId& aCalcCFNonId) 
{
  mPairManager = aCalcCFNonId.mPairManager;

  mNormP = aCalcCFNonId.mNormP;
  mNormN = aCalcCFNonId.mNormN;
  
  mFitLowKStar = aCalcCFNonId.mFitLowKStar;
  mFitHighKStar = aCalcCFNonId.mFitHighKStar;
  mNFitBins = aCalcCFNonId.mNFitBins;
  mBinWidth = aCalcCFNonId.mBinWidth;
  
  mNormLowKStar = aCalcCFNonId.mNormLowKStar;
  mNormHighKStar = aCalcCFNonId.mNormHighKStar;
}

short 
CalcCFNonId::InFitRange(Pair *aPair)
{
  return ((aPair->GetKStar() > mFitLowKStar) && (aPair->GetKStar() < mFitHighKStar));
}

short 
CalcCFNonId::InNormRange(Pair *aPair)
{
  return ((aPair->GetKStar() > mNormLowKStar) && (aPair->GetKStar() < mNormHighKStar));
}

int 
CalcCFNonId::FitRangeBin(Pair *aPair)
{
//   if (aPair->GetKStarTimeOutSign() > 0.0)
//     return (mNFitBins/2) + floor(((aPair->GetKStar()-mFitLowKStar) * (mNFitBins/2)) / (mFitHighKStar - mFitLowKStar));
//   else
//     return (mNFitBins/2) - floor(((aPair->GetKStar()-mFitLowKStar) * (mNFitBins/2)) / (mFitHighKStar - mFitLowKStar)) - 1;
  double tWhere = (aPair->GetKStar() - mFitLowKStar) / mBinWidth;
  return (aPair->GetKStarTimeOutSign() > 0.0) ? (mNFitBins/2) + trunc(tWhere) : (mNFitBins/2) - trunc(tWhere) - 1;
}

int 
CalcCFNonId::NormRangeBin(Pair *aPair)
{
  return 0;
  //  return (mNNormBins/2) + (aPair->GetKStarTimeOutSign() > 0.0) ? floor((aPair->GetKStar()-mNormLowKStar) * (mNNormBins/2) / (mNormHighKStar - mNormLowKStar)) : -floor((aPair->GetKStar()-mNormLowKStar) * mNNormBins / (mNormHighKStar - mNormLowKStar)) - 1;
}

void CalcCFNonId::Write()
{
}

//  void Read(TFile *aFile, TKey *aKey);
void CalcCFNonId::Normalize()
{
}

void CalcCFNonId::WriteParameters()
{
}

TGraphErrors *CalcCFNonId::GetGraph()
{
  STR tGName;
  
  Double_t *tX = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);
  Double_t tBinWidth = (mFitHighKStar - mFitLowKStar) / (mNFitBins / 2);
  for (int ti=0; ti<mNFitBins/2; ti++) {
    tX[mNFitBins/2 + ti] = mFitLowKStar + ti * tBinWidth + tBinWidth / 2;
    tX[mNFitBins/2 - 1 - ti] = -tX[mNFitBins/2 + ti];
  }
  TGraphErrors *tg = new TGraphErrors(mNFitBins, tX, mContent, 0, mErr2);
  tGName = "CalcGraph";

  tg->SetName((tGName+GetPairSystemName(mPairSystem)).Data());
  return tg;
}

TGraphErrors *CalcCFNonId::GetGraph(double aPurity)
{
  STR tGName;
  
  Double_t *tX    = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);
  Double_t *tCPur = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);
  Double_t *tEPur = (Double_t *) malloc(sizeof(Double_t) * mNFitBins);

  Double_t tBinWidth = (mFitHighKStar - mFitLowKStar) / (mNFitBins / 2);
  for (int ti=0; ti<mNFitBins/2; ti++) {
    tX[mNFitBins/2 + ti] = mFitLowKStar + ti * tBinWidth + tBinWidth / 2;
    tX[mNFitBins/2 - 1 - ti] = -tX[mNFitBins/2 + ti];
  }
  for (int ti=0; ti<mNFitBins; ti++) {
    tCPur[ti] = (mContent[ti] - 1.0)*aPurity + 1.0;
    tEPur[ti] = mErr2[ti]*aPurity;
  }
  tGName = "CalcGraph";

  TGraphErrors *tg = new TGraphErrors(mNFitBins, tX, tCPur, 0, tEPur);
  tg->SetName((tGName+GetPairSystemName(mPairSystem)).Data());
  return tg;
}

void CalcCFNonId::Calculate()
{
  //  if(mAlreadyCalcCF) return;

  double tMiu = 1.0;
  double tWsmiusum[mNFitBins];
   
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
    tWsmiusum[ti]=0.;
  }
  mNormN = 0.0;
  mNormP = 0.0;

  Pair *tPair;
  double tWeight;

  mWeightSumP->Reset("ICE");
  mWeightSumN->Reset("ICE");

  for(int ti=0; ti<mPairManager->GetFitPairCount(); ti++){
    tPair = mPairManager->GetFitPair(ti);
    tPair->SetPosition();
    tWeight = mWeightCalc->getWeight(*tPair);
    mContent[tPair->GetBin()] += tWeight;
    mErr2[tPair->GetBin()] += (tWeight*tWeight);
    if (tPair->GetBin() > mNFitBins/2)
      mWeightSumP->Fill(tPair->GetKStar(), tWeight);
    else
      mWeightSumN->Fill(tPair->GetKStar(), tWeight);
      
    tWsmiusum[tPair->GetBin()] += (tWeight - tMiu); 
  }
  for(int ti=0; ti<mPairManager->GetNormPairCount(); ti++){
    tPair = mPairManager->GetNormPair(ti);
    tPair->SetPosition();
    if (tPair->GetKStarTimeOutSign() > 0.0) {
      mNormP+=(mWeightCalc->getWeight(*tPair));
    }
    else {
      mNormN+=(mWeightCalc->getWeight(*tPair));
    }
  }

  mNormN /= mNNormNPair;
  mNormP /= mNNormPPair;

  double tMean;
  int tj;
  double tSigma;
  
  for(int ti=0; ti<mNFitBins/2; ti++){
    tMean = mContent[ti]/mNPairInFitBin[ti];
    //    mErr2[ti] = (mErr2[ti]/mNPairInFitBin[ti]-tMean*tMean)/mNPairInFitBin[ti];
    tSigma = mErr2[ti] + mNPairInFitBin[ti] * (1.0 + (1.0 - tMean)*(1.0 - tMean))
      - 2 * mContent[ti] + 2 * (1.0 - tMean) * tWsmiusum[ti];
    //    tSigma = TMath::Sqrt(tSigma);
    mErr2[ti] = tSigma / (mNPairInFitBin[ti] * mNPairInFitBin[ti]);

    if(mErr2[ti]<0){
      //PRINT_DEBUG( mErr2[ti] << " " << tMean;
    }
    //PRINT_DEBUG( std::endl;
    mContent[ti] = tMean/mNormN;

    tj = ti+mNFitBins/2;
    tMean = mContent[tj]/mNPairInFitBin[tj];
    mErr2[tj] = (mErr2[tj]/mNPairInFitBin[tj]-tMean*tMean)/mNPairInFitBin[tj];
    if(mErr2[tj]<0){
      PRINT_DEBUG("Cannot calculate error " << mErr2[tj] << " " 
	   << mNPairInFitBin[tj] << " " << tMean );
    }
    mContent[tj] = tMean/mNormP;
  }

  //  mAlreadyCalcCF=1;
  //  mCFStorage->store(this);

    
}

void CalcCFNonId::SetNormPBins(Int_t aBins) 
{
  mNNormPPair = aBins;
}

void CalcCFNonId::SetNormNBins(Int_t aBins)
{
  mNNormNPair = aBins;
}

int  
CalcCFNonId::ReadFromStorage()
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
    PRINT_MESSAGE("" << GetCFIdentifier() << " [Read from storage]");
  }
  catch (CF_Storage_Exception ce)
    {
      PRINT_MESSAGE("" << GetCFIdentifier() << " [Calculating]");
      return 0;
    }
  if ((tC) && (tE)) {
    if (tC->GetNoElements() == mNFitBins) {
      for (int ti=0; ti<tC->GetNoElements(); ti++) {
	mContent[ti] = (*tC)(ti);
	mErr2[ti]    = (*tE)(ti);
	PRINT_DEBUG_3("Bin " << ti << " " << mContent[ti] << "+/-" << mErr2[ti]);
      }
    }
    else {
      PRINT_MESSAGE("Error in stored function " << GetCFIdentifier());
      return 0;
    }
    
  }
  else {
    PRINT_MESSAGE("Error in stored function " << GetCFIdentifier());
    return 0;
  }
  
  
  return 1;
} 
 
void 
CalcCFNonId::WriteToStorage()
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

  for (int ti=0; ti<tC->GetNoElements(); ti++) {
    PRINT_DEBUG_3("Bin " << ti << " " << mContent[ti] << "+/-" << mErr2[ti]);
  }

  delete tC;
  delete tE;
}

const char *
CalcCFNonId::GetCFIdentifier()
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

void 
CalcCFNonId::WriteSums()
{
  mWeightSumN->Write();
  mWeightSumP->Write();
  
}

