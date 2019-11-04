#ifndef _CORRFIT_CALCCF_QUAD_H_
#define _CORRFIT_CALCCF_QUAD_H_

#include "TKey.h"
#include "TFile.h"
#include "CalcCF.h"
#include "TH2D.h"

//#define MAXFITPAIRS 200000
//#define MAXNORMPAIRS 40000

class CalcCFQuad: public CalcCF
{
 public:
  CalcCFQuad();
  CalcCFQuad(Double_t aFitLow, Double_t aFitHigh, Int_t aFitBin,
	      Double_t aNormLow, Double_t aNormHigh);
  ~CalcCFQuad();
  CalcCFQuad& operator=(CalcCFQuad& aCalcCFQuad);

  virtual void Write();
  //  virtual void Read(TFile *aFile, TKey *aKey);
  virtual void Normalize();
  virtual void WriteParameters();

  virtual TGraphErrors *GetGraph();
  virtual TGraphErrors *GetGraph(double aPurity);

  virtual short InFitRange(Pair *aPair);
  virtual short InNormRange(Pair *aPair);
  virtual int FitRangeBin(Pair *aPair);
  virtual int NormRangeBin(Pair *aPair);

  void SetNormPBins(Int_t aBins);
  void SetNormNBins(Int_t aBins);

  void WriteSums();
  
 protected:
  virtual int  ReadFromStorage();
  virtual void WriteToStorage();
  virtual void Calculate();
  const char  *GetCFIdentifier();
  
  Double_t mNormP;
  Double_t mNormN;

  Int_t mNNormNPair;
  Int_t mNNormPPair;
  
  Double_t mFitLowKStar;
  Double_t mFitHighKStar;
  
  Double_t mNormLowKStar;
  Double_t mNormHighKStar;

  Double_t mBinWidth;

  TH2D *mWeightSumP;
  TH2D *mWeightSumN;
  
};

#endif
