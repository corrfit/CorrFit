#ifndef _CORRFIT_CALCCF_1DHBT_H_
#define _CORRFIT_CALCCF_1DHBT_H_

#include "TKey.h"
#include "TFile.h"
#include "CalcCF.h"

#define MAXFITPAIRS 200000
#define MAXNORMPAIRS 40000

class CalcCF1DHBT: public CalcCF
{
 public:
  CalcCF1DHBT();
  CalcCF1DHBT(Double_t aFitLow, Double_t aFitHigh, Int_t aFitBin,
	      Double_t aNormLow, Double_t aNormHigh);
  ~CalcCF1DHBT();
  CalcCF1DHBT& operator=(CalcCF1DHBT& aCalcCF1DHBT);

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

  void SetNormBins(Int_t aBins);

 protected:
  virtual int  ReadFromStorage();
  virtual void WriteToStorage();
  virtual void Calculate();
  const char  *GetCFIdentifier();
  
  Double_t mNorm;

  Int_t mNNormPair;
  
  Double_t mFitLowKStar;
  Double_t mFitHighKStar;
  
  Double_t mNormLowKStar;
  Double_t mNormHighKStar;
};

#endif
