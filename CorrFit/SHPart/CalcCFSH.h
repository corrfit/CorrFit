#ifndef _CORRFIT_CALCCF_SH_H_
#define _CORRFIT_CALCCF_SH_H_

#include "TKey.h"
#include "TFile.h"
#include "CalcCF.h"
#include "TH2D.h"
#include <TRandom2.h>
#include "../SphericalHarmonics/CorrFctnDirectYlm.h"

//#define MAXFITPAIRS 200000
//#define MAXNORMPAIRS 40000

class CalcCFSH: public CalcCF
{
 public:
  CalcCFSH();
  CalcCFSH(Double_t aFitLow, Double_t aFitHigh, Int_t aFitBin,
	   Double_t aNormLow, Double_t aNormHigh);
  ~CalcCFSH();
  CalcCFSH& operator=(CalcCFSH& aCalcCFSH);

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

  void SetPrimaryFraction(Double_t aFraction, Int_t aSeed=0);

/*   void SetNormPBins(Int_t aBins); */
/*   void SetNormNBins(Int_t aBins); */

/*   void WriteSums(); */
  
 protected:
  virtual int  ReadFromStorage();
  virtual void WriteToStorage();
  virtual void Calculate();
  const char  *GetCFIdentifier();
  
  CorrFctnDirectYlm *mCFYlm;

/*   Double_t mNormP; */
/*   Double_t mNormN; */

/*   Int_t mNNormNPair; */
/*   Int_t mNNormPPair; */
  
  Double_t mFitLowKStar;
  Double_t mFitHighKStar;
  
  Double_t mNormLowKStar;
  Double_t mNormHighKStar;

  Double_t mBinWidth;

  Double_t mPrimaryFraction;

  TRandom2 *mRand;

/*   TH2D *mWeightSumP; */
/*   TH2D *mWeightSumN; */
  
};

#endif
