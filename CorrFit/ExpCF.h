#ifndef _CORRFIT_EXPCF_H_
#define _CORRFIT_EXPCF_H_

#include "CFGlobal.h"
#include "CF.h"
#include "TGraphErrors.h"

class ExpCF: public CF
{
 public:
  ExpCF();
  virtual ~ExpCF();

  virtual TGraphErrors* GetGraph(const char* aName) = 0;
  virtual TGraphErrors* GetGraph(const char* aName, double aPur) = 0;

  virtual double GetContent(int aIndex, double aPurity) const = 0;
  virtual double *GetContent();
  virtual double *GetPurity();
   
  virtual void WriteHisto() = 0;
  int GetNFitBin();
  double GetMomResScale();
  
 protected:
  double mNormN;
  double mNormP; 

  double mPurCorr;
  double mMomRes;

  int mNFitBin;

  double *mPurity;

  virtual void ApplyPurCorr(int aChargeType) = 0;
  virtual void ApplyMomResCorr(double aMomResCorr) = 0;
};

#endif
