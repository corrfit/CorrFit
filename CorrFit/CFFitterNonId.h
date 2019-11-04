#ifndef _CORRFIT_CFFITTERNONID_H_
#define _CORRFIT_CFFITTERNONID_H_

#include "CFGlobal.h"
#include "CFFitter.h"
#include "CalcCFNonId.h"
#include "ExpCFNonId.h"

class CFFitterNonId: public CFFitter
{
 public:
  CFFitterNonId();
  ~CFFitterNonId();
  
  //  virtual void   Initialize();
  virtual void   Fit();
  virtual void   Write();
  virtual void   ReadParameters();

 protected:
  virtual double GetChi2();
  virtual double GetChi2(double aPurity);
  virtual void   GenerateParameterStub();
  virtual void   GeneratePairs(int *aPairsPerBin, int aNBins);
  virtual void   InitializePairManager();
  
  int mPType;
  int mMinPairsPerBin;
  int mRandomSeed;
  int mGeneratePairs;

  TH2D *mKStarKOutFrac;
  TH2D *mKStarKSideFrac;
  TH2D *mKStarKLongFrac;
};



#endif
