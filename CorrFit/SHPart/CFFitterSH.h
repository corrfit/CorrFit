#ifndef _CORRFIT_CFFITTERSH_H_
#define _CORRFIT_CFFITTERSH_H_

#include "CFGlobal.h"
#include "CFFitter.h"
#include "CalcCFSH.h"
#include "ExpCFSH.h"

class CFFitterSH: public CFFitter
{
 public:
  CFFitterSH();
  ~CFFitterSH();
  
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
  double mPrimaryFraction;
};



#endif
