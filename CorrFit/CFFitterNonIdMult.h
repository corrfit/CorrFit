#ifndef _CORRFIT_CFFITTERNONIDMULT_H_
#define _CORRFIT_CFFITTERNONIDMULT_H_

#include "CFGlobal.h"
#include "CFFitter.h"
#include "CalcCFNonId.h"
#include "ExpCFNonId.h"

class CFFitterNonIdMult: public CFFitter
{
 public:
  CFFitterNonIdMult();
  ~CFFitterNonIdMult();
  
  virtual void   Initialize();
  virtual void   Fit();
  virtual void   Write();
  virtual void   ReadParameters();

 protected:
  virtual void   InitializePairManager();
  virtual double GetChi2();
  virtual double GetChi2(double aPurity);
  virtual void   GenerateParameterStub();
  virtual void   GenerateCF();
  virtual void   GeneratePairFileNameStub(fstream *os);
  virtual void   GeneratePairs(int *aPairsPerBin, int aNBins);
  
  ExpCFNonId  **mExpCFs;
  CalcCFNonId **mCalcCFs;
  
  int  mPType;
  int  mPTypes[4];
  int  mNPTypes;
  int  mRandomSeed;
  int  mMinPairsPerBin;
  STR  mPFileNames[4];
  PairManager *mPairManagers[4];
  
};

#endif
