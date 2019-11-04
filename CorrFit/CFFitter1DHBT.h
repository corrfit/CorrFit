#ifndef _CORRFIT_CFFITTER1DHBT_H_
#define _CORRFIT_CFFITTER1DHBT_H_

#include "CFGlobal.h"
#include "CFFitter.h"
#include "CalcCF1DHBT.h"
#include "ExpCF1DHBT.h"

class CFFitter1DHBT: public CFFitter
{
 public:
  CFFitter1DHBT();
  ~CFFitter1DHBT();
  
  virtual void Initialize();
  virtual void Fit();
  virtual void Write();
  virtual void ReadParameters();

 protected:
  virtual double GetChi2();
  virtual double GetChi2(double aPurity);
  virtual void GenerateParameterStub();
  
  int mPType;
  int mMinPairsPerBin;
  int mRandomSeed;
  int mGenerate;
  
  void GeneratePairs(int *aPairsPerBin, int aNBins);
  virtual void InitializePairManager();
  

};



#endif
