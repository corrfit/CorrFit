#ifndef _CORRFIT_CALCCF_H_
#define _CORRFIT_CALCCF_H_

#include "TGraphErrors.h"
#include "WeightCalculator.h"
#include "StandAloneFsiLednicky.h"
#include "StandAloneFSIPratt.h"
#include "StandAloneFsiKisiel.h"
#include "WeightCalculator.h"
#include "CF.h"
#include "PairManager.h"
#ifdef MYSQLSTORAGE
#include "CFStorageMySQL.h"
#else
#include "CFStorage.h"
#endif
#include "PairSystems.h"

class CalcCF: public CF 
{
 public:
  CalcCF();
  CalcCF(PairManager *aPairManager, Int_t aNFitBins);
  ~CalcCF();
  CalcCF& operator=(CalcCF& aCalcCF);

  void SetPairManager(PairManager *aPairManager);
  virtual void Generate();

  virtual TGraphErrors  *GetGraph() = 0;
  virtual TGraphErrors  *GetGraph(double aPurity) = 0;
  virtual short          InFitRange(Pair *aPair) = 0;
  virtual short          InNormRange(Pair *aPair) = 0; 
  virtual int            FitRangeBin(Pair *aPair) = 0;
  virtual int            NormRangeBin(Pair *aPair) = 0;
  virtual void           WriteParameters() = 0;
  void                   SetPairSystem(CorrFitPairSystem aPairSystem);
  CorrFitPairSystem      GetPairSystem();
  void                   SetFitBinCount(Int_t aBin, Int_t count);
  Int_t                  GetNFitBins();
  virtual void           ReadParameters();
  static const char     *GetPairSystemName(CorrFitPairSystem aSystem);
  static void            GenerateParameterStub();
  
 protected:
  PairManager      *mPairManager;
  CorrFitPairType   mPairType;
  CorrFitPairSystem mPairSystem;
  int               mNFitBins;
  int              *mNPairInFitBin;
  STR               mStorageFileName;
  STR               mWeightCalcName;
  
  int               mStrong;
  int               mCoulomb;
  int               mQuantum;
  int               mPairNumber;

  virtual int  ReadFromStorage() = 0;
  virtual void WriteToStorage() = 0;
  virtual void Calculate() = 0;
  virtual void InitWeight();
  virtual void InitStorage();
  
  static WeightCalculator *mWeightCalc;
#ifdef MYSQLSTORAGE
  static CFStorageMySQL        *mStorage;
#else
  static CFStorage             *mStorage;
#endif
};

#endif
