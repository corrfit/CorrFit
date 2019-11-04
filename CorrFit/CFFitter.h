#ifndef _CORRFIT_CFFITTER_H_
#define _CORRFIT_CFFITTER_H_

#include "CFGlobal.h"
#include "CalcCF.h"
#include "ExpCF.h"
#include "PairManager.h"
#include "ReadPar.h"
#include "Chi2Map.h"

class CFFitter
{
 public:
  CFFitter();
  ~CFFitter();
  
  virtual void Initialize();
  virtual void Fit();
  virtual void Write();
  static  void SetSourceModel(SourceModel *aSourceModel);
  virtual void ReadParameters();
 protected:
  CalcCF      *mCalcCF;
  ExpCF       *mExpCF;
  PairManager *mPairManager;
  static SourceModel *mSourceModel;
  Chi2Map     *mChi2Map;
  Chi2Map     *mChi2DetMap;
  double      *mBestPar;
  double       mBestPur;
  double       mBestChi;
  int         *mDetSize;
  double      *mDetStep;
  double      *mParMins;
  double      *mParMaxs;
  int         *mParBins;
  STR          mPairFileName;
  double       mPurityMin;
  double       mPurityMax;
  int          mPurityBin;

  static int maxpairperbin;

  virtual void     InitializePairManager();
  virtual double   GetChi2() = 0;
  virtual double   GetChi2(double aPurity) = 0;
  virtual Chi2Map *GenerateMap(int aDim, double *aParMins, double *aParMaxs, int* aParNBins, double **aBestPars, double *aBestPur, double *aBestChi);
  virtual void     GenerateParameterStub();
  virtual void     GenerateCF();
  virtual void     GeneratePairFileNameStub(fstream *os);
  virtual void     GeneratePairs(int *aPairsPerBin, int aNBins);
};

#endif
