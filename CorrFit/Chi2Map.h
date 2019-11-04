#ifndef _CORRFIT_CHI2MAP_
#define _CORRFIT_CHI2MAP_

#include "CFGlobal.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#ifdef RANDOMC2MORDER
#include <list>
#endif

enum _CorrFit_Chi2Map_Exception {
  CHI2MAP_EXCEPTION_NO_MAP,
  CHI2MAP_EXCEPTION_END_OF_MAP,
  CHI2MAP_EXCEPTION_PARAMETERS_NOT_SET
};

typedef enum _CorrFit_Chi2Map_Exception CorrFit_Chi2Map_Exception;

// This class stores the Chi2 values in a map
// of n dimensions, where n is dependent on the
// source model used.
// It can save the map in a format of table of doubles
// or root histograms

class Chi2Map
{
 public:
  Chi2Map();
  Chi2Map(Int_t aDimension, Int_t *aNBins, Double_t *aBinMins, Double_t *aBinMaxs);
  ~Chi2Map();
  
  // Initialization functions
  void      SetPurityParams(Int_t aPurityBins, Double_t aPurityMin, Double_t aPurityMax);
  void      SetNormParams(Int_t aNormBins, Double_t aNormMin, Double_t aNormMax);
  void      SetMapParams(Int_t aDimension, Int_t *aNBins, Double_t *aBinMins, Double_t *aBinMaxs);

  // Functions used when filling and calculating the map
  void      InitFill() throw(CorrFit_Chi2Map_Exception);
  void      SetCellContent(Int_t *aCoords, Double_t *values) throw (CorrFit_Chi2Map_Exception);
  void      GetBestParams(Double_t *aParams, Double_t *aPurity, Double_t *aNorm, Double_t *aChiMin) throw(CorrFit_Chi2Map_Exception);
  void      GetNextCellParams(Double_t *aPars) throw(CorrFit_Chi2Map_Exception);
  void      SetCurrentCellContent(Double_t *values) throw(CorrFit_Chi2Map_Exception);
  void      MinuitMinimize() throw(CorrFit_Chi2Map_Exception);

  // Input/output functions
  void      WriteHisto();
  int       GetMapSize();
  
 private:
  Double_t *mMapStorage;
  Int_t    *mParBins;
  Double_t *mParMins;
  Double_t *mParMaxs;
  Double_t *mParSteps;
  Int_t     mPurityBins;
  Double_t  mPurityMin;
  Double_t  mPurityMax;
  Double_t  mPurityStep;
  Int_t     mNormBins;
  Double_t  mNormMin;
  Double_t  mNormMax;
  Double_t  mNormStep;
  Int_t    *mCurCoords;
  Int_t     mCurPurity;
  Int_t     mCurNorm;
  Int_t     mDimension;
  Int_t    *mCurCoordsInc;
  Double_t *mBestPars;
  Double_t  mBestPurity;
  Double_t  mBestNorm;
  Double_t  mBestChi;
  Double_t  mChiMin;
  Int_t     mMapSize;
  Int_t    *mSubMapSizes;

#ifdef RANDOMC2MORDER
  list<Int_t> mOrder;
  list<Int_t>::iterator mOrderCur;
#endif
};

#endif
