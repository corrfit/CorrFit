#ifndef _CORRFIT_EXPCF1DHBT_H_
#define _CORRFIT_EXPCF1DHBT_H_

#include "CFGlobal.h"
#include "ExpCF.h"
#include "ReadPar.h"
#include "TH1D.h"
#include "TFile.h"
#include <TF1.h>
#include "ExpCFNonId.h"

class ExpCF1DHBT: public ExpCF
{
 public:
  ExpCF1DHBT();
  ExpCF1DHBT(int aChargeType,
	     TFile* aOutFile);
  virtual ~ExpCF1DHBT();
  
  virtual TGraphErrors* GetGraph(const char* aName);
  virtual TGraphErrors* GetGraph(const char* aName, double aPur);
  
  virtual double        GetContent(int aIndex, double aPurity) const;
   
  virtual void          WriteHisto();
  virtual void          Write();
  double                GetLowKStarFit();
  double                GetHighKSstarFit();
  double                GetLowKStarNorm();
  double                GetHighKStarNorm();
  static void           GenerateParameterStub();
  
 protected:
  virtual void ApplyPurCorr(int aChargeType);
  virtual void ApplyMomResCorr(double aMomResCorr);
  void         SetParameters(int aNFitBin, double aLowKStarFit, double aHighKStarFit,
			     double aLowKStarNorm, double aHighKStarNorm);

  TH1D* mHCF;  
  
  TF1*  mHDenFit;

  double mLowKStarFit;
  double mHighKStarFit; 
  double mLowKStarNorm;
  double mHighKStarNorm; 
  double* mKStarSigned;

  // Parts of the normal histogram name
  int mOrder[EXPCF_PART_COUNT];
  STR mTypes[2];
  STR mSystems[EXPCF_SYSTEM_COUNT];
  STR mGeneral;
  
  // Parts of the purity histogram name
  int mPOrder[EXPCF_PURITY_PART_COUNT];
  STR mPPrefix;
  STR mPSystems[EXPCF_SYSTEM_COUNT];
  STR mPSuffix;

  int mDataKStar;
  int mPurityKStar;
    
  int          ReadCF(int aChargeType, TFile* aOutFile);
  int          ReadCFFromCorrectedFile(int aChargeType, TFile* aOutFile);	     
  void         Build(int aChargeType);
  virtual void Normalize();
  void         FillFitArray();
  STR          MakeHistName(int aType, int aSystem);
  STR          MakePurityHistName(int aSystem);
  STR          MakeFileName(int aChargeType);
  STR          MakePurityFileName(int aChargeType);
  void         ReadParameters();

  friend class CFFitter1DHBT;

};

#endif
