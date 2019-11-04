#ifndef _CORRFIT_EXPCFQUAD_H_
#define _CORRFIT_EXPCFQUAD_H_

#include "CFGlobal.h"
#include "ExpCF.h"
#include "ReadPar.h"
#include "TH1D.h"
#include "TFile.h"
#include "ExpCFNonId.h"

class ExpCFQuad: public ExpCF
{
 public:
  ExpCFQuad();
  ExpCFQuad(int aChargeType,
	     TFile* aOutFile);
  virtual ~ExpCFQuad();

  virtual TGraphErrors* GetGraph(const char* aName);
  virtual TGraphErrors* GetGraph(const char* aName, double aPur);
  
  virtual double GetContent(int aIndex, double aPurity) const;
  
  double GetPurity(int aIndex);
  virtual void WriteHisto();
  virtual void Write();
  double GetLowKStarFit();
  double GetHighKSstarFit();
  double GetLowKStarNorm();
  double GetHighKStarNorm();
  static void   GenerateParameterStub();
  
 protected:
  virtual void ApplyPurCorr(int aChargeType);
  virtual void ApplyMomResCorr(double aMomResCorr);
  void SetParameters(int aNFitBin, double aLowKStarFit, double aHighKStarFit,
		     double aLowKStarNorm, double aHighKStarNorm);

  TH1D* mHHCF; // unsigned correlation function 
  TH1D* mHCF;  
  double *mPurityContent;

  TF1*  mHDenOutPFit;
  TF1*  mHDenOutNFit;  

  double mLowKStarFit;
  double mHighKStarFit; 
  double mLowKStarNorm;
  double mHighKStarNorm; 
  double* mKStarSigned;

  // Parts of the normal histogram name
  int mOrder[EXPCF_PART_COUNT];
  STR mTypes[2];
  STR mProjs[3];
  STR mSigns[3];
  STR mSystems[EXPCF_SYSTEM_COUNT];
  STR mGeneral;
  
  // Parts of the purity histogram name
  ExpCF_Purity_Type mPType;
  int mPOrder[EXPCF_PURITY_PART_COUNT];
  STR mPPrefix;
  STR mPProjs[3];
  STR mPSigns[3];
  STR mPSystems[EXPCF_SYSTEM_COUNT];
  STR mPSuffix;

  int mDataKStar;
  int mPurityKStar;
  int mBinominalErrors;
  
  int          ReadCF(int aChargeType, TFile* aOutFile);
  int          ReadCFFromCorrectedFile(int aChargeType, TFile* aOutFile);	     
  void         Build(int aChargeType);
  virtual void Normalize();
  void         FillFitArray();
  STR          MakeHistName(int aType, int aProj, int aSign, int aSystem);
  STR          MakePurityHistName(int aProj, int aSign, int aSystem);
  STR          MakeFileName(int aChargeType);
  STR          MakePurityFileName(int aChargeType);
  void         ReadParameters();

  friend class CFFitterQuad;
};

#endif
