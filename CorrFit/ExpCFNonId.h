#ifndef _CORRFIT_EXPCFNONID_H_
#define _CORRFIT_EXPCFNONID_H_

#include "CFGlobal.h"
#include "ExpCF.h"
#include "ReadPar.h"
#include "TH1D.h"
#include "TFile.h"

#define EXPCF_SYSTEM_COUNT 9
#define EXPCF_PART_COUNT 5
#define EXPCF_PURITY_PART_COUNT 5

enum _ExpCF_System {
  EXPCF_SYSTEM_PP = 0,
  EXPCF_SYSTEM_PM = 1,
  EXPCF_SYSTEM_MP = 2,
  EXPCF_SYSTEM_MM = 3,
  EXPCF_SYSTEM_ZZ = 4,
  EXPCF_SYSTEM_ZP = 5,
  EXPCF_SYSTEM_ZM = 6,
  EXPCF_SYSTEM_PZ = 7,
  EXPCF_SYSTEM_MZ = 8,
};

typedef _ExpCF_System ExpCF_System;

enum _ExpCF_Type {
  EXPCF_TYPE_NUM = 0,
  EXPCF_TYPE_DEN = 1
};

typedef _ExpCF_Type ExpCF_Type;

enum _ExpCF_Proj {
  EXPCF_PROJ_OUT = 0,
  EXPCF_PROJ_SIDE = 1,
  EXPCF_PROJ_LONG = 2
};

typedef _ExpCF_Proj ExpCF_Proj;

enum _ExpCF_Sign {
  EXPCF_SIGN_POS = 0,
  EXPCF_SIGN_NEG = 1
};

typedef _ExpCF_Sign ExpCF_Sign;

enum _ExpCF_Part {
  EXPCF_PART_TYPE = 1,
  EXPCF_PART_PROJ = 2,
  EXPCF_PART_SIGN = 3,
  EXPCF_PART_SYSTEM = 4,
  EXPCF_PART_GENERAL = 5
};

typedef _ExpCF_Part ExpCF_Part;

enum _ExpCF_Purity_Part {
  EXPCF_PURITY_PART_PREFIX = 1,
  EXPCF_PURITY_PART_SYSTEM = 2,
  EXPCF_PURITY_PART_PROJ = 3,
  EXPCF_PURITY_PART_SIGN = 4,
  EXPCF_PURITY_PART_SUFFIX = 5
};

typedef _ExpCF_Purity_Part ExpCF_Purity_Part;

enum _ExpCF_Purity_Type {
  EXPCF_PURITY_TYPE_SINGLE = 0,
  EXPCF_PURITY_TYPE_DOUBLE = 1
};

typedef _ExpCF_Purity_Type ExpCF_Purity_Type;

class ExpCFNonId: public ExpCF
{
 public:
  ExpCFNonId();
  ExpCFNonId(int aChargeType,
	     TFile* aOutFile);
  virtual ~ExpCFNonId();

  virtual TGraphErrors* GetGraph(const char* aName);
  virtual TGraphErrors* GetGraph(const char* aName, double aPur);

  virtual double GetContent(int aIndex, double aPurity) const;
  virtual double GetContent(int aIndex, double aPurity, double aNorm) const;
   
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
  int mCFReduceEdgeEffects;
  
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

  friend class CFFitterNonId;
  friend class CFFitterNonIdMult;
  friend class CFFitterNonIdOS;
};

#endif
