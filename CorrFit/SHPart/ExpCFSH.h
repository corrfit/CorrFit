#ifndef _CORRFIT_EXPCFSH_H_
#define _CORRFIT_EXPCFSH_H_

#include "CFGlobal.h"
#include "ExpCF.h"
#include "ReadPar.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TFile.h"

#define EXPCFSH_SYSTEM_COUNT 9
#define EXPCFSH_PART_COUNT 4
#define EXPCFSH_PURITY_PART_COUNT 4

enum _ExpCFSH_System {
  EXPCFSH_SYSTEM_PP = 0,
  EXPCFSH_SYSTEM_PM = 1,
  EXPCFSH_SYSTEM_MP = 2,
  EXPCFSH_SYSTEM_MM = 3,
  EXPCFSH_SYSTEM_ZZ = 4,
  EXPCFSH_SYSTEM_ZP = 5,
  EXPCFSH_SYSTEM_ZM = 6,
  EXPCFSH_SYSTEM_PZ = 7,
  EXPCFSH_SYSTEM_MZ = 8,
};

typedef _ExpCFSH_System ExpCFSH_System;

enum _ExpCFSH_Part {
  EXPCFSH_PART_TYPE = 1,
  EXPCFSH_PART_PROJ = 2,
  EXPCFSH_PART_SYSTEM = 3,
  EXPCFSH_PART_GENERAL = 4
};

typedef _ExpCFSH_Part ExpCFSH_Part;

enum _ExpCFSH_Purity_Part {
  EXPCFSH_PURITY_PART_PREFIX = 1,
  EXPCFSH_PURITY_PART_SYSTEM = 2,
  EXPCFSH_PURITY_PART_PROJ = 3,
  EXPCFSH_PURITY_PART_SUFFIX = 4
};

typedef _ExpCFSH_Purity_Part ExpCFSH_Purity_Part;

enum _ExpCFSH_Purity_Type {
  EXPCFSH_PURITY_TYPE_SINGLE = 0,
  EXPCFSH_PURITY_TYPE_DOUBLE = 1
};

typedef _ExpCFSH_Purity_Type ExpCFSH_Purity_Type;

class ExpCFSH: public ExpCF
{
 public:
  ExpCFSH();
  ExpCFSH(int aChargeType,
	     TFile* aOutFile);
  virtual ~ExpCFSH();

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
  
  virtual double* GetCovariance();
  virtual const double* GetCovariance() const;

 protected:
  virtual void ApplyPurCorr(int aChargeType);
  virtual void ApplyMomResCorr(double aMomResCorr);
  void SetParameters(int aNFitBin, double aLowKStarFit, double aHighKStarFit,
		     double aLowKStarNorm, double aHighKStarNorm);

/*   TH1D* mHHCF; // unsigned correlation function  */
/*   TH1D* mHCF;   */

  TH1D *mRC00;
  TH1D *mRC11;
  TH3D *mCov;

  double *mCovariance;

/*   TF1*  mHDenOutPFit; */
/*   TF1*  mHDenOutNFit;   */

  double mLowKStarFit;
  double mHighKStarFit; 
  double mLowKStarNorm;
  double mHighKStarNorm; 
  //  double* mKStarSigned;

  // Parts of the normal histogram name
  int mOrder[EXPCFSH_PART_COUNT];
  STR mTypes[2];
  STR mProjs[3];
  STR mSigns[3];
  STR mSystems[EXPCFSH_SYSTEM_COUNT];
  STR mGeneral;
  
  // Parts of the purity histogram name
  ExpCFSH_Purity_Type mPType;
  int mPOrder[EXPCFSH_PURITY_PART_COUNT];
  STR mPPrefix;
  STR mPProjs[3];
  STR mPSigns[3];
  STR mPSystems[EXPCFSH_SYSTEM_COUNT];
  STR mPSuffix;

  int mDataKStar;
  int mPurityKStar;
  int mBinominalErrors;
  //  int mCFReduceEdgeEffects;
  
  int          ReadCF(int aChargeType, TFile* aOutFile);
  int          ReadCFFromCorrectedFile(int aChargeType, TFile* aOutFile);	     
  void         Build(int aChargeType);
  virtual void Normalize();
  void         FillFitArray();
  STR          MakeHistName(int aType, int aProj, int aSystem);
  STR          MakePurityHistName(int aProj, int aSign, int aSystem);
  STR          MakeFileName(int aChargeType);
  STR          MakePurityFileName(int aChargeType);
  void         ReadParameters();

  friend class CFFitterSH;
};

#endif
