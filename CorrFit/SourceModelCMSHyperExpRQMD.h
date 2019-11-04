#ifndef _CORRFIT_SOURCEMODELCMSHYPERRQMDEXP_H_
#define _CORRFIT_SOURCEMODELCMSHYPERRQMDEXP_H_

#include <TRandom.h>
#include <math.h>
#include <TMath.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelCMSHyperExpRQMD: public SourceModel
{
 public:
  SourceModelCMSHyperExpRQMD();
  ~SourceModelCMSHyperExpRQMD();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aPair, VectorPair *aPos, double *aRandVar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual const char* GetModelIdentifier();

  static TRandom *mRandom;
 protected:
  double GetHyper(Double_t aMu, Double_t aAlfa, Double_t aSigma, Double_t aHalfWide);
  double GetHyperDouble(Double_t aMu1, Double_t aAlfa1, Double_t aSigma1, Double_t aMu2, Double_t aAlfa2, Double_t aSigma2, Double_t aHalfWide, Double_t aSwitch);
  double GetHyperExp(Double_t aMu, Double_t aAlfa, Double_t aSigma, Double_t aHalfWide, Double_t aExpScale, Double_t aExpSlope);
  double GetHyperDoubleExp(Double_t aMu1, Double_t aAlfa1, Double_t aSigma1, Double_t aMu2, Double_t aAlfa2, Double_t aSigma2, Double_t aHalfWide, Double_t aSwitch, Double_t aExpScale, Double_t aExpSlope);
  double ScaleHyper(Double_t aScale, Double_t aX);
  
  double mParameters[5];
  double mWideProbability;
  double mHalfWides[4];
};

#endif
