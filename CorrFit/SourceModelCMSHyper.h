#ifndef _CORRFIT_SOURCEMODELCMSHYPER_H_
#define _CORRFIT_SOURCEMODELCMSHYPER_H_

#include <TRandom.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelCMSHyper: public SourceModel
{
 public:
  SourceModelCMSHyper();
  ~SourceModelCMSHyper();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aPair, VectorPair *aPos, double *aRandVar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual const char* GetModelIdentifier();

  static TRandom *mRandom;
 protected:
  double GetHyper(Double_t aMu, Double_t aAlfa, Double_t aSigma);
  double GetHyperDouble(Double_t aMu1, Double_t aAlfa1, Double_t aSigma1, Double_t aMu2, Double_t aAlfa2, Double_t aSigma2);
  double ScaleHyper(Double_t aScale, Double_t aX);
  
  double mParameters[5];
  double mWideProbability;
};

#endif
