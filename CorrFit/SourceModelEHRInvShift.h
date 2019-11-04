#ifndef _CORRFIT_SOURCEMODELEHRINVSHIFT_H_
#define _CORRFIT_SOURCEMODELEHRINVSHIFT_H_

#include <TRandom.h>
#include <math.h>
#include <TMath.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelEHRInvShift: public SourceModel
{
 public:
  SourceModelEHRInvShift();
  ~SourceModelEHRInvShift();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetModelIdentifier();

 protected:
  double GetHyper(Double_t aMu, Double_t aAlfa, Double_t aSigma);

  double mParameters[5];
  double mNegProbability;
  static TRandom *mRandom;
};

#endif
