#ifndef _CORRFIT_SOURCEMODELGAUSRINV_H_
#define _CORRFIT_SOURCEMODELGAUSRINV_H_

#include <TRandom.h>
#include <math.h>
#include <TMath.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelGausRInv: public SourceModel
{
 public:
  SourceModelGausRInv();
  ~SourceModelGausRInv();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetModelIdentifier();

 protected:
  double mParameters[5];
  double mNegProbability;
  static TRandom *mRandom;
};

#endif
