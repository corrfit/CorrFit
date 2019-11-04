#ifndef _CORRFIT_SOURCEMODELEXPGUSRINV_H_
#define _CORRFIT_SOURCEMODELEXPGAUSRINV_H_

#include <TRandom.h>
#include <math.h>
#include <TMath.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelExpGausRInv: public SourceModel
{
 public:
  SourceModelExpGausRInv();
  ~SourceModelExpGausRInv();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetModelIdentifier();

 protected:
  double GetXFromY(double aY, double aPar);
  double FunExp(double aArg, double aPar);
  double mParameters[3];
  double mExpProbability;
  static TRandom *mRandom;
};

#endif
