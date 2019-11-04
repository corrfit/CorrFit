#ifndef _CORRFIT_SOURCEMODELGAUSROUT_H_
#define _CORRFIT_SOURCEMODELGAUSROUT_H_

#include <TRandom.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelGausROut: public SourceModel
{
 public:
  SourceModelGausROut();
  ~SourceModelGausROut();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aPair, VectorPair *aPos, double *aRandVar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual const char* GetModelIdentifier();

  static TRandom *mRandom;
 protected:
  double mParameters[2];
};

#endif
