#ifndef _CORRFIT_SOURCEMODELGAUSROUTTWO_H_
#define _CORRFIT_SOURCEMODELGAUSROUTTWO_H_

#include <TRandom.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelGausROutTwo: public SourceModel
{
 public:
  SourceModelGausROutTwo();
  ~SourceModelGausROutTwo();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aPair, VectorPair *aPos, double *aRandVar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual const char* GetModelIdentifier();

 protected:
  double mParameters[3];
  static TRandom *mRandom;
};

#endif
