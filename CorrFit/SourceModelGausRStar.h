#ifndef _CORRFIT_SOURCEMODELGAUSRSTAR_H_
#define _CORRFIT_SOURCEMODELGAUSRSTAR_H_

#include <TRandom.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelGausRStar: public SourceModel
{
 public:
  SourceModelGausRStar();
  ~SourceModelGausRStar();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aPair, VectorPair *aPos, double *aRandVar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual const char* GetModelIdentifier();

  static TRandom *mRandom;
 protected:
  double mParameters[1];
};

#endif
