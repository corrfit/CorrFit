#ifndef _CORRFIT_SOURCEMODELGAUSROUTSOURCECMS_H_
#define _CORRFIT_SOURCEMODELGAUSROUTSOURCECMS_H_

#include <TRandom.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelGausROutSourceCMS: public SourceModel
{
 public:
  SourceModelGausROutSourceCMS();
  ~SourceModelGausROutSourceCMS();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aPair, VectorPair *aPos, double *aRandVar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual const char* GetModelIdentifier();

  static TRandom *mRandom;
 protected:
  double mParameters[5];
  double mWideProbability;
};

#endif
