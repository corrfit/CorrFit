#ifndef _CORRFIT_SOURCEMODELGAUSROUTFILE_H_
#define _CORRFIT_SOURCEMODELGAUSROUTFILE_H_

#include <TRandom.h>
#include "CFGlobal.h"
#include "SourceModel.h"

class SourceModelGausROutFile: public SourceModel
{
 public:
  SourceModelGausROutFile();
  ~SourceModelGausROutFile();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aPair, VectorPair *aPos, double *aRandVar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual const char* GetModelIdentifier();

 protected:
  double mParameters[2];
  static TRandom *mRandom;
};

#endif
