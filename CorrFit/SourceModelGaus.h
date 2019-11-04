#ifndef _CORRFIT_SOURCEMODELGAUS_H_
#define _CORRFIT_SOURCEMODELGAUS_H_

#include <TRandom.h>
#include "CFGlobal.h"
#include "SourceModel.h"

enum _SourceModelGausType {
  SOURCEMODELGAUS_ONEGAUS_ONERADIUS,
  SOURCEMODELGAUS_ONEGAUS_THREERADII,
  SOURCEMODELGAUS_TWOGAUS_ONERADIUS,
  SOURCEMODELGAUS_TWOGAUS_THREERADII,
};

typedef enum _SourceModelGausType SourceModelGausType;

class SourceModelGaus: public SourceModel
{
 public:
  SourceModelGaus();
  SourceModelGaus(SourceModelGausType aType);
  ~SourceModelGaus();
  
  virtual void        SetModelParameters(double *aParTable);
  virtual double     *GetModelParameters();
  virtual void        GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar);
  virtual const char* GetParameterName(int aPar);
  virtual void        InitRandVar(double** aRandVar);
  virtual const char* GetModelIdentifier();

 protected:
  double mParameters[8];
  SourceModelGausType mType;
  static TRandom *mRandom;
};

#endif
