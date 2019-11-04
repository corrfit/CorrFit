#ifndef _CORRFIT_SOURCEMODEL_H_
#define _CORRFIT_SOURCEMODEL_H_

#include "CFGlobal.h"

class SourceModel 
{
 public:
  SourceModel();
  ~SourceModel();
  
  virtual void        SetModelParameters(double *aParTable) = 0;
  virtual double*     GetModelParameters() = 0;
  virtual void        GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar) = 0;
  virtual int         GetNParameters();
  virtual int         GetNRandVar();
  virtual void        InitRandVar(double** aRandVar) = 0;
  virtual const char* GetParameterName(int aPar) = 0;
  virtual const char* GetModelIdentifier() = 0;
  
 protected:
  int mNPar;
  int mNRandVar;
};

#endif
