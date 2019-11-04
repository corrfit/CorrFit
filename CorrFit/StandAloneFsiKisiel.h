/***************************************************************************
 *
 * $Id: 
 *
 * Author: Adam Kisiel
 ***************************************************************************
 *
 * Description : Calculate correlation weight using function symmetrization
 *               and exact Coulomb solution
 *
 ***************************************************************************
 *
 * $Log: 
 *
 ***************************************************************************/

#ifndef StandAloneFsiKisiel_hh
#define StandAloneFsiKisiel_hh

#include <math.h>
#include <TMath.h>
#include "Pair.h"
#include "PairSystems.h"
#include "WeightCalculator.h"
#include "ReadPar.h"

struct _dcomplex 
{
  long double re;
  long double im;
};

typedef struct _dcomplex dcomplex;

class StandAloneFsiKisiel: public WeightCalculator {
 public: 
  // --- Constructor
  StandAloneFsiKisiel(); // call setDefaultCalcPar and setPipPip 
  // --- Destructor : nothing to explicitly delete
  ~StandAloneFsiKisiel() {};

  // --- Function to be called in the correlation function
  virtual double getWeight(Pair& aPair);

  // >>> For other pairs, use the function below
  virtual void setPairType(const int aPairType); 

  // >>> Calculation mode
  virtual void setDefaultCalcPar(); // Default is CoulOn, QuantumOn, StrongOn, 3BodyOff
  virtual void setCoulOn();
  virtual void setCoulOff();
  virtual void setQuantumOn();
  virtual void setQuantumOff();
  virtual void setStrongOn();
  virtual void setStrongOff();
  virtual void set3BodyOn();
  virtual void set3BodyOff();
  static void ReadParameters();
  static void GenerateParameterStub();
  void PrintNext();

protected:
  void   InitializeGamow();
  long double Gamow(double arg);
  void PairKinematics(Pair &aPair);
  void GetFFdouble(dcomplex *ffp, dcomplex *ffm, Pair &aPair);
  void GetFFsingle(dcomplex *ffp, int sign, Pair &aPair);
  long double geth(long double eta);
  long double chiim(long double eta);
  void bfunpfun(long double eta, long double rho, long double &bret, long double &pret);
  dcomplex GetG(long double eta, long double rho, long double hfun);
  void Getfc(long double kstar, long double eta, long double hfun, dcomplex &fcs, dcomplex &fct);

  double GetCoulomb(Pair &aPair);
  double GetQuantumCoulomb(Pair &aPair);
  double GetQuantum(Pair &aPair);
  double GetQuantumCoulombStrong(Pair &aPair);

  long double mBohrRadius;
  long double mGamovMesh[4000];
  long double euler;

  dcomplex d0s;
  dcomplex f0s;
  dcomplex d0t;
  dcomplex f0t;

  int twospin;

  double mROS, mRSS, mRLS, mRTS, mRSt;

  int mPrint;
};

#endif
