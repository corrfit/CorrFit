/***************************************************************************
 *
 * $Id: 
 *
 * Author: Laurent Conin, Fabrice Retiere, Subatech, France
 ***************************************************************************
 *
 * Description : Calculate correlation weight using R.Lednicky's code 
 *  Use the fortran files : FsiLednickyWeight.F and FsiTools.F
 *
 ***************************************************************************
 *
 * $Log: 
 *
 ***************************************************************************/

#ifndef WeightCalculator_hh
#define WeightCalculator_hh

#include "Pair.h"

class WeightCalculator {
 public: 
  // --- Constructor
  WeightCalculator(); // call setDefaultCalcPar and setPipPip 
  // --- Destructor : nothing to explicitly delete
  ~WeightCalculator() {};

  // --- Function to be called in the correlation function
  virtual double getWeight(Pair& aPair) = 0;

  // --- Setting
  // >>> Set the type pair  
  // >>> For other pairs, use the function below
  virtual void setPairType(const int aPairType) = 0; 
  // >>> The pair type can be chosen according to the following numbering
  // aPairType: 1  2  3  4    5   6   7   8  9 10  11  12  13  14 15 16 17
  //   part. 1: n  p  n  alfa pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
  //   part. 2: n  p  p  alfa pi- pi0 pi+ d  d K-  K+  p   p   K- K+ p  p
  // aPairType: 18 19   20  21   22  23  24 25 26
  //   part. 1: d  d    t   t    K0  K0  d  p  p
  //   part. 2: d  alfa t   alfa K0  K0b t  t  alfa
  // >>> Calculation mode
  virtual void setDefaultCalcPar() = 0; // Default is CoulOn, QuantumOn, StrongOn, 3BodyOff
  virtual void setCoulOn() = 0;
  virtual void setCoulOff() = 0;
  virtual void setQuantumOn() = 0;
  virtual void setQuantumOff() = 0;
  virtual void setStrongOn() = 0;
  virtual void setStrongOff() = 0;
  virtual void set3BodyOn() = 0;
  virtual void set3BodyOff() = 0;
  virtual int  getPairType();
  static void ReadParameters() {}
  static void GenerateParameterStub() {}
  
protected:
  // Setting parameters
  int mLL;
  int mIch;
  int mIqs;
  int mIsi;
  int mI3c;

};

#endif
