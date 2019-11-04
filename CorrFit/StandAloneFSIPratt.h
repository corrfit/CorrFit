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

#ifndef StandAloneFSIPratt_hh
#define StandAloneFSIPratt_hh

#include "Pair.h"
#include "PairSystems.h"
#include "WeightCalculator.h"
#include "wf/common/common.h"
#include "wf/wf_assorted.h"
#include "ReadPar.h"

static double  gPrattDelQ;
static int     gPrattNQmax;
static double  gPrattEpsilon;
static int     gPrattEllMax;
static int     gPrattStrong;
static int     gPrattCoulomb;
static STR     gPrattQArrayFilename;

class StandAloneFSIPratt: public WeightCalculator {
 public: 
  // --- Constructor
  StandAloneFSIPratt(); // call setDefaultCalcPar and setPipPip 
  // --- Destructor : nothing to explicitly delete
  ~StandAloneFSIPratt() {};

  // --- Function to be called in the correlation function
  virtual double getWeight(Pair& aPair);

  // --- Setting
  // >>> Set the type pair  
  virtual void setPipPip(); 
  virtual void setPimPim();
  virtual void setPipPim();
  virtual void setKpKp();
  virtual void setKmKm();
  virtual void setKpKm();
  virtual void setProtonProton();
  virtual void setPipProton();
  virtual void setPimProton();
// >>> For other pairs, use the function below
  virtual void setPairType(const int aPairType); 
// >>> The pair type can be chosen according to the following numbering
  // aPairType: 1  2  3  4    5   6   7   8  9 10  11  12  13  14 15 16 17
  //   part. 1: n  p  n  alfa pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
  //   part. 2: n  p  p  alfa pi- pi0 pi+ d  d K-  K+  p   p   K- K+ p  p
  // aPairType: 18 19   20  21   22  23  24 25 26
  //   part. 1: d  d    t   t    K0  K0  d  p  p
  //   part. 2: d  alfa t   alfa K0  K0b t  t  alfa
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

protected:
  // Setting parameters
  Twavefunction *mWaveFunction;
  
};

#endif
