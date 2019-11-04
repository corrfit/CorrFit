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

#ifndef StandAloneFsiLednicky_hh
#define StandAloneFsiLednicky_hh

#include "WeightCalculator.h"
#include "Pair.h"

#include "fortranc.h"

// --- Prototype of the function used in the weight calculator 
//     (in FsiWeightLedinicky.F)
#define fsiini F77_NAME(fsiini,FSIINI)
extern "C" {void type_of_call F77_NAME(fsiini,FSIINI)(const int &itest,const int &ll,const int &ns,const int &ich, const int &iqs, const int &isi,const int &i3c);}
#define fsinucl F77_NAME(fsinucl,FSINUCL)
extern "C" {void type_of_call  F77_NAME(fsinucl,FSINUCL)(const double &mn,const double &cn);}
#define fsimomentum F77_NAME(fsimomentum,FSIMOMENTUM)
extern "C" {void type_of_call F77_NAME(fsimomentum,FSIMOMENTUM)(double &p1,double &p2);}
#define fsiposition F77_NAME(fsiposition,FSIPOSITION)
extern "C" {void type_of_call F77_NAME(fsiposition,FSIPOSITION)(double &x1,double &x2);}
#define fsiw F77_NAME(fsiw,FSIW)
extern "C" {void type_of_call F77_NAME(fsiw,FSIW)(const int &i,double &weif,
						  double &wei,double &wein);}
#define ltran12 F77_NAME(ltran12,LTRAN12)
extern "C" {void type_of_call ltran12_();}

// --- Additional prototyping of some CERN functions (in FsiTool.F)
typedef float   REAL;
typedef struct { REAL re; REAL im; } COMPLEX;
#define cgamma F77_NAME(cgamma,CGAMMA)
extern "C" {COMPLEX type_of_call cgamma_(COMPLEX*);}

class StandAloneFsiLednicky : public WeightCalculator {
 public: 
// --- Constructor
  StandAloneFsiLednicky(); // call setDefaultCalcPar and setPipPip 
  // avoid this constructor, kept for backward compatibility
  StandAloneFsiLednicky(int aPairType,
		   bool aSphere,
		   bool aNormMode,
		   bool aCoul,
		   bool aQuantum,
		   bool aStrong,
		   bool a3Body);
// --- Destructor : nothing to explicitly delete
  ~StandAloneFsiLednicky() {};

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
  void setSphere();
  void setSquare();
  virtual void setQuantumOn();
  virtual void setQuantumOff();
  virtual void setStrongOn();
  virtual void setStrongOff();
  virtual void set3BodyOn();
  virtual void set3BodyOff();
  void setNuclCharge(const double aNuclCharge);
  void setNuclMass(const double aNuclMass);
  virtual int  getPairType();
  static void ReadParameters() {}
  static void GenerateParameterStub() {}
  
private:
  // Fsi weight output
  double mWei;
  double mWein;
  double mWeif;
  double mWeightDen;
  // Setting parameters
  int mItest;
  int mNs;
  double mNuclMass;
  double mNuclCharge;
  short mNuclChargeSign;
  // Interface to the fortran functions
  void FsiInit();
  void FsiNucl();

};

inline  void StandAloneFsiLednicky::setPairType(const int aPairType) {mLL=aPairType;mNuclChargeSign=1;FsiInit();};
inline  void StandAloneFsiLednicky::setNuclCharge(const double aNuclCharge) {mNuclCharge=aNuclCharge;FsiNucl();};
inline  void StandAloneFsiLednicky::setNuclMass(const double aNuclMass){mNuclMass=aNuclMass;FsiNucl();};
inline  void StandAloneFsiLednicky::setSphere(){mNs=1;FsiInit();FsiNucl();};
inline  void StandAloneFsiLednicky::setSquare(){mNs=4;FsiInit();FsiNucl();};
inline  void StandAloneFsiLednicky::setDefaultCalcPar()
{mItest=1;mIch=1;mIqs=1;mIsi=1;mI3c=0;FsiInit();};
inline  void StandAloneFsiLednicky::setCoulOn(){mItest=1;mIch=1;FsiInit();};
inline  void StandAloneFsiLednicky::setCoulOff(){mItest=1;mIch=0;FsiInit();};
inline  void StandAloneFsiLednicky::setQuantumOn(){mItest=1;mIqs=1;FsiInit();};
inline  void StandAloneFsiLednicky::setQuantumOff(){mItest=1;mIqs=0;FsiInit();};
inline  void StandAloneFsiLednicky::setStrongOn(){mItest=1;mIsi=1;FsiInit();};
inline  void StandAloneFsiLednicky::setStrongOff(){mItest=1;mIsi=0;FsiInit();};
inline  void StandAloneFsiLednicky::set3BodyOn()
{mItest=1;mI3c=1;FsiInit();FsiNucl();};
inline  void StandAloneFsiLednicky::set3BodyOff()
{mItest=1;mI3c=0;FsiInit();mWeightDen=1.;FsiNucl();};

inline  void StandAloneFsiLednicky::setPipPip()      {mLL=7; mNuclChargeSign=1;FsiInit();FsiNucl();} 
inline  void StandAloneFsiLednicky::setPimPim()      {mLL=7; mNuclChargeSign=-1;FsiInit();FsiNucl();} 
inline  void StandAloneFsiLednicky::setPipPim()      {mLL=5; mNuclChargeSign=1;FsiInit();FsiNucl();} 
inline  void StandAloneFsiLednicky::setKpKp()        {mLL=15; mNuclChargeSign=1;FsiInit();FsiNucl();}
inline  void StandAloneFsiLednicky::setKmKm()        {mLL=15; mNuclChargeSign=-1;FsiInit();FsiNucl();}
inline  void StandAloneFsiLednicky::setKpKm()        {mLL=14; mNuclChargeSign=1;FsiInit();FsiNucl();}
inline  void StandAloneFsiLednicky::setProtonProton(){mLL=2; mNuclChargeSign=1;FsiInit();FsiNucl();}
inline  void StandAloneFsiLednicky::setPipProton()   {mLL=12; mNuclChargeSign=1;FsiInit();FsiNucl();}
inline  void StandAloneFsiLednicky::setPimProton()   {mLL=13; mNuclChargeSign=1;FsiInit();FsiNucl();}

inline double StandAloneFsiLednicky::getWeight(Pair& aPair){
  aPair.SetPosMom();
  fsimomentum(*Pair::mom1,*Pair::mom2);
  fsiposition(*Pair::pos1,*Pair::pos2);
  ltran12();
  fsiw(1,mWeif,mWei,mWein);
  //if (mI3c==0) return mWein;
  //mWeightDen=mWeif;
  return mWei;
}
inline int  StandAloneFsiLednicky::getPairType() { return mLL; }

#endif
