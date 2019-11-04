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

#include "StandAloneFsiLednicky.h"

StandAloneFsiLednicky::StandAloneFsiLednicky() {
  mLL=7; // Default PipPip
  mNs=4;
  mItest=0;
  mIch=1;
  mIqs=1;
  mIsi=1;
  mI3c=0;
  mNuclMass=1.;
  mNuclCharge=0.;
  mNuclChargeSign=1.0;
  FsiInit();
  FsiNucl();
};

StandAloneFsiLednicky::StandAloneFsiLednicky(
   int aPairType,bool aSphere,bool aNormMode,
   bool aCoul,bool aQuantum,bool aStrong,
   bool a3Body){
  mLL=aPairType;
  aSphere   ? mNs=4    : mNs=1;
  aNormMode ? mItest=0 : mItest=1;
  aCoul     ? mIch=1   : mIch=0;
  aQuantum  ? mIqs=1   : mIqs=0;
  aStrong   ? mIsi=1   : mIsi=0;
  a3Body    ? mI3c=1   : mI3c=0;
  mNuclMass=1.;
  mNuclCharge=0.;
  mNuclChargeSign=1.0;
  FsiInit();
  FsiNucl();
};

void StandAloneFsiLednicky::FsiInit(){
  fsiini(mItest,mLL,mNs,mIch,mIqs,mIsi,mI3c);
};

void StandAloneFsiLednicky::FsiNucl(){
  mNuclCharge*=mNuclChargeSign;
  fsinucl(mNuclMass,mNuclCharge);
};


