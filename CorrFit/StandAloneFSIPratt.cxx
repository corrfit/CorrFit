#include "StandAloneFSIPratt.h"

#include "fstream"

extern ReadPar *sRPInstance;
extern STR      sRPFileName;

const double mPi = 139.0;
const double mK  = 493.0;
const double mP  = 938.0;
const double mLa = 1115.0;

StandAloneFSIPratt::StandAloneFSIPratt()
{
  mLL = 0;
  mWaveFunction = 0;
}

double StandAloneFSIPratt::getWeight(Pair& aPair)
{
  if (mWaveFunction) {
    aPair.SetPosMomEFirstNotation();
    return mWaveFunction->getpsisquared(aPair.mom2,aPair.pos2,aPair.mom1,aPair.pos1);
  }
  else{
    PRINT_MESSAGE("Trying to get weight without initializing the Wave Function!");
    return 0;
  }
}

void StandAloneFSIPratt::setPipPip() {
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mPi,mPi,1.0);
}
void StandAloneFSIPratt::setPimPim(){
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mPi,mPi,1.0);
}
void StandAloneFSIPratt::setPipPim(){
  mWaveFunction = new Twavefunction_pipluspiminus("parameters/wfparameters_pipluspiminus.dat");
}
void StandAloneFSIPratt::setKpKp(){
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mK,mK,1.0);
}
void StandAloneFSIPratt::setKmKm(){
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mK,mK,1.0);
}
void StandAloneFSIPratt::setKpKm(){
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_kpluskminus.dat",-1,mK,mK,0.5);
}
void StandAloneFSIPratt::setProtonProton(){
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mP,mP,1.0);
}
void StandAloneFSIPratt::setPipProton(){
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mPi,mP,0.5);
}
void StandAloneFSIPratt::setPimProton(){
  mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",-1,mPi,mP,0.5);
}

void StandAloneFSIPratt::setPairType(const int aPairType){ 
  mLL = aPairType;
  switch (aPairType) {
  case PiPlusPiPlus:
    mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mPi,mPi,1.0);
    break;
  case PiPlusPiMinus:
    mWaveFunction = new Twavefunction_pipluspiminus("parameters/wfparameters_pipluspiminus.dat");
    break;
  case KaonPlusKaonPlus:
    mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mK,mK,1.0);
    break;
  case KaonPlusKaonMinus:
    mWaveFunction = new Twavefunction_generic("parameters/wfparameters_kpluskminus.dat",-1,mK,mK,0.5);
    break;
  case ProtonProton:
    mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",1,mP,mP,1.0);
    break;
  case PiPlusProton:
    mWaveFunction = new Twavefunction_ppi("parameters/wfparameters_ppi.dat");
    break;
  case PiMinusProton:
    mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",-1,mPi,mP,0.5);
    break;
  case PiPlusKaonMinus:
    mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",-1,mPi,mK,0.5);
    break;
  case PiPlusKaonPlus:
    mWaveFunction = new Twavefunction_kpi("parameters/wfparameters_kpi.dat");
    break;
  case KaonPlusProton:
    mWaveFunction = new Twavefunction_pk("parameters/wfparameters_pk.dat");
    break;
  case KaonMinusProton:
    mWaveFunction = new Twavefunction_generic("parameters/wfparameters_generic.dat",-1,mK,mP,0.5);
    break;
  default:
    PRINT_MESSAGE("Scott Pratt's code does not support that pair type: " << aPairType);
    break;
  }
  
}

void StandAloneFSIPratt::setDefaultCalcPar(){ 
  mIch=1;
  mIqs=1;
  mIsi=1;
  mI3c=0;
}

void StandAloneFSIPratt::setCoulOn(){ 
  mIch = 1;
  gPrattCoulomb = 1;
}

void StandAloneFSIPratt::setCoulOff(){ 
  mIch = 0;
  gPrattCoulomb = 0;
}

void StandAloneFSIPratt::setQuantumOn(){
  mIqs = 1;
}

void StandAloneFSIPratt::setQuantumOff(){
  mIqs = 0;
}

void StandAloneFSIPratt::setStrongOn(){
  mIsi = 1;
  gPrattStrong = 1;
}

void StandAloneFSIPratt::setStrongOff(){
  mIsi = 0;
  gPrattStrong = 0;
}

void StandAloneFSIPratt::set3BodyOn(){
  mI3c = 1;
}

void StandAloneFSIPratt::set3BodyOff(){
  mI3c = 0;
}

void StandAloneFSIPratt::ReadParameters()
{
  STR pName;

  try {
    gPrattDelQ           = atof((sRPInstance->getPar("FsiPrattDelQ")).Data());
    gPrattNQmax          = atoi((sRPInstance->getPar("FsiPrattNQMax")).Data());
    gPrattEpsilon        = atof((sRPInstance->getPar("FsiPrattEpsilon")).Data());
    gPrattEllMax         = atoi((sRPInstance->getPar("FsiPrattEllMax")).Data());
    gPrattQArrayFilename = sRPInstance->getPar("FsiPrattQArrayFileName");
  }
  catch (STR e)
    {
      PRINT_MESSAGE("StandAloneFSIPratt: parameter " << e.Data() << " undefined");
      exit(0);
    }
  
}

void StandAloneFSIPratt::GenerateParameterStub()
{
  fstream *os;
  os = new fstream(sRPFileName.Data(),ios::ate|ios::out|ios::in);

  (*os) 
    << "# Parameters for the Scott Pratt's weight calculator" << endl 
    << "# The parameters of the Q mesh" << endl
    << "FsiPrattNQMax = 20" << endl 
    << "FsiPrattDelQ = 10.0" << endl 
    << "# If You put a filename here, a mesh points from this file will be taken" << endl
    << "# instead of evenly-spaced one" << endl
    << "FsiPrattQArrayFileName = " << endl
    << "# Other parameters (see instructions.txt for details)" <<  endl
    << "FsiPrattEpsilon = 1.0" << endl
    << "FsiPrattEllMax = 2" << endl
    << endl;
  os->close();
  
}


#include "wf/common/wf.cc"
#include "wf/common/misc.cc"
#include "wf/common/planewave.cc"
#include "wf/common/partwave.cc"
#include "wf/common/wfsubs.cc"
#include "wf/wf_generic.cc"
#include "wf/wf_kpi.cc"
#include "wf/wf_pipluspiminus.cc"
#include "wf/wf_pipluspiplus.cc"
#include "wf/wf_pk.cc"
#include "wf/wf_ppi.cc"
