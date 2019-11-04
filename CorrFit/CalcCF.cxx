#include "fstream"
#include "CalcCF.h"
#include "ReadPar.h"

WeightCalculator *CalcCF::mWeightCalc = new StandAloneFsiLednicky;
#ifdef MYSQLSTORAGE
CFStorageMySQL             *CalcCF::mStorage    = NULL;
#else
CFStorage             *CalcCF::mStorage    = NULL;
#endif
extern ReadPar *sRPInstance;
extern STR      sRPFileName;

using namespace std;

CalcCF::CalcCF()
{
  mPairManager = 0;
  mPairType = PiPlusPiPlus;
  mNFitBins = 0;
  mNPairInFitBin = NULL;
  ReadParameters();
  InitWeight();
  InitStorage();
}

CalcCF::CalcCF(PairManager *aPairManager, Int_t aNFitBins)
{
  mPairManager = aPairManager;
  mPairType = PiPlusPiPlus;
  mNFitBins = aNFitBins;

  mNPairInFitBin = (Int_t *) malloc(sizeof(Int_t) * mNFitBins);
  for (int ti=0; ti<mNFitBins; ti++)
    mNPairInFitBin[ti] = 0;
  ReadParameters();
  InitWeight();
  InitStorage();
}

CalcCF::~CalcCF()
{
  mPairManager = 0;
}

void CalcCF::SetPairManager(PairManager *aPairManager)
{
  mPairManager = aPairManager;
}

CalcCF& CalcCF::operator=(CalcCF& aCalcCF)
{
  mPairManager = aCalcCF.mPairManager;
  mPairType = aCalcCF.mPairType;
}

void                   
CalcCF::SetPairSystem(CorrFitPairSystem aPairSystem)
{
  PRINT_DEBUG("Setting Pair System to " << aPairSystem);
  mPairSystem = aPairSystem;
  InitWeight();
  InitStorage();
}

CorrFitPairSystem
CalcCF::GetPairSystem()
{
  return mPairSystem;
}

Int_t                  
CalcCF::GetNFitBins()
{
  return mNFitBins;
}

void                   
CalcCF::SetFitBinCount(Int_t aBin, Int_t count)
{
  mNPairInFitBin[aBin] = count;
}

void 
CalcCF::InitWeight()
{
  if (mWeightCalcName == "FsiLednicky")
    mWeightCalc = new StandAloneFsiLednicky;
  else if (mWeightCalcName == "FsiPratt") {
    StandAloneFSIPratt::ReadParameters();
    mWeightCalc = new StandAloneFSIPratt;
  }
  else if (mWeightCalcName == "FsiKisiel") {
    mWeightCalc = new StandAloneFsiKisiel();
  }
  else
    {
      PRINT_MESSAGE("Unknown weight calculator: " << mWeightCalcName.Data() << endl <<
		    "Select an existing one");
      exit(0);
    }
  
  if (mCoulomb)
    mWeightCalc->setCoulOn();
  else
    mWeightCalc->setCoulOff();
  if (mStrong)
    mWeightCalc->setStrongOn();
  else
    mWeightCalc->setStrongOff();
  if (mQuantum)
    mWeightCalc->setQuantumOn();
  else
    mWeightCalc->setQuantumOff();
  mWeightCalc->set3BodyOff();
  mWeightCalc->ReadParameters();

  //  mWeightCalc->setPairType(PiPlusKaonMinus);
  switch (mPairSystem) {
  case PiplusKplus:
  case PiminusKminus:
    mPairType = PiPlusKaonPlus;
    mWeightCalc->setPairType(PiPlusKaonPlus);
    break;
  case PiplusKminus:
  case PiminusKplus:
    mPairType = PiPlusKaonMinus;
    mWeightCalc->setPairType(PiPlusKaonMinus);
    break;
  case PiplusP:
  case PiminusAntiP:
    mPairType = PiPlusProton;
    mWeightCalc->setPairType(PiPlusProton);
    break;
  case PiplusAntiP:
  case PiminusP:
    mPairType = PiMinusProton;
    mWeightCalc->setPairType(PiMinusProton);
    break;
  case KplusP:
  case KminusAntiP:
    mPairType = KaonPlusProton;
    mWeightCalc->setPairType(KaonPlusProton);
    break;
  case KplusAntiP:
  case KminusP:
    mPairType = KaonMinusProton;
    mWeightCalc->setPairType(KaonMinusProton);
    break;
  case PiplusPiplus:
    mPairType = PiPlusPiPlus;
    mWeightCalc->setPairType(PiPlusPiPlus);
    break;
  case PiplusPiminus:
    mPairType = PiPlusPiMinus;
    mWeightCalc->setPairType(PiPlusPiMinus);
    break;
  case KplusKplus:
    mPairType = KaonPlusKaonPlus;
    mWeightCalc->setPairType(KaonPlusKaonPlus);
    break;
  case KplusKminus:
    mPairType = KaonPlusKaonMinus;
    mWeightCalc->setPairType(KaonPlusKaonMinus);
    break;
  case KzeroKzero:
    mPairType = KaonZeroKaonZero;
    mWeightCalc->setPairType(KaonZeroKaonZero);
    break;
  case PP:
    mPairType = ProtonProton;
    mWeightCalc->setPairType(ProtonProton);
    break;
  case PAntiP:
    mPairType = ProtonAntiProton;
    mWeightCalc->setPairType(ProtonAntiProton);
    break;
  case NN:
    mPairType = NeutronNeutron;
    mWeightCalc->setPairType(NeutronNeutron);
    break;
  case NP:
    mPairType = NeutronProton;
    mWeightCalc->setPairType(NeutronProton);
    break;
  case PD:
    mPairType = ProtonDeuteron;
    mWeightCalc->setPairType(ProtonDeuteron);
    break;
  case ND:
    mPairType = NeutronDeuteron;
    mWeightCalc->setPairType(NeutronDeuteron);
    break;
  case PT:
    mPairType = ProtonTriton;
    mWeightCalc->setPairType(ProtonTriton);
    break;
  case PA:
    mPairType = ProtonAlfa;
    mWeightCalc->setPairType(ProtonAlfa);
    break;
  case DD:
    mPairType = DeuteronDeuteron;
    mWeightCalc->setPairType(DeuteronDeuteron);
    break;
  case DT:
    mPairType = DeuteronTriton;
    mWeightCalc->setPairType(DeuteronTriton);
    break;
  case DA:
    mPairType = DeuteronAlfa;
    mWeightCalc->setPairType(DeuteronAlfa);
    break;
  case TT:
    mPairType = TritonTriton;
    mWeightCalc->setPairType(TritonTriton);
    break;
  case TA:
    mPairType = TritonAlfa;
    mWeightCalc->setPairType(TritonAlfa);
    break;
  case AA:
    mPairType = AlfaAlfa;
    mWeightCalc->setPairType(AlfaAlfa);
    break;
  case KzeroSKzeroS:
    mPairType = KaonZeroAntiKaonZero;
    mWeightCalc->setPairType(KaonZeroAntiKaonZero);
    break;
  default:
    PRINT_MESSAGE("No such PairSystem: " << mPairSystem);
    exit(0);
    break;
  }

  PRINT_DEBUG("Theoretical calculations for " << mCoulomb << " " << mStrong << " " << mQuantum << " " << mPairType);

}

void           
CalcCF::GenerateParameterStub()
{
  STR     tWeightCalcName;
  fstream *os;
  os = new fstream(sRPFileName.Data(),ios::ate|ios::out|ios::in);

  (*os) << "Possible pair types are: " << endl;
  for (int iter = 0; iter <= LastSystem; iter++)
    if (strlen(GetPairSystemName((CorrFitPairSystem) iter)) > 1) {
      (*os) << "# " << GetPairSystemName((CorrFitPairSystem) iter) << " = " << iter << endl;
    }
  (*os) << endl;
    
  (*os) 
    << "# Theoretical function calculator parameters" << endl 
    << "# Which interactions should be included in the calculation" << endl
    << "StrongInteraction = " << endl 
    << "CoulombInteraction = " << endl 
    << "QuantumStatistics = " << endl << endl
    << "# The name of the file, where the calculated correlation functions will be stored" << endl
    << "OutCalculationFileName = " << endl
    << endl;
  os->close();

  tWeightCalcName = sRPInstance->getPar("FsiCalculator");
  
  if (tWeightCalcName == "FsiLednicky")
    StandAloneFsiLednicky::GenerateParameterStub();
  else if (tWeightCalcName == "FsiPratt")
    StandAloneFSIPratt::GenerateParameterStub();
  else if (tWeightCalcName == "FsiKisiel")
    StandAloneFsiKisiel::GenerateParameterStub();
}

void           
CalcCF::ReadParameters()
{
  STR fVar;
  
  try 
    {
      fVar = sRPInstance->getPar("StrongInteraction");
      if ((fVar == "1") || (fVar == "On") || (fVar == "on"))
	mStrong = 1;
      else
	mStrong = 0;
      fVar = sRPInstance->getPar("CoulombInteraction");
      if ((fVar == "1") || (fVar == "On") || (fVar == "on"))
	mCoulomb = 1;
      else
	mCoulomb = 0;
      fVar = sRPInstance->getPar("QuantumStatistics");
      if ((fVar == "1") || (fVar == "On") || (fVar == "on"))
	mQuantum = 1;
      else
	mQuantum = 0;
      mStorageFileName = sRPInstance->getPar("OutCalculationFileName");
      mWeightCalcName = sRPInstance->getPar("FsiCalculator");
    }
  
  catch (STR e)
    {
      PRINT_MESSAGE("Error reading parameters in CalcCF::ReadParameters: " << e);
    }

  try
    {
      mPairSystem = (CorrFitPairSystem) atoi((sRPInstance->getPar("PairType")).Data());
    }
  catch (STR e)
    {
      PRINT_MESSAGE("Default PairType not found!");
      mPairSystem = PiplusKplus;
    }
  
}

void 
CalcCF::Generate()
{
  if (!ReadFromStorage())
    {
      Calculate();
//       for (int ti=0; ti<mNFitBins; ti++) {
// 	PRINT_DEBUG("Bin " << ti << " : " << mContent[ti] << " +/- " << mErr2[ti]);
//       }
      
      WriteToStorage();
    }
}

void 
CalcCF::InitStorage()
{
  if (!mStorage) {
#ifdef MYSQLSTORAGE
    mStorage = new CFStorageMySQL();
#else
    mStorage = new CFStorage(mStorageFileName.Data());
#endif
  }
  mPairNumber = mStorage->SetPairName(GetPairSystemName(mPairSystem));
  mStorage->SetCalculationParameters(mPairNumber, mStrong, mCoulomb, mQuantum);
}

const char *
CalcCF::GetPairSystemName(CorrFitPairSystem aSystem)
{
  switch (aSystem) {
  case PiplusKplus    : return "PiPlusKPlus";
  case PiminusKminus  : return "PiMinusKMinus";
  case PiplusKminus   : return "PiPlusKMinus";
  case PiminusKplus   : return "PiMinusKPlus";
  case PiplusP        : return "PiPlusP";
  case PiminusAntiP   : return "PiMinusAntiP";
  case PiplusAntiP    : return "PiPlusAntiP";
  case PiminusP       : return "PiMinusP";
  case KplusP         : return "KPlusP";
  case KminusAntiP    : return "KMinusAntiP";
  case KplusAntiP     : return "KPlusAntiP";
  case KminusP        : return "KMinusP";
  case PiplusPiplus   : return "PiPlusPiPlus";
  case PiminusPiminus : return "PiMinusPiMinus";
  case PiplusPiminus  : return "PiPlusPiMinus";
  case KplusKplus     : return "KPlusKPlus";
  case KminusKminus   : return "KMinusKMinus";
  case KplusKminus    : return "KPlusKMinus";
  case KzeroKzero     : return "KZeroKZero";
  case PP             : return "ProtonProton";
  case PAntiP         : return "ProtonAntiProton";
  case NN             : return "NeutronNeutron";
  case NP             : return "NeutronProton";
  case PD             : return "ProtonDeuteron";
  case ND             : return "NeutronDeuteron";
  case PT             : return "ProtonTriton";
  case PA             : return "ProtonAlfa";
  case DD             : return "DeuteronDeuteron";
  case DT             : return "DeuteronTriton";
  case DA             : return "DeuteronAlfa";
  case TT             : return "TritonTriton";
  case AA             : return "AlfaAlfa";
  case KzeroSKzeroS   : return "KZeroShortKZeroShort";
  }
  return "";
}

