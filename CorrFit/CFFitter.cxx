#include "fstream"
#include "CFFitter.h"

extern ReadPar *sRPInstance;
extern STR      sRPFileName;

using namespace std;

SourceModel *CFFitter::mSourceModel = NULL;
int CFFitter::maxpairperbin = 20000;

CFFitter::CFFitter()
{
  //  ReadParameters();
}

CFFitter::~CFFitter()
{
  if (mChi2Map) delete mChi2Map;
  if (mChi2DetMap) delete mChi2DetMap;
}

void     
CFFitter::InitializePairManager()
{
  Int_t  tPCount;
  Pair  *tPair;
  Int_t *mNFPBuf;
  Int_t  tPairUsed;
  
  mPairManager = new PairManager(mPairFileName.Data());
  mCalcCF->SetPairManager(mPairManager);

  mNFPBuf = (Int_t *) malloc(sizeof(Int_t) * mCalcCF->GetNFitBins());
  for (int ti=0; ti<mCalcCF->GetNFitBins(); ti++)
    mNFPBuf[ti] = 0;
  
  for (tPCount=0; tPCount<mPairManager->GetPairCount(); tPCount++){
    tPair = mPairManager->ReadPair();
    tPairUsed=0;
    
    if (mCalcCF->InFitRange(tPair)) {
      tPair->SetBin(mCalcCF->FitRangeBin(tPair));
      if (mNFPBuf[mCalcCF->FitRangeBin(tPair)] < maxpairperbin) {
	mNFPBuf[mCalcCF->FitRangeBin(tPair)]++;
	mPairManager->StoreFitPair(tPair);
	tPairUsed = 1;
	//      PRINT_DEBUG("Pair " << tPCount << " in fit range, Bin " << (mCalcCF->FitRangeBin(tPair)));
      }
    }
    if (mCalcCF->InNormRange(tPair)) {
      if (mPairManager->GetNormPairCount() < maxpairperbin*3) {
	tPair->SetBin(mCalcCF->NormRangeBin(tPair));
	mPairManager->StoreNormPair(tPair);
	tPairUsed = 1;
      }
      //      PRINT_DEBUG("Pair " << tPCount << " in norm range");
    }
    //    D_("Read pair " << tPair->GetKStar());
    
    if (!tPairUsed) delete tPair;
  }
  
  mPairManager->CloseFile();

  for (int ti=0; ti<mCalcCF->GetNFitBins(); ti++) {
    mCalcCF->SetFitBinCount(ti, mNFPBuf[ti]);
    if (!mNFPBuf[ti]) {
      PRINT_MESSAGE("Bin number " << ti << " has no pairs. Aborting");
      //      exit(1);
    }
    else {
      PRINT_DEBUG("Bin number " << ti << " has " << mNFPBuf[ti] << " pairs");
    }
  }
  
  if (!(mPairManager->GetNormPairCount())) {
    PRINT_MESSAGE("Normalization bin has no pairs. Aborting");
    exit(1);
  }
  
  if (!(mCalcCF->GetNFitBins())) {
    PRINT_MESSAGE("There are zero fitting bins. Check the input file and the debug information. Aborting");
    exit(1);
  }

  D_("CFFitter::Initialize");
  D_("Number of Fit Pairs read: " << mPairManager->GetFitPairCount());
  D_("Number of Norm Pairs read: " << mPairManager->GetNormPairCount());


  GeneratePairs(mNFPBuf, mCalcCF->GetNFitBins());

  for (tPCount=0; tPCount<mPairManager->GetFitPairCount(); tPCount++)
    {
      tPair = mPairManager->GetFitPair(tPCount);
      tPair->InitRandVar();
    }
  
  for (tPCount=0; tPCount<mPairManager->GetNormPairCount(); tPCount++)
    {
      tPair = mPairManager->GetNormPair(tPCount);
      tPair->InitRandVar();
    }
}


void 
CFFitter::Initialize()
{
  ReadParameters();
  InitializePairManager();
}

void 
CFFitter::Fit()
{
  int     tDim;
  int*    tParNBins;
  double* tParCoarseMin;
  double* tParMins;
  double* tParMaxs;
  double  tBestCoarsePur;
  double  tBestCoarseChi;
  TString tOutput;
  
  tDim = mSourceModel->GetNParameters();
  
  tParCoarseMin = (double *) malloc(sizeof(double) * tDim);
  tParMins      = (double *) malloc(sizeof(double) * tDim);
  tParMaxs      = (double *) malloc(sizeof(double) * tDim);
  tParNBins     = (int *)    malloc(sizeof(int) * tDim);

  // Generate a coarse Chi2 Map

  mChi2Map = GenerateMap(tDim, mParMins, mParMaxs, mParBins, &tParCoarseMin, &tBestCoarsePur, &tBestCoarseChi);
  
  tOutput = "\nMinimum from the coarse map: [";
  for (int ti=0; ti<tDim; ti++) {
    if (ti) tOutput += ", ";
    tOutput += Form("%.3lf", tParCoarseMin[ti]);
  }
  tOutput += "] for Purity: ";
  tOutput += Form("%.2lf", tBestCoarsePur);
  PRINT_MESSAGE(tOutput);

  // Generate a detailed Chi2 Map
  // around the mean of the coarse Map

  for (int ti=0; ti<tDim; ti++) {
    tParNBins[ti] = mDetSize[ti];
    tParMins[ti] = tParCoarseMin[ti] - ((mDetSize[ti] - 1.0) / 2.0) * mDetStep[ti];
    tParMaxs[ti] = tParCoarseMin[ti] + ((mDetSize[ti] - 1.0) / 2.0) * mDetStep[ti];
  }

  mChi2DetMap = GenerateMap(tDim, tParMins, tParMaxs, tParNBins, &mBestPar, &mBestPur, &mBestChi);

  tOutput = "\nFinal Minimum from the map: [";
  for (int ti=0; ti<tDim; ti++) {
    if (ti) tOutput += ", ";
    tOutput += Form("%.3lf", mBestPar[ti]);
  }
  tOutput += "] for Purity: ";
  tOutput += Form("%.2lf", mBestPur);
  tOutput += " and Chi2Min: ";
  tOutput += Form("%.2lf", mBestChi); 
  PRINT_MESSAGE(tOutput);
}

void 
CFFitter::Write()
{
  TString tOutput; 
  
  // Write ExpCF characteristics
  mExpCF->Write();
  
  // Write the best fit values
  int tDim = mSourceModel->GetNParameters();
  TVectorD *tResult = new TVectorD(tDim+1);
  tOutput = "\nFinal output parameters: [";
  for (int ti=0; ti<tDim; ti++) {
    if (ti) tOutput += ", ";
    tOutput += Form("%.5lf", mBestPar[ti]);
    (*tResult)(ti) = mBestPar[ti];
  }
  tOutput += "]";
  (*tResult)(tDim) = mBestPur;
  PRINT_MESSAGE(tOutput);
  tResult->Write();

  // Write additional histograms
  Pair *tpair;
  TH1D *pbeta = new TH1D("pairbetat","Pair \\beta_t",100,0.0,1.0);
  for (int iter=0; iter<mPairManager->GetFitPairCount(); iter++)
    {
      tpair = mPairManager->GetFitPair(iter);
      pbeta->Fill(tpair->GetBetat());
    }
  pbeta->Write();
  
}

void 
CFFitter::SetSourceModel(SourceModel *aSourceModel)
{
  mSourceModel = aSourceModel;
  Pair::SetSourceModel(aSourceModel);
  PRINT_DEBUG("New SourceModel is " << Pair::GetSourceModel());
}

// Procedure to generate a chi2 map, using a 
// CalcCF generated simulated functions
// and ExpCF experimental function
// The dimensions and the size od the map
// is arbitrary
Chi2Map *
CFFitter::GenerateMap(int aDim, double *aParMins, double *aParMaxs, int* aParNBins, double **aBestPars, double *aBestPur, double *aBestChi)
{
  double  *tPar;
  double  *tValues;
  Chi2Map *tChi2Map;
  double  tPurStep;
  int     tCellCount = 0;
  double  tBestNorm;

  PRINT_DEBUG("Generating Chi2 map which has " << aDim << " dimension(s): ");
  
  // Create the helper tables
  tPar          = (Double_t *) malloc(sizeof(double) * aDim);
  tValues       = (Double_t *) malloc(sizeof(double) * mPurityBin);
  *aBestPars    = (double *) malloc(sizeof(double) * aDim);

  tPurStep      = (mPurityBin > 1) ? (mPurityMax - mPurityMin) / (mPurityBin - 1.0) : 0.0;

  tChi2Map = new Chi2Map(aDim, aParNBins, aParMins, aParMaxs);
  tChi2Map->SetPurityParams(mPurityBin, mPurityMin, mPurityMax);
  
  try {
    tChi2Map->InitFill();
  }
  catch(CorrFit_Chi2Map_Exception e) {
    PRINT_DEBUG("CFFitter::GenerateMap - Caught " << e << " when initializing map");
    exit(0);
  }
  try {
    while (1) {
      tChi2Map->GetNextCellParams(tPar);
      mSourceModel->SetModelParameters(tPar);
      GenerateCF();
      for (int tP=0; tP<mPurityBin; tP++) {
	tValues[tP] = GetChi2(mPurityMin+tPurStep*tP);
      }
      tChi2Map->SetCurrentCellContent(tValues);
      PRINT_DEBUG_2("Chi2Map Cell No. " << ++tCellCount << " per " << tChi2Map->GetMapSize());
    }
  }
  catch (CorrFit_Chi2Map_Exception e) {
    if (e!=CHI2MAP_EXCEPTION_END_OF_MAP) {
      PRINT_DEBUG("CFFitter::GenerateMap - Caught " << e << " when calculating the map");
      exit(0);
    }
  }
  try {
    tChi2Map->GetBestParams(*aBestPars, aBestPur, &tBestNorm, aBestChi);
  }
  catch(CorrFit_Chi2Map_Exception e) {
    PRINT_DEBUG("CFFitter::GenerateMap - Caught " << e << " when initializing map");
    exit(0);
  }
  free(tPar);
  free(tValues);

  return tChi2Map;
}

void    
CFFitter::GenerateCF()
{
  mCalcCF->Generate();
}

void     
CFFitter::GeneratePairFileNameStub(fstream *os)
{
  (*os) << endl
    << "# This is the name of the file containing input pairs" << endl
    << "InPairCalcName = " << endl << endl;
}

void 
CFFitter::GeneratePairs(int *aPairsPerBin, int aNBins)
{
  PRINT_DEBUG("Normal CFFitter - GeneratePairs is a placeholder");
}

void 
CFFitter::GenerateParameterStub()
{
  int tDim = 0;
  fstream *os;
  if (mSourceModel)
    tDim = mSourceModel->GetNParameters();
  
  PRINT_MESSAGE("Parameter file not complete. Some parameters are missing.");
  PRINT_MESSAGE("Updating a stub parameter file \"" << sRPFileName.Data() << "\"");
  PRINT_MESSAGE("Please edit it and fill all the fileds.");
  os = new fstream(sRPFileName.Data(),ios::ate|ios::out|ios::in);

  (*os) 
    << endl
    << "# Source model parameters" << endl 
    << "# Below is the list of model parameter names and values" << endl;
  for (int ti=0; ti<tDim; ti++) {
    STR tSNum = "";
    tSNum += (ti+1);
    
    (*os) << "# Parameter no " << (ti+1) << ": " << mSourceModel->GetParameterName(ti) << endl;
    (*os)
      << ("Parameter"+tSNum+"Min").Data() << " = " << endl
      << ("Parameter"+tSNum+"Max").Data() << " = " << endl
      << ("Parameter"+tSNum+"Bin").Data() << " = " << endl;
  }
  
  (*os) << endl
	<< "# This are the parameters of the datailed chi2 map" << endl;
  
  for (int ti=0; ti<tDim; ti++) {
    (*os) << "DetailedMapNBins"  << (ti+1) << " = " << endl
	  << "DetailedMapStep" << (ti+1) << " = " << endl << endl;
  }
  
  (*os) 
    << "# and the optional purity parameters" << endl
    << "PurityMin = " << endl
    << "PurityMax = " << endl
    << "PurityBins = " << endl << endl;
  (*os) 
    << "# The maximum number of pairs in each fit bin" << endl
    << "MaxPairsPerFitBin = 20000" << endl << endl;
  
  GeneratePairFileNameStub(os);
  os->close();
  
}

void 
CFFitter::ReadParameters()
{
  STR pName, tDetName;
  int tDim = 0;
  double tMomRes;
  int tPairType;
  
  if (mSourceModel) 
    tDim = mSourceModel->GetNParameters();

  mParMins = (double *) malloc(sizeof(double) * tDim);
  mParMaxs = (double *) malloc(sizeof(double) * tDim);
  mParBins = (int *)    malloc(sizeof(int) * tDim);
  mDetSize = (int *)    malloc(sizeof(double) * tDim);
  mDetStep = (double *) malloc(sizeof(double) * tDim);
  
  pName = "Parameter";
  tDetName = "DetailedMap";
  
  try {
    for (int ti=0; ti<tDim; ti++) {
      STR tSNum = "";
      tSNum += ti+1;
      
      D_("Reading " << (pName+tSNum+"Min").Data());
      
      mParMins[ti] = atof((sRPInstance->getPar((pName+tSNum+"Min").Data())).Data());
      mParMaxs[ti] = atof((sRPInstance->getPar((pName+tSNum+"Max").Data())).Data());
      mParBins[ti] = atoi((sRPInstance->getPar((pName+tSNum+"Bin").Data())).Data());
      mDetSize[ti] = atoi((sRPInstance->getPar((tDetName+"NBins"+tSNum).Data())).Data());
      mDetStep[ti] = atof((sRPInstance->getPar((tDetName+"Step"+tSNum).Data())).Data());
    }
    mPairFileName = sRPInstance->getPar("InPairCalcName");

    tMomRes  = atof((sRPInstance->getPar("MomentumResolutionCorrection")).Data())/100.0;
    tPairType = atoi((sRPInstance->getPar("PairType")).Data());
    
    if (tMomRes > 0.0001)
      Pair::SetMomRes(tMomRes, tPairType);

    if (mPairFileName == "")
      {
	PRINT_MESSAGE("The name of the file with pair momentum distribution not defined.");
	PRINT_MESSAGE("Please fill in the values in \"" << sRPFileName.Data() << "\" file.");
	PRINT_MESSAGE("Aborting");
	exit(0);
      }
    maxpairperbin = atoi((sRPInstance->getPar("MaxPairsPerFitBin")).Data());
    if (!maxpairperbin) { maxpairperbin = 20000; }
  }
  catch (STR e)
    {
      PRINT_MESSAGE("Error reading parameters in CFFitters::ReadParameters: " << e);
      PRINT_MESSAGE("Generating stub parameter file");
      GenerateParameterStub();
    }
  try {
    mPurityMin = atof((sRPInstance->getPar("PurityMin")).Data());
    mPurityMax = atof((sRPInstance->getPar("PurityMax")).Data());
    mPurityBin = atoi((sRPInstance->getPar("PurityBins")).Data());
  }
  catch (STR e)
    {
      D_("No purity parameters. Using one purity bin = 1.0");
      mPurityMin = 1.0;
      mPurityMax = 1.0;
      mPurityBin = 1;
    }
}
