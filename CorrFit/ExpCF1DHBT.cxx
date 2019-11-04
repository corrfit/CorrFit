#include "fstream"
#include "stdlib.h"
#include "ExpCF1DHBT.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TKey.h"
#include "TROOT.h"
#include "TError.h"
#include "ReadPar.h"

extern ReadPar *sRPInstance;
extern STR      sRPFileName;

using namespace std;

ExpCF1DHBT::ExpCF1DHBT() : ExpCF()
{
}

ExpCF1DHBT::ExpCF1DHBT(int aChargeType,
		       TFile* aOutFile) : ExpCF()
{
  ReadParameters();
  int tValid = ReadCF(aChargeType,aOutFile);
  if(tValid){
    int tiFirstFitBin = mHCF->FindBin(mLowKStarFit);
    int tiLastFitBin  = mHCF->FindBin(mHighKStarFit);
    int tiFirstScale = mHCF->FindBin(mLowKStarNorm);
    if(tiLastFitBin==tiFirstScale) tiFirstScale++;
    int tiLastScale = mHCF->FindBin(mHighKStarNorm);
    PRINT_DEBUG("Setting parameters:");
    PRINT_DEBUG("Fit bins. First: " << tiFirstFitBin << " Last: " << tiLastFitBin);
    PRINT_DEBUG("Norm bins. First: " << tiFirstScale << " Last: " << tiLastScale);
    SetParameters((tiLastFitBin-tiFirstFitBin+1),
		  mHCF->GetBinLowEdge(tiFirstFitBin),
		  mHCF->GetBinLowEdge(tiLastFitBin+1),
		  mHCF->GetBinLowEdge(tiFirstScale),
		  mHCF->GetBinLowEdge(tiLastScale+1));
    Build(aChargeType);
  }
}

ExpCF1DHBT::~ExpCF1DHBT()
{
}

TGraphErrors* ExpCF1DHBT::GetGraph(const char* aName)
{
}

TGraphErrors* ExpCF1DHBT::GetGraph(const char* aName, double aPur)
{
}

double ExpCF1DHBT::GetContent(int aIndex, double aPurity) const
{
  double tVar = mContent[aIndex];
  return (tVar-1.0) * aPurity +1.0;
}

void ExpCF1DHBT::WriteHisto()
{
}

void ExpCF1DHBT::Write()
{
  mHCF->Write();
}

void ExpCF1DHBT::ApplyPurCorr(int aChargeType)
{
  if(mPurCorr==0) return;

  TFile* tFIn = new TFile(MakePurityFileName(aChargeType).Data());
  if(!tFIn->IsOpen()) return;

  TH1D* tHCF = mHCF;

  TH2D* tHPur;

  // First check if the purity histogram is a profile
  // or a TH2D
  TKey *ktemp = tFIn->FindKey(MakePurityHistName(aChargeType).Data());
  if (ktemp)
    {
      TProfile *ttPPur;
      if (!strcmp(ktemp->GetClassName(),"TProfile"))
	{
	  // It's a TProfile!
	  // No need to do any fancy rewriting, just read it
	  
	  ttPPur = (TProfile *) tFIn->Get(MakePurityHistName(aChargeType).Data());
	  
	  if (!ttPPur) return;
	}
      else if (!strcmp(ktemp->GetClassName(),"TH2D"))
	{
	  // It's a TH2D 
	  // We need to rewrite it and cast it to a TProfile
	  TH2D *ttHPur;
	    
	  tHPur= (TH2D*) tFIn->Get(MakePurityHistName(aChargeType).Data());
	  PRINT_DEBUG("Try to get " << MakePurityHistName(aChargeType) << " " << tHPur);;

	  ttPPur = tHPur->ProfileX("ttPPur",1,ttHPur->GetNbinsY(),"E");
	}
      else 
	{
	  // It's neither a TProfile nor TH2D
	  // We bail out then
	  PRINT_MESSAGE("The key " << MakePurityHistName(aChargeType) << 
	     " is neither TProfile nor TH2D - cannot use" << endl <<
	     "Not applying Purity correction");
	}
      
      if(!mPurityKStar){
	PRINT_DEBUG("Rebinning KStar");
	int tNBin = ttPPur->GetNbinsX();
	double tXMin = ttPPur->GetXaxis()->GetXmin()/2.;
	double tXMax = ttPPur->GetXaxis()->GetXmax()/2.;
	ttPPur->GetXaxis()->Set(tNBin,tXMin,tXMax);
      }
      
      double tCorr;
      double tSignal;
      
      for(int ti=1;ti<=mHCF->GetNbinsX();ti++){
	tCorr = mPurCorr*ttPPur->GetBinContent(ti);
	if(tCorr>1.) tCorr=1.;
	tSignal = mHCF->GetBinContent(ti)-1.;
	mHCF->SetBinContent(ti,tSignal/tCorr+1);
	mHCF->SetBinError(ti,mHCF->GetBinError(ti)/tCorr);
	PRINT_DEBUG("Applying purity correction " << tCorr);
      }
      ttPPur->Delete();  
    }
  else
    {
      PRINT_MESSAGE("Purity histogram " << MakePurityHistName(aChargeType).Data() << " not found");
    }
  tFIn->Delete();
}

void ExpCF1DHBT::ApplyMomResCorr(double aMomResCorr)
{
}

int ExpCF1DHBT::ReadCF(int aChargeType, TFile* aOutFile)
{
  // A new code using the input parameter class

  TFile *tFIn=0;
  STR fName;
  
  PRINT_DEBUG("Reading the experimental CF");
  
  fName = MakeFileName(aChargeType);
  
  PRINT_DEBUG("Got FileName: " << fName);
  
  if(!tFIn){
    PRINT_DEBUG("Opening " << fName << " " << tFIn );
    tFIn = new TFile(fName.Data());
    PRINT_DEBUG("Opened " << tFIn);
  }
  if(!tFIn->IsOpen()) return 0;

  TH1D* tHCF=0;
  TH1D* tHDen=0;

  int ti=0;
  tHCF=0;

  char *temp;
  
  TString tDenname(MakeHistName(EXPCF_TYPE_DEN,aChargeType));
  TString tHCFname(MakeHistName(EXPCF_TYPE_NUM,aChargeType));
  
  tHDen = (TH1D*) tFIn->Get(tDenname);
  if (tHDen == 0x0 )
   {
     ::Fatal("","Can not find denominator histogram named %s",tDenname.Data());
     return 0;
   }
  tHCF  = (TH1D*) tFIn->Get(tHCFname);
  if (tHDen == 0x0 )
   {
     ::Fatal("","Can not find numerator histogram named %s",tHCFname.Data());
     return 0;
   }
  
  Double_t tFitMin = tHDen->GetBinLowEdge(tHDen->FindBin(mLowKStarFit));
  Double_t tFitMax = tHDen->GetBinLowEdge(tHDen->FindBin(mHighKStarFit)+1);
  TF1 *poll = new TF1("poll","[0]*x*x*x+[1]*x*x+[2]*x+[3]");

  tHDen->Fit("poll","QNO","",tFitMin,tFitMax);
  mHDenFit = new TF1(*poll);

  PRINT_DEBUG("Try to get " << MakeHistName(EXPCF_TYPE_NUM,aChargeType) << " " << tHCF);
    
  if(!tHCF){ return 0;}
  
  aOutFile->cd();
  
  tHCF->Divide(tHDen);

  STR tName;
  tName = "HCF";
  if (mSystems[aChargeType] != "")
    tName += mSystems[aChargeType];
  else
    switch (aChargeType)
      {
      case EXPCF_SYSTEM_PP:
	tName += "_PP_" ;
	break;
      case EXPCF_SYSTEM_MM:
	tName += "_MM_" ;
	break;
      }
  tName += '\0';

  tHCF->SetTitle(tName.Data());
  tHCF->SetName(tName.Data());
  
  mHCF = new TH1D(*tHCF);
  
  tFIn->Delete();

  if(!mDataKStar){
    int tNBin = mHCF->GetNbinsX();
    double tXMin = mHCF->GetXaxis()->GetXmin()/2.;
    double tXMax = mHCF->GetXaxis()->GetXmax()/2.;
    mHCF->GetXaxis()->Set(tNBin,tXMin,tXMax);
  }
  
  return 1;
}

void ExpCF1DHBT::Build(int aChargeType)
{
  Normalize();
  ApplyPurCorr(aChargeType);
  //  applyMomResCorr(aMomRes);
  FillFitArray();
}

void ExpCF1DHBT::Normalize()
{
  int tiFirstScale = mHCF->FindBin(mLowKStarNorm);
  int tiLastScale = mHCF->FindBin(mHighKStarNorm);
  double tScale = mHCF->Integral(tiFirstScale,tiLastScale)/
    (tiLastScale-tiFirstScale+1);
  mHCF->Scale(1./tScale);
}

void ExpCF1DHBT::FillFitArray()
{
  PRINT_DEBUG(" --- Fill variables used in chi2 function");
  PRINT_DEBUG("Process histo: " << mHCF->GetName());;

  double tHalfStep = (mHighKStarFit-mLowKStarFit)/(2.0*mNFitBin);
  int tFirstFitBin = mHCF->FindBin(mLowKStarFit+tHalfStep);
  int tLastFitBin = mHCF->FindBin(mHighKStarFit-tHalfStep);
  for(int ti=0; ti<mNFitBin; ti++){
    mContent[ti] = mHCF->GetBinContent(ti+tFirstFitBin);
    mErr2[ti] = mHCF->GetBinError(ti+tFirstFitBin)*
      mHCF->GetBinError(ti+tFirstFitBin);
  }
}

STR ExpCF1DHBT::MakeHistName(int aType, int aSystem)
{
  STR buff;

  buff = "";

  for (int iter=0; (mOrder[iter] && (iter<4)); iter++)
    {
      switch(mOrder[iter])
	{
	case EXPCF_PART_TYPE:
	  buff += mTypes[aType];
	  break;
	case EXPCF_PART_SYSTEM:
	  buff +=mSystems[aSystem];
	  break;
	case EXPCF_PART_GENERAL:
	  buff +=mGeneral;
	  break;
	}
    }
  buff += '\0';
  return buff;
}

STR ExpCF1DHBT::MakePurityHistName(int aSystem)
{
  STR buff;
  
  for (int iter=0; (mPOrder[iter] && (iter<4)); iter++)
    {
      switch(mPOrder[iter])
	{
	case EXPCF_PURITY_PART_PREFIX:
	  buff +=mPPrefix;
	  break;
	case EXPCF_PURITY_PART_SYSTEM:
	  buff +=mPSystems[aSystem];
	  break;
	case EXPCF_PURITY_PART_SUFFIX:
	  buff +=mPSuffix;
	  break;
	}
    }
  buff += '\0';
  return buff;  
}

STR ExpCF1DHBT::MakeFileName(int aChargeType)
{
  STR fName;
  fName="";
  try
    {
      fName = sRPInstance->getPar("InFileNamePrefix");
      switch (aChargeType)
	{
	case EXPCF_SYSTEM_PP:
	  fName += sRPInstance->getPar("InFileNamePositivePositivePart");
	  break;
	case EXPCF_SYSTEM_MM:
	  fName += sRPInstance->getPar("InFileNameNegativeNegativePart");
	  break;
	case EXPCF_SYSTEM_ZZ:
	  fName += sRPInstance->getPar("InFileNameNeutralNeutralPart");
	  break;
	}
      fName += sRPInstance->getPar("InFileNameSuffix");
      fName += '\0';
    }
  catch (STR e)
    {
      PRINT_MESSAGE("ExpCf <MakeFileName>: Error creating filename " << e);
    }
  return fName;
}

STR ExpCF1DHBT::MakePurityFileName(int aChargeType)
{
  STR fName;
  try
    {
      fName = sRPInstance->getPar("InPurityFileNamePrefix");
      switch (aChargeType)
	{
	case EXPCF_SYSTEM_PP:
	  fName += sRPInstance->getPar("InPurityFileNamePositivePositivePart");
	  break;
	case EXPCF_SYSTEM_MM:
	  fName += sRPInstance->getPar("InPurityFileNameNegativeNegativePart");
	  break;
	case EXPCF_SYSTEM_ZZ:
	  fName += sRPInstance->getPar("InPurityFileNeutralNeutralPart");
	  break;
	}
      fName += sRPInstance->getPar("InPurityFileNameSuffix");
      fName += '\0';
      PRINT_DEBUG("Got purity filename " << fName.Data());
    }
  catch (STR e)
    {
      PRINT_MESSAGE("ExpCF <MakePurityFileName>: Purity filename parameters not found " << e);
      PRINT_MESSAGE("Using regular input filename instead");
      return MakeFileName(aChargeType);
    }
  return fName;
}
void           
ExpCF1DHBT::GenerateParameterStub()
{
  fstream *os;
  os = new fstream(sRPFileName.Data(),ios::ate|ios::out|ios::in);

  (*os) 
    << "# The parameters for the experimental correlation functions" << endl
    << "#   for 1D identical particle particle HBT correlations" << endl
    << "# The percentage of the momentum resolution correction" << endl
    << "MomentumResolutionCorrection	= 0" << endl
    << "" << endl
    << "# The percentage of the purity correction (0 - do not apply correction)" << endl
    << "PurityCorrection	= 0" << endl
    << "" << endl
    << "# Uncomment and fill this lines if the file containing purity" << endl
    << "# histograms is different than InFileName" << endl
    << "# InPurityFileNamePrefix	= " << endl
    << "# InPurityFileNameSuffix	= " << endl
    << "# InPurityFileNamePositivePositivePart	= " << endl
    << "# InPurityFileNameNegativeNegativePart	= " << endl
    << "# InPurityFileNameNeutralNeutralPart	= " << endl
    << endl
    << "# The section with the name of the purity histograms" << endl
    << "# The purity histogram name consists of four parts" << endl
    << "#  - the prefix (identical for all)" << endl
    << "#  - the system name" << endl
    << "#  - the suffix (identical for all)" << endl
    << "# in an arbitrary order" << endl
    << endl
    << "PurityPartOrder = Prefix,System,Suffix" << endl
    << endl
    << "PurityPrefixPart = " << endl
    << "PuritySuffixPart = " << endl
    << endl
    << "PurityPositivePositivePart	= " << endl
    << "PurityNegativeNegativePart	= " << endl
    << "PurityNeutralNeutralPart	= " << endl
    << endl
    << "# The variable in which the purity histogram was made" << endl
    << "# possible values are:" << endl
    << "# KStar, 2KStar, Qinv" << endl
    << "PurityHistogramVariable	= Qinv" << endl
    << endl
    << "# The section describing the name of the file containig" << endl
    << "# all the neccessary histograms" << endl
    << "# The Filename consists of three parts:" << endl
    << "# - the prefix (identical for all files)" << endl
    << "# - the pair specific part - corresponding to" << endl
    << "#   the part of the name that changes from " << endl
    << "#   one charge combination to the other" << endl
    << "#   If all pair types are in one file," << endl
    << "#   than all these parts will be the same" << endl
    << "# - the suffix (identical for all files)" << endl
    << "InFileNamePrefix = " << endl
    << "InFileNameSuffix = " << endl
    << "InFileNamePositivePositivePart = " << endl
    << "InFileNameNegativeNegativePart = " << endl
    << "InFileNameNeutralNeutralPart = " << endl
    << endl
    << "# The section describing the histogram names in the input file" << endl
    << "# The name consists of five sections:" << endl
    << "#	General - genereal part of the name, the same for all histograms" << endl
    << "#	Type	- type of the histogram (numerator or denomninator) " << endl
    << "#	System	- the system name (e.g. PimKp, PP, etc.)" << endl
    << "# Example:	RatPipPip1DHBT" << endl
    << "#		^^^      ^^^^^" << endl
    << "#		Type     General " << endl
    << "#		    ^^^^^" << endl
    << "#		    System" << endl
    << "# First we say what is the order in which the above parts go" << endl
    << "PartOrder	= Type,System,General" << endl
    << endl
    << "#Next we give the values for each of the parts" << endl
    << "# First the Type parts" << endl
    << "NumeratorPart	= Num" << endl
    << "DenominatorPart	= Den" << endl
    << endl
    << "# the System Part" << endl
    << "PositivePositivePart	= " << endl
    << "NegativeNegativePart	= " << endl
    << "NeutralNeutralPart	= " << endl
    << endl
    << "# the General Part" << endl
    << "GeneralPart	= " << endl
    << endl
    << "# The variable in which the data histograms were made" << endl
    << "# possible values are:" << endl
    << "# KStar, 2KStar, Qinv" << endl
    << "DataHistogramVariable	= Qinv" << endl << endl
    << "# Fitting range parameters" << endl
    << "# The range of the fitting (in k*)" << endl
    << "FitRangeMin	= " << endl
    << "FitRangeMax	= " << endl
    << endl
    << "# The range where the function is normalized" << endl
    << "NormalizationRangeMin = " << endl
    << "NormalizationRangeMax = " << endl
    << endl;
  
  os->close();
}

void ExpCF1DHBT::ReadParameters()
{
  STR fVar;

  try 
    {
      // For now - only the non-id 3D case is supported
      // First read the order
      char *fOrder = strdup(sRPInstance->getPar("PartOrder").Data());
      PRINT_DEBUG("Got order: " << mOrder);
    
      char *next = strstr(fOrder,",");
      char *curr = fOrder;
      int count=0;
      // First clear the order table
      for (int yyy=0; yyy<5; yyy++)
	mOrder[yyy] = 0;

      while (curr != next)
	{
	  if (!strncmp(curr,"Type",4))
	    mOrder[count] = EXPCF_PART_TYPE;
	  else if (!strncmp(curr,"System",6))
	    mOrder[count] = EXPCF_PART_SYSTEM;
	  else if (!strncmp(curr,"General",7))
	    mOrder[count] = EXPCF_PART_GENERAL;
	  else
	    {
	      PRINT_MESSAGE("ExpCF <readParameters>: unknown order part: " << curr);
	      exit(1);
	    }
	  if (next)
	    {
	      curr = next+1;
	      next = strstr(curr,",");
	    }
	  else
	    curr = next;
	  count++;
	}
      PRINT_DEBUG("Dishing out the order: ");
      for (int yyy=0; yyy<5; yyy++)
	PRINT_DEBUG(mOrder[yyy] << " ");
      
      mTypes[EXPCF_TYPE_NUM] = sRPInstance->getPar("NumeratorPart");
      mTypes[EXPCF_TYPE_DEN] = sRPInstance->getPar("DenominatorPart");
      mSystems[EXPCF_SYSTEM_PP] = sRPInstance->getPar("PositivePositivePart");
      mSystems[EXPCF_SYSTEM_MM] = sRPInstance->getPar("NegativeNegativePart");
      mSystems[EXPCF_SYSTEM_ZZ] = sRPInstance->getPar("NeutralNeutralPart");
      mGeneral = sRPInstance->getPar("GeneralPart");
      mPurCorr = atof((sRPInstance->getPar("PurityCorrection")).Data())/100.0;
      mMomRes  = atof((sRPInstance->getPar("MomentumResolutionCorrection")).Data())/100.0;

      if (mPurCorr!=0)
	{
	  // First clear the purity order
	  for (int yyy=0; yyy<5; yyy++)
	    mPOrder[yyy]= 0;
	  
	  fOrder = strdup(sRPInstance->getPar("PurityPartOrder").Data());
	  PRINT_DEBUG("Got order: " << fOrder);
	  
	  next = strstr(fOrder,",");
	  curr = fOrder;
	  count=0;
	  while (curr != next)
	    {
	      if (!strncmp(curr,"Prefix",6))
		mPOrder[count] = EXPCF_PURITY_PART_PREFIX;
	      else if (!strncmp(curr,"System",6))
		mPOrder[count] = EXPCF_PURITY_PART_SYSTEM;
	      else if (!strncmp(curr,"Suffix",6))
		mPOrder[count] = EXPCF_PURITY_PART_SUFFIX;
	      else
		{
		  PRINT_MESSAGE("ExpCF <readParameters>: unknown purity order part: " << curr);
		  exit(1);
		}
	      if (next)
		{
		  curr = next+1;
		  next = strstr(curr,",");
		}
	      else
		curr = next;
	      count++;
	    }
	  PRINT_DEBUG("Dishing out the purity order: ");
	  for (int yyy=0; yyy<5; yyy++)
	    PRINT_DEBUG(mPOrder[yyy] << " ");
	  
	  mPSystems[EXPCF_SYSTEM_PP] = sRPInstance->getPar("PurityPositivePositivePart");
	  mPSystems[EXPCF_SYSTEM_MM] = sRPInstance->getPar("PurityNegativeNegativePart");
	  mPSystems[EXPCF_SYSTEM_ZZ] = sRPInstance->getPar("PurityNeutralNeutralPart");
	  mPPrefix = sRPInstance->getPar("PurityPrefixPart");
	  mPSuffix = sRPInstance->getPar("PuritySuffixPart");
	}
      mLowKStarFit   = atof((sRPInstance->getPar("FitRangeMin")).Data());
      mHighKStarFit  = atof((sRPInstance->getPar("FitRangeMax")).Data());
      mLowKStarNorm  = atof((sRPInstance->getPar("NormalizationRangeMin")).Data());
      mHighKStarNorm = atof((sRPInstance->getPar("NormalizationRangeMax")).Data());
      
      fVar = sRPInstance->getPar("PurityHistogramVariable");
      if (fVar == "KStar")
	mPurityKStar = 1;
      else
	mPurityKStar = 0;
      if ((fVar != "KStar") && (fVar != "2KStar") && (fVar != "Qinv"))
	PRINT_MESSAGE("Unknown purity histogram variable " << fVar << " assuming KStar\n");
      
      
      fVar = sRPInstance->getPar("DataHistogramVariable");
      if (fVar == "KStar")
	mDataKStar = 1;
      else
	mDataKStar = 0;
      if ((fVar != "KStar") && (fVar != "2KStar") && (fVar != "Qinv"))
	PRINT_MESSAGE("Unknown data histogram variable " << fVar << " assuming KStar\n");
      
    }
  
  catch (STR e)
    {
      PRINT_MESSAGE("Error reading parameters in ExpCF <readParameters>: " << e);
    }
}

void ExpCF1DHBT::SetParameters(int aNFitBin, double aLowKStarFit, double aHighKStarFit,
			       double aLowKStarNorm, double aHighKStarNorm)
{
  if(mContent){
    PRINT_DEBUG("parameters can only be set once");
  }
  else{    
    mNFitBin=aNFitBin;
    mContent = new double[mNFitBin];
    mErr2 = new double[mNFitBin];
    mLowKStarFit=aLowKStarFit;
    mHighKStarFit=aHighKStarFit;
    mLowKStarNorm=aLowKStarNorm;
    mHighKStarNorm=aHighKStarNorm;
    mKStarSigned = new double[mNFitBin];
    double tStep = (mHighKStarFit-mLowKStarFit)/mNFitBin;
    for(int ti=0; ti<mNFitBin; ti++){
      mKStarSigned[ti] = mLowKStarFit + (ti+0.5)*tStep;
    }
    PRINT_DEBUG( "Set parameters. Fit: " << mNFitBin << " , " <<  mLowKStarFit
	<< " < k* < " << mHighKStarFit << " | Norm: " << mLowKStarNorm
	<< " < k* < " << mHighKStarNorm);
  }
}

double ExpCF1DHBT::GetLowKStarFit()
{
  return mLowKStarFit;
}

double ExpCF1DHBT::GetHighKSstarFit()
{
  return mHighKStarFit;
}

double ExpCF1DHBT::GetLowKStarNorm()
{
  return mLowKStarNorm;
}

double ExpCF1DHBT::GetHighKStarNorm()
{
  return mHighKStarNorm;
}

