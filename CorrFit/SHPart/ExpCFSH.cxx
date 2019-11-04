#include "fstream"
#include "stdlib.h"
#include "ExpCFSH.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TKey.h"
#include "TROOT.h"
#include "ReadPar.h"
#include "TF1.h"

extern ReadPar *sRPInstance;
extern STR      sRPFileName;

using namespace std;


ExpCFSH::ExpCFSH() : ExpCF()
{
}

ExpCFSH::ExpCFSH(int aChargeType,
		 TFile* aOutFile) : ExpCF()
{
  ReadParameters();
  int tValid = ReadCF(aChargeType,aOutFile);
  if(tValid){
    int tiFirstFitBin = mRC00->GetXaxis()->FindFixBin(mLowKStarFit);
    int tiLastFitBin  = mRC00->GetXaxis()->FindFixBin(mHighKStarFit);
    int tiFirstScale = mRC00->FindBin(mLowKStarNorm);
    if(tiLastFitBin==tiFirstScale) tiFirstScale++;
    int tiLastScale = mRC00->FindBin(mHighKStarNorm);
    PRINT_MESSAGE("Setting parameters:");
    PRINT_MESSAGE("Fit bins. First: " << tiFirstFitBin << " Last: " << tiLastFitBin);
    PRINT_MESSAGE("Norm bins. First: " << tiFirstScale << " Last: " << tiLastScale);
    SetParameters((tiLastFitBin-tiFirstFitBin+1),
		  mRC00->GetBinLowEdge(tiFirstFitBin),
		  mRC00->GetBinLowEdge(tiLastFitBin+1),
		  mRC00->GetBinLowEdge(tiFirstScale),
		  mRC00->GetBinLowEdge(tiLastScale+1));
    Build(aChargeType);
  }
}

ExpCFSH::~ExpCFSH()
{
  if (mPurity) 
    delete [] mPurity;
}

TGraphErrors* ExpCFSH::GetGraph(const char* aName)
{
}

TGraphErrors* ExpCFSH::GetGraph(const char* aName, double aPur)
{
}

double ExpCFSH::GetContent(int aIndex, double aPurity) const
{
  double tVar = mContent[aIndex];
  return (tVar-1.0) * aPurity + 1.0;
}

double ExpCFSH::GetContent(int aIndex, double aPurity, double aNorm) const
{
  double tVar = mContent[aIndex]*aNorm;
  return (tVar-1.0) / (aPurity * mPurity[aIndex]) + 1.0;
}

void ExpCFSH::WriteHisto()
{
}

void ExpCFSH::Write()
{
//   mHHCF->Write();
//   mHCF->Write();

//   PRINT_DEBUG_2("Experimental fitted best-fit CF");
//   for (int iter=1; iter<=mHCF->GetNbinsX(); iter++)
//     PRINT_DEBUG_2((mHCF->GetXaxis()->GetBinCenter(iter)) << " " << (mHCF->GetBinContent(iter)) << " +/- " << (mHCF->GetBinError(iter)));
  
  mRC00->Write();
  mRC11->Write();
}

void ExpCFSH::ApplyPurCorr(int aChargeType)
{
  if(mPurCorr==0) return;

//   TFile* tFIn = new TFile(MakePurityFileName(aChargeType).Data());
//   if(!tFIn->IsOpen()) return;

//   int tIAxis = EXPCFSH_PROJ_OUT;
  
//   //  for(int ti=0; ti<3; ti++){
//     TH1D* tHCF = mHCF;

//     TH2D* tHPur;

//     // First check if the purity histogram is a profile
//     // or a TH2D
//     TKey *ktemp = tFIn->FindKey(MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data());
//     if (ktemp)
//       {
// 	TProfile *ttPPur;
// 	if (!strcmp(ktemp->GetClassName(),"TProfile"))
// 	  {
// 	    // It's a TProfile!
// 	    // No need to do any fancy rewriting, just read it
	    
// 	    if (mPType == EXPCFSH_PURITY_TYPE_SINGLE)
// 	      {
// 		ttPPur = (TProfile *) tFIn->Get(MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data());
// 	      }
	    
// 	    else if (mPType == EXPCFSH_PURITY_TYPE_DOUBLE)
// 	      {
// 		TProfile* tHP = (TProfile *) tFIn->Get(MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data());
// 		TProfile* tHN = (TProfile *) tFIn->Get(MakePurityHistName(tIAxis,EXPCFSH_SIGN_NEG,aChargeType).Data());
// 		PRINT_DEBUG_2("Getting purity profiles " << 
// 			      MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType) << " " << tHP << " " <<
// 			      MakePurityHistName(tIAxis,EXPCFSH_SIGN_NEG,aChargeType) << " " << tHN);
			      
// 		int tNX = tHN->GetNbinsX();
// 		ttPPur = new TProfile("ttPPur","ttPPur",
// 				      2*tNX,-tHN->GetXaxis()->GetXmax(),
// 				      tHN->GetXaxis()->GetXmax(),
// 				      0.,1.);
// 		for(int ti=1; ti<=tNX; ti++){
//  		  ttPPur->SetBinEntries(ti+tNX,tHP->GetBinEntries(ti));
//  		  ttPPur->SetBinEntries(tNX-ti+1,tHN->GetBinEntries(ti));
// 		  ttPPur->SetBinContent(ti+tNX,tHP->GetBinContent(ti)*tHP->GetBinEntries(ti));
// 		  ttPPur->SetBinContent(tNX-ti+1,tHN->GetBinContent(ti)*tHN->GetBinEntries(ti));
// 		  ttPPur->SetBinError(ti+tNX,tHP->GetBinError(ti));
// 		  ttPPur->SetBinError(tNX-ti+1,tHN->GetBinError(ti));
// 		  PRINT_DEBUG_3("Setting bin " << ti << " to " << ttPPur->GetBinContent(ti+tNX) << " " << ttPPur->GetBinContent(tNX-ti+1));
// 		}
// 	      }
	    
// 	    //     }
// 	    //	    if(!tHPur){ continue;}
// 	  }
// 	else if (!strcmp(ktemp->GetClassName(),"TH2D"))
// 	  {
// 	    // It's a TH2D 
// 	    // We need to rewrite it and cast it to a TProfile
// 	    TH2D *ttHPur;
	    
// 	    if (mPType == EXPCFSH_PURITY_TYPE_SINGLE)
// 	      {
// 		tHPur= (TH2D*) tFIn->Get(MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data());
// 		PRINT_DEBUG_2("Try to get " << MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType) << " " << tHPur);;
//  		int tNBin = tHPur->GetNbinsX()/2;
// // 		ttHPur = new TH2D("ttHPur","ttHPur",tNBin,0.,
// // 				  tHPur->GetXaxis()->GetXmax(),tHPur->GetNbinsY(),0.,1.);
// // 		for(int ti=1; ti<=tNBin; ti++){
// // 		  for(int tj=1; tj<=ttHPur->GetNbinsY(); tj++){
// // 		    ttHPur->SetCellContent(ti,tj,tHPur->GetCellContent(tNBin+ti,tj)+
// // 					   tHPur->GetCellContent(tNBin-ti+1,tj));
// // 		  }
// // 		}
//  		ttHPur = new TH2D(*tHPur);
// 		ttHPur->SetTitle("ttHPur");
// 		ttHPur->SetName("ttHPur");
		
// 	      }
	    
// 	    else if (mPType == EXPCFSH_PURITY_TYPE_DOUBLE)
// 	      {
// 		TH2D* tHN = (TH2D*) tFIn->Get(MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data());
// 		TH2D* tHP = (TH2D*) tFIn->Get(MakePurityHistName(tIAxis,EXPCFSH_SIGN_NEG,aChargeType).Data());
// 		ttHPur = (TH2D*) gROOT->FindObject("ttHPur");
// 		if(!ttHPur){	
// 		  int tNX = tHN->GetNbinsX();
// 		  ttHPur = new TH2D("ttHPur","ttHPur",
// 				    2*tNX,-tHN->GetXaxis()->GetXmax(),
// 				    tHN->GetXaxis()->GetXmax(),
// 				    tHN->GetNbinsY(),0.,1.);
// 		  for(int ti=1; ti<=tNX; ti++){
// 		    for(int tj=1; tj<=tHN->GetNbinsY(); tj++){
// 		      ttHPur->SetCellContent(ti+tNX,tj,tHP->GetCellContent(ti,tj));
// 		      ttHPur->SetCellContent(tNX-ti+1,tj,tHN->GetCellContent(ti,tj));
// 		    }
// 		  }
// 		}
// 	      }
	    
// 	    //     }
// 	    if(!ttHPur){ return;}
// 	    // ---
// 	    ttPPur = ttHPur->ProfileX("ttPPur",1,ttHPur->GetNbinsY(),"E");
// 	  }
// 	else 
// 	  {
// 	    // It's neither a TProfile nor TH2D
// 	    // We bail out then
// 	    PRINT_MESSAGE("The key " << MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType) << 
// 			  " is neither TProfile nor TH2D - cannot use" << endl <<
// 			  "Not applying Purity correction");
// 	    return;
// 	  }
	
// 	if(!mPurityKStar){
// 	  PRINT_DEBUG_2("Rebinning KStar");
// 	  int tNBin = ttPPur->GetNbinsX();
// 	  double tXMin = ttPPur->GetXaxis()->GetXmin()/2.;
// 	  double tXMax = ttPPur->GetXaxis()->GetXmax()/2.;
// 	  ttPPur->GetXaxis()->Set(tNBin,tXMin,tXMax);
// 	}
  
// 	double tCorr;
// 	double tSignal;

// // 	for(int ti=1;ti<=tHCF->GetNbinsX();ti++){
// // 	  tCorr = mPurCorr*ttPPur->GetBinContent(ti);
// // 	  PRINT_DEBUG_3("Purity bin value: " << ttPPur->GetBinContent(ti));

// // 	  if(tCorr>1.) tCorr=1.;
// // 	  tSignal = tHCF->GetBinContent(ti)-1.0;
// // 	  PRINT_DEBUG_3("Signal in bin " << ti << " is " << tSignal);

// // 	  tHCF->SetBinContent(ti,(tSignal/tCorr)+1.0);
// // 	  PRINT_DEBUG_3("Setting bin " << ti << " to " << tHCF->GetBinContent(ti));

// // 	  tHCF->SetBinError(ti,tHCF->GetBinError(ti)/tCorr);
// // 	  PRINT_DEBUG_2("Applying purity correction " << tCorr);
	  
// // 	}

// // 	for(int ti=1;ti<=mHHCF->GetNbinsX();ti++){
// // 	  tCorr = mPurCorr*ttPPur->GetBinContent(ti);
// // 	  if(tCorr>1.) tCorr=1.;
// // 	  tSignal = mHHCF->GetBinContent(ti)-1.;
// // 	  mHHCF->SetBinContent(ti,tSignal/tCorr+1);
// // 	  mHHCF->SetBinError(ti,mHHCF->GetBinError(ti)/tCorr);
// // 	  PRINT_DEBUG_2("Applying purity correction " << tCorr);
// // 	  PRINT_DEBUG_3("Setting bin " << ti << " to " << mHHCF->GetBinContent(ti));
	  
// // 	}

// 	double tHalfStep = (mHighKStarFit-mLowKStarFit)/mNFitBin;
// 	int tFirstFitBin = tHCF->FindBin(-1.*mHighKStarFit+tHalfStep);
// 	int tLastFitBin = tHCF->FindBin(-1.*mLowKStarFit-tHalfStep);    
// 	for(int ti=0; ti<mNFitBin/2; ti++){
// 	  mPurity[ti] = ttPPur->GetBinContent(ti+tFirstFitBin);
// 	}
 
// 	tFirstFitBin = tHCF->FindBin(mLowKStarFit+tHalfStep);
// 	tLastFitBin = tHCF->FindBin(mHighKStarFit-tHalfStep);
// 	for(int ti=0; ti<mNFitBin/2; ti++){
// 	  mPurity[ti+mNFitBin/2] = ttPPur->GetBinContent(ti+tFirstFitBin);
// 	}

// 	ttPPur->Delete();  
//       }
//     else
//       {
// 	PRINT_MESSAGE("Purity histogram " << MakePurityHistName(tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data() << " not found");
//       }
//     //  }
//   tFIn->Delete();
}

void ExpCFSH::ApplyMomResCorr(double aMomResCorr)
{
}

int ExpCFSH::ReadCF(int aChargeType, TFile* aOutFile)
{
  // A new code using the input parameter class

  TFile *tFIn=0;
  STR fName;
  
  PRINT_DEBUG_2("Reading the experimental CF");
  
  fName = MakeFileName(aChargeType);
  
  PRINT_DEBUG_2("Got FileName: " << fName);
  
  if(!tFIn){
    PRINT_DEBUG_2("Opening " << fName << " " << tFIn );
    tFIn = new TFile(fName.Data());
    PRINT_DEBUG_2("Opened " << tFIn);
  }
  if(!tFIn->IsOpen()) return 0;

//   TH1D* tHCFN=0;
//   TH1D* tHCFP=0;
//   TH1D* tHDenN=0;
//   TH1D* tHDenP=0;
//   TH1D* tHNumN=0;
//   TH1D* tHNumP=0;

//   int tNAxis = 1;
//   int tIAxis;

  
  PRINT_DEBUG_2("Reading " << MakeHistName(0,0,aChargeType).Data());
  TKey *ktemp = tFIn->FindKey(MakeHistName(0,0,aChargeType).Data());
  if (ktemp) {
    mRC00 = new TH1D(*((TH1D *) tFIn->Get(MakeHistName(0,0,aChargeType).Data())));
  }
  else {
    return 0;
  }

  PRINT_DEBUG_2("Reading " << MakeHistName(0,1,aChargeType).Data());
  ktemp = tFIn->FindKey(MakeHistName(0,1,aChargeType).Data());
  if (ktemp) {
    mRC11 = new TH1D(*((TH1D *) tFIn->Get(MakeHistName(0,1,aChargeType).Data())));
  }
  else {
    return 0;
  }

  PRINT_DEBUG_2("Reading " << MakeHistName(1,2,aChargeType).Data());
  ktemp = tFIn->FindKey(MakeHistName(1,2,aChargeType).Data());
  if (ktemp) {
    mCov = new TH3D(*((TH3D *) tFIn->Get(MakeHistName(1,2,aChargeType).Data())));
  }
  else {
    return 0;
  }

//   for( tIAxis=0; tIAxis<tNAxis; tIAxis++){
//     int ti=0;
//     tHNumN=0;
//     tHNumP=0;
//     char *temp;
    
//     TKey *ktemp = tFIn->FindKey(MakeHistName(EXPCFSH_TYPE_DEN,tIAxis,EXPCFSH_SIGN_NEG,aChargeType).Data());
//     if (ktemp)
//       {
// 	if (!strcmp(ktemp->GetClassName(),"TProfile")) {
// 	  tHDenN = ((TProfile*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_DEN,tIAxis,EXPCFSH_SIGN_NEG,aChargeType).Data()))->ProjectionX();
// 	  tHDenP = ((TProfile*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_DEN,tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data()))->ProjectionX();
// 	  tHNumN  = ((TProfile*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_NUM,tIAxis,EXPCFSH_SIGN_NEG,aChargeType).Data()))->ProjectionX();
// 	  tHNumP  = ((TProfile*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_NUM,tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data()))->ProjectionX();
// 	}
// 	else {
// 	  tHDenN = (TH1D*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_DEN,tIAxis,EXPCFSH_SIGN_NEG,aChargeType).Data());
// 	  tHDenP = (TH1D*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_DEN,tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data());
// 	  tHNumN  = (TH1D*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_NUM,tIAxis,EXPCFSH_SIGN_NEG,aChargeType).Data());
// 	  tHNumP  = (TH1D*) tFIn->Get(MakeHistName(EXPCFSH_TYPE_NUM,tIAxis,EXPCFSH_SIGN_POS,aChargeType).Data());

// 	}
// 	if (tIAxis == EXPCFSH_PROJ_OUT){
// 	  Double_t tFitMin = tHDenP->GetBinLowEdge(tHDenP->FindBin(mLowKStarFit));
// 	  Double_t tFitMax = tHDenP->GetBinLowEdge(tHDenP->FindBin(mHighKStarFit)+1);
// 	  TF1 *poll = new TF1("poll","[0]*x*x*x+[1]*x*x+[2]*x+[3]");
// 	  tHDenP->Fit("poll","QNO","",tFitMin,tFitMax);
// 	  mHDenOutPFit = new TF1(*poll);
// 	  tHDenN->Fit("poll","QNO","",tFitMin,tFitMax);
// 	  mHDenOutNFit = new TF1(*poll);
// 	}
//       }


//     PRINT_DEBUG_2("Try to get " << MakeHistName(EXPCFSH_TYPE_NUM,tIAxis,EXPCFSH_SIGN_NEG,aChargeType) << " " << tHNumN << " " << MakeHistName(EXPCFSH_TYPE_NUM,tIAxis,EXPCFSH_SIGN_POS,aChargeType) << " " << tHNumP);
    
//     ti++;
//     if (tHNumN==0 && tHNumP==0) {
//       PRINT_DEBUG_1("Function " << MakeHistName(EXPCFSH_TYPE_NUM,tIAxis,EXPCFSH_SIGN_NEG,aChargeType) << " not found!");
//       // *** FIXME ***
//       // Determine which projection we want to fit from the parameters
//       // for now assume Out - mostly correct
//       if (tIAxis == 0) {
// 	PRINT_MESSAGE("Cannot fit without " << mProjs[EXPCFSH_PROJ_OUT] << " projection!");
// 	return 0;
//       }
//     }
    
//     if(!tHNumN || !tHNumP){ return 0;}
    
//     aOutFile->cd();

//     TH1D *tHHNum;
//     if(tIAxis==0){ // calc cf
//       mHHCF = new TH1D(*tHNumN);
//       tHHNum = new TH1D(*tHNumN);
//       STR tName;
//       tName = "HCF";
//       if (mSystems[aChargeType] != "")
// 	tName += mSystems[aChargeType];
//       else
// 	switch (aChargeType)
// 	  {
// 	  case EXPCFSH_SYSTEM_PP:
// 	    tName += "_PP_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_PM:
// 	    tName += "_PM_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_MP:
// 	    tName += "_MP_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_MM:
// 	    tName += "_MM_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_ZZ:
// 	    tName += "_ZZ_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_ZP:
// 	    tName += "_ZP_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_ZM:
// 	    tName += "_ZM_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_PZ:
// 	    tName += "_PZ_";
// 	    break;
// 	  case EXPCFSH_SYSTEM_MZ:
// 	    tName += "_MZ_";
// 	    break;
// 	  }
//       tName += '\0';
//       mHHCF->SetName(tName.Data());
//       tHHNum->Add(tHNumP);
//       TH1D tHDen(*tHDenN);
//       tHDen.SetName("tHDen");
//       tHDen.Add(tHDenP);
//       if (mBinominalErrors)
// 	mHHCF->Divide(tHHNum,&tHDen,1.0,1.0,"B");
//       else
// 	mHHCF->Divide(tHHNum,&tHDen,1.0,1.0);
//     }

//     tHCFN = new TH1D(*tHNumN);
//     tHCFP = new TH1D(*tHNumP);
//     if (mBinominalErrors) {
//       tHCFN->Divide(tHNumN,tHDenN,1.0,1.0,"B");
//       tHCFP->Divide(tHNumP,tHDenP,1.0,1.0,"B");
//     }
//     else {
//       tHCFN->Divide(tHNumN,tHDenN,1.0,1.0);
//       tHCFP->Divide(tHNumP,tHDenP,1.0,1.0);
//     }

//     int tNBin = tHCFN->GetNbinsX();
//     double tXmax = tHCFN->GetXaxis()->GetXmax();
//     STR tName;
//     tName = "HCF";
//     if (mSystems[aChargeType] != "")
//       tName += mSystems[aChargeType];
//     else
//       switch (aChargeType)
// 	{
// 	case EXPCFSH_SYSTEM_PP:
// 	  tName += "_PP_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_PM:
// 	  tName += "_PM_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_MP:
// 	  tName += "_MP_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_MM:
// 	  tName += "_MM_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_ZZ:
// 	  tName += "_ZZ_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_ZP:
// 	  tName += "_ZP_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_ZM:
// 	  tName += "_ZM_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_PZ:
// 	  tName += "_PZ_" ;
// 	  break;
// 	case EXPCFSH_SYSTEM_MZ:
// 	  tName += "_MZ_" ;
// 	  break;
// 	}
//     if (mProjs[tIAxis] != "")
//       tName += mProjs[tIAxis];
//     else
//       switch (aChargeType)
// 	{
// 	case EXPCFSH_PROJ_OUT:
// 	  tName += "_OUT_" ;
// 	  break;
// 	case EXPCFSH_PROJ_SIDE:
// 	  tName += "_SIDE_" ;
// 	  break;
// 	case EXPCFSH_PROJ_LONG:
// 	  tName += "_LONG_" ;
// 	  break;
// 	}
//     tName += '\0';
      
//     TH1D* tHCF = new TH1D(tName.Data(),tName.Data(),2*tNBin,-tXmax,tXmax);
      
//     for(int ti=1;ti<=tNBin;ti++){
//       tHCF->SetBinContent(ti,tHCFN->GetBinContent(tNBin-ti+1));
//       tHCF->SetBinError(ti,tHCFN->GetBinError(tNBin-ti+1));
//       tHCF->SetBinContent(ti+tNBin,tHCFP->GetBinContent(ti));
//       tHCF->SetBinError(ti+tNBin,tHCFP->GetBinError(ti));
//     }
//     mHCF = new TH1D(*tHCF);
//   }
//  tFIn->Delete();
  
//   if(!mDataKStar){
//     int tNBin = mHHCF->GetNbinsX();
//     double tXMin = mHHCF->GetXaxis()->GetXmin()/2.;
//     double tXMax = mHHCF->GetXaxis()->GetXmax()/2.;
//     mHHCF->GetXaxis()->Set(tNBin,tXMin,tXMax);

//     tNBin = mHCF->GetNbinsX();
//     tXMin = mHCF->GetXaxis()->GetXmin()/2.;
//     tXMax = mHCF->GetXaxis()->GetXmax()/2.;

//     mHCF->GetXaxis()->Set(tNBin,tXMin,tXMax);
//   }

  return 1;
}

void ExpCFSH::Build(int aChargeType)
{
  Normalize();
  ApplyPurCorr(aChargeType);
  //  applyMomResCorr(aMomRes);
  FillFitArray();
}

void ExpCFSH::Normalize()
{
//   {
//     int tiFirstScale = mHHCF->FindBin(mLowKStarNorm);
//     int tiLastScale = mHHCF->FindBin(mHighKStarNorm);
//     double tScale = mHHCF->Integral(tiFirstScale,tiLastScale)/
//       (tiLastScale-tiFirstScale+1);
//     mHHCF->Scale(1./tScale);
//   }
  
//   int tiFirstScale = mHCF->FindBin(mLowKStarNorm);
//   int tiLastScale = mHCF->FindBin(mHighKStarNorm)-1;
// //   mNormP = mHCF->Integral(tiFirstScale,tiLastScale)/
// //     (tiLastScale-tiFirstScale+1);
//   mNormP = mHCF->Integral(tiFirstScale,tiLastScale);

//   tiFirstScale = mHCF->FindBin(-1.*mHighKStarNorm);
//   tiLastScale = mHCF->FindBin(-1.*mLowKStarNorm)-1;
// //   mNormN = mHCF->Integral(tiFirstScale,tiLastScale)/
// //     (tiLastScale-tiFirstScale+1);
//   mNormN = mHCF->Integral(tiFirstScale, tiLastScale);

// //   int tNBinOver2 = mHCF->GetNbinsX()/2;
// //   for(int ti=1;ti<=tNBinOver2;ti++){
// //     mHCF->SetBinContent(ti,mHCF->GetBinContent(ti)/mNormN);
// //     mHCF->SetBinError(ti,mHCF->GetBinError(ti)/mNormN);
// //     mHCF->SetBinContent(ti+tNBinOver2,
// // 			mHCF->GetBinContent(ti+tNBinOver2)/mNormP);
// //     mHCF->SetBinError(ti+tNBinOver2,
// // 		      mHCF->GetBinError(ti+tNBinOver2)/mNormP);
// //   }
  

//   if (mCFReduceEdgeEffects)
//     {
//       PRINT_MESSAGE("Reducing edge effects in the CF");
//       TF1 *funpol = new TF1("funpol","[0]+[1]*x+[2]*x*x",0.0,0.5);
      
//       mHCF->Fit("funpol","QNO","",mLowKStarNorm,mHCF->GetXaxis()->GetXmax());
//       int nbinhalf = mHCF->GetNbinsX()/2;
//       for (int iter=1; iter <= nbinhalf; iter++)
// 	{
// 	  mHCF->SetBinContent(nbinhalf + iter, mHCF->GetBinContent(nbinhalf + iter) - (funpol->Eval(mHCF->GetXaxis()->GetBinCenter(nbinhalf + iter)) - 1.0));
// 	}
      
//       mHCF->Fit("funpol","QNO","",-mHCF->GetXaxis()->GetXmax(),-mLowKStarNorm);
//       for (int iter=0; iter < nbinhalf; iter++)
// 	{
// 	  mHCF->SetBinContent(nbinhalf - iter, mHCF->GetBinContent(nbinhalf - iter) - (funpol->Eval(mHCF->GetXaxis()->GetBinCenter(nbinhalf - iter)) - 1.0));
// 	}
//     }

}

void ExpCFSH::FillFitArray()
{
  PRINT_DEBUG_2(" --- Fill variables used in chi2 function");
  PRINT_DEBUG_2("Process histo: " << mRC00->GetName());;

  double tHalfStep = (mHighKStarFit-mLowKStarFit)/(2*mNFitBin);
  int tFirstFitBin = mRC00->GetXaxis()->FindFixBin(mLowKStarFit+tHalfStep);
  int tLastFitBin = mRC00->GetXaxis()->FindFixBin(mHighKStarFit-tHalfStep);    
  //  int tmb = (int) (mHighKStarNorm/((mHighKStarFit - mLowKStarFit)/mNFitBin));

  PRINT_DEBUG_3("Read in experimental function " << mNFitBin);
  for(int ti=0; ti<mNFitBin; ti++){
    mContent[ti] = mRC00->GetBinContent(ti+tFirstFitBin);
    mErr2[ti] = (mRC00->GetBinError(ti+tFirstFitBin)*
		 mRC00->GetBinError(ti+tFirstFitBin));

    PRINT_DEBUG_3("Bin val err " << ti << " " << mContent[ti] << " " << mErr2[ti]);
  }
 
//   tFirstFitBin = mHCF->FindBin(mLowKStarFit+tHalfStep);
//   tLastFitBin = mHCF->FindBin(mHighKStarFit-tHalfStep);
  for(int ti=0; ti<mNFitBin; ti++){
    mContent[mNFitBin+ti] = mRC11->GetBinContent(tFirstFitBin+ti);
    mErr2[mNFitBin+ti] = (mRC11->GetBinError(tFirstFitBin+ti)*
			  mRC11->GetBinError(tFirstFitBin+ti));
    PRINT_DEBUG_3("Bin val err " << ti+mNFitBin << " " << mContent[ti+mNFitBin] << " " << mErr2[ti+mNFitBin]);
  }

  for(int ti=0; ti<mNFitBin; ti++){
    mCovariance[ti] = (mCov->GetBinContent(tFirstFitBin+ti,1,7)/
		       (mRC11->GetBinError(tFirstFitBin+ti)*
			mRC00->GetBinError(tFirstFitBin+ti)));
    PRINT_DEBUG_3("Bin cov " << ti << " " << mCovariance[ti]);
  }
}

STR ExpCFSH::MakeHistName(int aType, int aProj, int aSystem)
{
  STR buff;

  buff = "";

  for (int iter=0; (mOrder[iter] && (iter<4)); iter++)
    {
      switch(mOrder[iter])
	{
	case EXPCFSH_PART_TYPE:
	  buff += mTypes[aType];
	  break;
	case EXPCFSH_PART_PROJ:
	  buff +=mProjs[aProj];
	  break;
	case EXPCFSH_PART_SYSTEM:
	  buff +=mSystems[aSystem];
	  break;
	case EXPCFSH_PART_GENERAL:
	  buff +=mGeneral;
	  break;
	}
    }
  buff += '\0';
  return buff;
}

STR ExpCFSH::MakePurityHistName(int aProj, int aSign, int aSystem)
{
  STR buff;
  
  for (int iter=0; (mPOrder[iter] && (iter<5)); iter++)
    {
      switch(mPOrder[iter])
	{
	case EXPCFSH_PURITY_PART_PREFIX:
	  buff +=mPPrefix;
	  break;
	case EXPCFSH_PURITY_PART_SYSTEM:
	  buff +=mPSystems[aSystem];
	  break;
	case EXPCFSH_PURITY_PART_PROJ:
	  buff +=mPProjs[aProj];
	  break;
// 	case EXPCFSH_PURITY_PART_SIGN:
// 	  if (mPType == EXPCFSH_PURITY_TYPE_DOUBLE)
// 	    buff +=mPSigns[aSign];
// 	  break;
	case EXPCFSH_PURITY_PART_SUFFIX:
	  buff +=mPSuffix;
	  break;
	}
    }
  buff += '\0';
  return buff;  
}

STR ExpCFSH::MakeFileName(int aChargeType)
{
  STR fName;
  fName="";
  try
    {
      fName = sRPInstance->getPar("InFileNamePrefix");
      switch (aChargeType)
	{
	case EXPCFSH_SYSTEM_PP:
	  fName += sRPInstance->getPar("InFileNamePositivePositivePart");
	  break;
	case EXPCFSH_SYSTEM_MM:
	  fName += sRPInstance->getPar("InFileNameNegativeNegativePart");
	  break;
	case EXPCFSH_SYSTEM_PM:
	  fName += sRPInstance->getPar("InFileNamePositiveNegativePart");
	  break;
	case EXPCFSH_SYSTEM_MP:
	  fName += sRPInstance->getPar("InFileNameNegativePositivePart");
	  break;
	case EXPCFSH_SYSTEM_ZZ:
	  fName += sRPInstance->getPar("InFileNameNeutralNeutralPart");
	  break;
	case EXPCFSH_SYSTEM_ZP:
	  fName += sRPInstance->getPar("InFileNameNeutralPositivePart");
	  break;
	case EXPCFSH_SYSTEM_ZM:
	  fName += sRPInstance->getPar("InFileNameNeutralNegativePart");
	  break;
	case EXPCFSH_SYSTEM_PZ:
	  fName += sRPInstance->getPar("InFileNamePositiveNeutralPart");
	  break;
	case EXPCFSH_SYSTEM_MZ:
	  fName += sRPInstance->getPar("InFileNameNegativeNeutralPart");
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

STR ExpCFSH::MakePurityFileName(int aChargeType)
{
  STR fName;
  try
    {
      fName = sRPInstance->getPar("InPurityFileNamePrefix");
      switch (aChargeType)
	{
	case EXPCFSH_SYSTEM_PP:
	  fName += sRPInstance->getPar("InPurityFileNamePositivePositivePart");
	  break;
	case EXPCFSH_SYSTEM_MM:
	  fName += sRPInstance->getPar("InPurityFileNameNegativeNegativePart");
	  break;
	case EXPCFSH_SYSTEM_PM:
	  fName += sRPInstance->getPar("InPurityFileNamePositiveNegativePart");
	  break;
	case EXPCFSH_SYSTEM_MP:
	  fName += sRPInstance->getPar("InPurityFileNameNegativePositivePart");
	  break;
	case EXPCFSH_SYSTEM_ZZ:
	  fName += sRPInstance->getPar("InPurityFileNameNeutralNeutralPart");
	  break;
	case EXPCFSH_SYSTEM_ZP:
	  fName += sRPInstance->getPar("InPurityFileNameNeutralPositivePart");
	  break;
	case EXPCFSH_SYSTEM_ZM:
	  fName += sRPInstance->getPar("InPurityFileNameNeutralNegativePart");
	  break;
	case EXPCFSH_SYSTEM_PZ:
	  fName += sRPInstance->getPar("InPurityFileNamePositiveNeutralPart");
	  break;
	case EXPCFSH_SYSTEM_MZ:
	  fName += sRPInstance->getPar("InPurityFileNameNegativeNeutralPart");
	  break;
	}
      fName += sRPInstance->getPar("InPurityFileNameSuffix");
      fName += '\0';
      PRINT_DEBUG_2("Got purity filename " << fName.Data());
    }
  catch (STR e)
    {
      PRINT_MESSAGE("ExpCf <MakePurityFileName>: Purity filename parameters not found " << e);
      PRINT_MESSAGE("Using regular input filename instead");
      return MakeFileName(aChargeType);
    }
  return fName;
}
void           
ExpCFSH::GenerateParameterStub()
{
  fstream *os;
  os = new fstream(sRPFileName.Data(),ios::ate|ios::out|ios::in);

  (*os) 
    << "# The parameters for the experimental correlation functions" << endl
    << "# in spherical harmonics" << endl
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
    << "# InPurityFileNamePositiveNegativePart	= " << endl
    << "# InPurityFileNameNegativePositivePart	= " << endl
    << "# InPurityFileNameNegativeNegativePart	= " << endl
    << "# InPurityFileNameNeutralNeutralPart	= " << endl
    << "# InPurityFileNameNeutralPositivePart	= " << endl
    << "# InPurityFileNameNeutralNegativePart	= " << endl
    << "# InPurityFileNamePositiveNeutralPart	= " << endl
    << "# InPurityFileNameNegativeNeutralPart	= " << endl
    << endl
    << "# The type of the purity histogram -" << endl
    << "#  - single - one histogram for both signs of the projection" << endl
    << "#  - double - two histograms, one for positive" << endl
    << "#             one for negative part." << endl
    << "#             In this case the sign part of the name is ignored" << endl
    << "PurityHistType = single" << endl
    << endl
    << "# The section with the name of the purity histograms" << endl
    << "# The purity histogram name consists of four parts" << endl
    << "#  - the prefix (identical for all)" << endl
    << "#  - the system name" << endl
    << "#  - the sign of the projection" << endl
    << "#  - the suffix (identical for all)" << endl
    << "# in an arbitrary order" << endl
    << endl
    << "PurityPartOrder = Prefix,System,Proj,Suffix" << endl
    << endl
    << "PurityPrefixPart = " << endl
    << "PuritySuffixPart = " << endl
    << endl
    << "PurityPositivePositivePart	= " << endl
    << "PurityPositiveNegativePart	= " << endl
    << "PurityNegativePositivePart	= " << endl
    << "PurityNegativeNegativePart	= " << endl
    << "PurityNeutralNeutralPart	= " << endl
    << "PurityNeutralPositivePart	= " << endl
    << "PurityNeutralNegativePart	= " << endl
    << "PurityPositiveNeutralPart	= " << endl
    << "PurityNegativeNeutralPart	= " << endl
    << endl
    << "PurityPositivePart	= P" << endl
    << "PurityNegativePart	= N" << endl
    << endl
    << "# the Proj Parts" << endl
    << "PurityOutPart	= " << endl
    << "PuritySidePart	= " << endl
    << "PurityLongPart	= " << endl
    << "" << endl
    << "# The variable in which the purity histogram was made" << endl
    << "# possible values are:" << endl
    << "# KStar, 2KStar, Qinv" << endl
    << "PurityHistogramVariable	= KStar" << endl
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
    << "InFileNamePositiveNegativePart = " << endl
    << "InFileNameNegativePositivePart = " << endl
    << "InFileNameNegativeNegativePart = " << endl
    << "InFileNameNeutralNeutralPart = " << endl
    << "InFileNameNeutralPositivePart = " << endl
    << "InFileNameNeutralNegativePart = " << endl
    << "InFileNamePositiveNeutralPart = " << endl
    << "InFileNameNegativeNeutralPart = " << endl
    << endl
    << "# The section describing the histogram names in the input file" << endl
    << "# The name consists of five sections:" << endl
    << "#	General - genereal part of the name, the same for all histograms" << endl
    << "#	Type	- type of the histogram (correlation function) " << endl
    << "#	Proj	- the projection (l=0,m=0 or l=1,m=1)" << endl
    << "#	System	- the system name (e.g. PimKp, PP, etc.)" << endl
    << "# First we say what is the order in which the above parts go" << endl
    << "PartOrder	= Type,Proj,System,General" << endl
    << endl
    << "#Next we give the values for each of the parts" << endl
    << "# First the Type parts" << endl
    << "CFPart	= Cfn" << endl
    << "CovPart	= CovCfc" << endl
    << endl
    << "# the Proj Parts" << endl
    << "ReL0M0Part	= ReYlm00" << endl
    << "ReL1M1Part	= ReYlm11" << endl
    << endl
    << "# the System Part" << endl
    << "PositivePositivePart	= " << endl
    << "PositiveNegativePart	= " << endl
    << "NegativePositivePart	= " << endl
    << "NegativeNegativePart	= " << endl
    << "NeutralNeutralPart	= " << endl
    << "NeutralPositivePart	= " << endl
    << "NeutralNegativePart	= " << endl
    << "PositiveNeutralPart	= " << endl
    << "NegativeNeutralPart	= " << endl
    << endl
    << "# the General Part" << endl
    << "GeneralPart	= " << endl
    << endl
    << "# The variable in which the data histograms were made" << endl
    << "# possible values are:" << endl
    << "# KStar, 2KStar, Qinv" << endl
    << "DataHistogramVariable	= KStar" << endl << endl
    << "# Turn this on for binominal error calculation" << endl
    << "# (usefull for theoretically generated CFs)" << endl
    << "BinominalErrors        = 0" << endl << endl
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

void ExpCFSH::ReadParameters()
{
  STR fVar;

  try 
    {
      // For now - only the non-id 3D case is supported
      // First read the order
      char *fOrder = strdup(sRPInstance->getPar("PartOrder").Data());
      PRINT_DEBUG_2("Got order: " << mOrder);
    
      char *next = strstr(fOrder,",");
      char *curr = fOrder;
      int count=0;
      // First clear the order table
      for (int yyy=0; yyy<5; yyy++)
	mOrder[yyy] = 0;

      while (curr != next)
	{
	  if (!strncmp(curr,"Type",4))
	    mOrder[count] = EXPCFSH_PART_TYPE;
	  else if (!strncmp(curr,"Proj",4))
	    mOrder[count] = EXPCFSH_PART_PROJ;
	  else if (!strncmp(curr,"System",6))
	    mOrder[count] = EXPCFSH_PART_SYSTEM;
	  else if (!strncmp(curr,"General",7))
	    mOrder[count] = EXPCFSH_PART_GENERAL;
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
      PRINT_DEBUG_2("Dishing out the order: ");
      for (int yyy=0; yyy<5; yyy++)
	PRINT_DEBUG_2(mOrder[yyy] << " ");
      
      mTypes[0] = sRPInstance->getPar("CFPart");
      mTypes[1] = sRPInstance->getPar("CovPart");
      mProjs[0] = sRPInstance->getPar("ReL0M0Part");
      mProjs[1] = sRPInstance->getPar("ReL1M1Part");
      mProjs[2] = "";
      mSystems[EXPCFSH_SYSTEM_PP] = sRPInstance->getPar("PositivePositivePart");
      mSystems[EXPCFSH_SYSTEM_PM] = sRPInstance->getPar("PositiveNegativePart");
      mSystems[EXPCFSH_SYSTEM_MP] = sRPInstance->getPar("NegativePositivePart");
      mSystems[EXPCFSH_SYSTEM_MM] = sRPInstance->getPar("NegativeNegativePart");
      mSystems[EXPCFSH_SYSTEM_ZZ] = sRPInstance->getPar("NeutralNeutralPart");
      mSystems[EXPCFSH_SYSTEM_ZP] = sRPInstance->getPar("NeutralPositivePart");
      mSystems[EXPCFSH_SYSTEM_ZM] = sRPInstance->getPar("NeutralNegativePart");
      mSystems[EXPCFSH_SYSTEM_PZ] = sRPInstance->getPar("PositiveNeutralPart");
      mSystems[EXPCFSH_SYSTEM_MZ] = sRPInstance->getPar("NegativeNeutralPart");
      mGeneral = sRPInstance->getPar("GeneralPart");
      mPurCorr = atof((sRPInstance->getPar("PurityCorrection")).Data())/100.0;
      mMomRes  = atof((sRPInstance->getPar("MomentumResolutionCorrection")).Data())/100.0;

//       if (mPurCorr>0.0001)
// 	{
// 	  // First clear the purity order
// 	  for (int yyy=0; yyy<5; yyy++)
// 	    mPOrder[yyy]= 0;

// 	  fOrder = strdup(sRPInstance->getPar("PurityPartOrder").Data());
// 	  PRINT_DEBUG_2("Got order: " << fOrder);
	  
// 	  if (sRPInstance->getPar("PurityHistType") == "single")
// 	    mPType = EXPCFSH_PURITY_TYPE_SINGLE;
// 	  else if (sRPInstance->getPar("PurityHistType") == "double")
// 	    mPType = EXPCFSH_PURITY_TYPE_DOUBLE;
// 	  else
// 	    {
// 	      PRINT_MESSAGE("ExpCF <readParameters>: wrong Purity histogram type");
// 	      exit(1);
// 	    }

// 	  next = strstr(fOrder,",");
// 	  curr = fOrder;
// 	  count=0;
// 	  while (curr != next)
// 	    {
// 	      if (!strncmp(curr,"Prefix",6))
// 		mPOrder[count] = EXPCFSH_PURITY_PART_PREFIX;
// 	      else if (!strncmp(curr,"System",6))
// 		mPOrder[count] = EXPCFSH_PURITY_PART_SYSTEM;
// 	      else if (!strncmp(curr,"Proj",4))
// 		mPOrder[count] = EXPCFSH_PURITY_PART_PROJ;
// 	      else if (!strncmp(curr,"Sign",4))
// 		mPOrder[count] = EXPCFSH_PURITY_PART_SIGN;
// 	      else if (!strncmp(curr,"Suffix",6))
// 		mPOrder[count] = EXPCFSH_PURITY_PART_SUFFIX;
// 	      else
// 		{
// 		  PRINT_MESSAGE("ExpCF <readParameters>: unknown purity order part: " << curr);
// 		  exit(1);
// 		}
// 	      if (next)
// 		{
// 		  curr = next+1;
// 		  next = strstr(curr,",");
// 		}
// 	      else
// 		curr = next;
// 	      count++;
// 	    }
// 	  PRINT_DEBUG_2("Dishing out the purity order: ");
// 	  for (int yyy=0; yyy<5; yyy++)
// 	    PRINT_DEBUG_2(mPOrder[yyy] << " ");
      
// 	  mPProjs[EXPCFSH_PROJ_OUT] = sRPInstance->getPar("PurityOutPart");
// 	  mPProjs[EXPCFSH_PROJ_SIDE] = sRPInstance->getPar("PuritySidePart");
// 	  mPProjs[EXPCFSH_PROJ_LONG] = sRPInstance->getPar("PurityLongPart");
// 	  mPSigns[EXPCFSH_SIGN_POS] = sRPInstance->getPar("PurityPositivePart");
// 	  mPSigns[EXPCFSH_SIGN_NEG] = sRPInstance->getPar("PurityNegativePart");
// 	  mPSystems[EXPCFSH_SYSTEM_PP] = sRPInstance->getPar("PurityPositivePositivePart");
// 	  mPSystems[EXPCFSH_SYSTEM_PM] = sRPInstance->getPar("PurityPositiveNegativePart");
// 	  mPSystems[EXPCFSH_SYSTEM_MP] = sRPInstance->getPar("PurityNegativePositivePart");
// 	  mPSystems[EXPCFSH_SYSTEM_MM] = sRPInstance->getPar("PurityNegativeNegativePart");
// 	  mPSystems[EXPCFSH_SYSTEM_ZZ] = sRPInstance->getPar("PurityNeutralNeutralPart");
// 	  mPSystems[EXPCFSH_SYSTEM_ZP] = sRPInstance->getPar("PurityNeutralPositivePart");
// 	  mPSystems[EXPCFSH_SYSTEM_ZM] = sRPInstance->getPar("PurityNeutralNegativePart");
// 	  mPSystems[EXPCFSH_SYSTEM_PZ] = sRPInstance->getPar("PurityPositiveNeutralPart");
// 	  mPSystems[EXPCFSH_SYSTEM_MZ] = sRPInstance->getPar("PurityNegativeNeutralPart");
// 	  mPPrefix = sRPInstance->getPar("PurityPrefixPart");
// 	  mPSuffix = sRPInstance->getPar("PuritySuffixPart");
// 	}
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
      
      fVar = sRPInstance->getPar("BinominalErrors");
      if ((fVar == "on") || (fVar == "On") || (fVar == "1"))
	mBinominalErrors = 1;
      else
	mBinominalErrors = 0;

//       fVar = sRPInstance->getPar("ReduceEdgeEffects");
//       if ((fVar == "on") || (fVar == "On") || (fVar == "1"))
// 	mCFReduceEdgeEffects = 1;
//       else
// 	mCFReduceEdgeEffects = 0;
    }
  
  catch (STR e)
    {
      PRINT_MESSAGE("Error reading parameters in ExpCF <readParameters>: " << e);
    }
}

void ExpCFSH::SetParameters(int aNFitBin, double aLowKStarFit, double aHighKStarFit,
			       double aLowKStarNorm, double aHighKStarNorm)
{
  if(mContent){
    PRINT_MESSAGE("parameters can only be set once");
  }
  else{    
    mNFitBin=aNFitBin;
    //    int tmb = (int) (aHighKStarNorm/((aHighKStarFit - aLowKStarFit)/aNFitBin));
    mContent = new double[2*mNFitBin];
    mErr2 = new double[2*mNFitBin];
    mCovariance = new double[mNFitBin];
    mLowKStarFit=aLowKStarFit;
    mHighKStarFit=aHighKStarFit;
    mLowKStarNorm=aLowKStarNorm;
    mHighKStarNorm=aHighKStarNorm;
    //    mKStarSigned = new double[mNFitBin];
//     double tStep = (mHighKStarFit-mLowKStarFit)/mNFitBin*2;
//     for(int ti=0; ti<mNFitBin/2; ti++){
//       mKStarSigned[ti] = -1.* mHighKStarFit + (ti+0.5)*tStep;
//       mKStarSigned[ti+aNFitBin/2] = mLowKStarFit + (ti+0.5)*tStep;
//     }
    PRINT_DEBUG_1( "Set parameters. Fit: " << mNFitBin << " , " <<  mLowKStarFit
		 << " < k* < " << mHighKStarFit << " | Norm: " << mLowKStarNorm
		 << " < k* < " << mHighKStarNorm);
  }
}

double ExpCFSH::GetLowKStarFit()
{
  return mLowKStarFit;
}

double ExpCFSH::GetHighKSstarFit()
{
  return mHighKStarFit;
}

double ExpCFSH::GetLowKStarNorm()
{
  return mLowKStarNorm;
}

double ExpCFSH::GetHighKStarNorm()
{
  return mHighKStarNorm;
}

double* ExpCFSH::GetCovariance()
{
  return mCovariance;
}

const double* ExpCFSH::GetCovariance() const
{
  return mCovariance;
}

