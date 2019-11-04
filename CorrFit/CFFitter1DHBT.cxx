#include "CFFitter1DHBT.h"
#include "TH2D.h"
#include <TRandom.h>
#include <TMath.h>
#include "fstream"
#include "ExpCFNonId.h"

extern TFile   *sOutFile;
extern ReadPar *sRPInstance;
extern STR      sRPFileName;

using namespace std;

CFFitter1DHBT::CFFitter1DHBT() : CFFitter()
{
  ReadParameters();
  
  mPType = mPType % 10;
  switch (mPType){
  case 1:
    mExpCF       = new ExpCF1DHBT(EXPCF_SYSTEM_PP, sOutFile);
    break;
  case 2:
    mExpCF       = new ExpCF1DHBT(EXPCF_SYSTEM_MM, sOutFile);
    break;
  case 5:
    mExpCF       = new ExpCF1DHBT(EXPCF_SYSTEM_ZZ, sOutFile);
    break;
  default:
    PRINT_MESSAGE("Icorrect pair type charge " << mPType << " . Aborting");
    exit(1);
  }
  mCalcCF      = new CalcCF1DHBT(((ExpCF1DHBT *) mExpCF)->GetLowKStarFit(),
				 ((ExpCF1DHBT *) mExpCF)->GetHighKSstarFit(),
				 ((ExpCF1DHBT *) mExpCF)->GetNFitBin(),
				 ((ExpCF1DHBT *) mExpCF)->GetLowKStarNorm(),
				 ((ExpCF1DHBT *) mExpCF)->GetHighKStarNorm());
}

CFFitter1DHBT::~CFFitter1DHBT()
{
  if (mPairManager)
    delete mPairManager;
  if (mCalcCF)
    delete mCalcCF;
  if (mExpCF)
    delete mExpCF;
  if (mChi2Map)
    delete mChi2Map;
  if (mChi2DetMap)
    delete mChi2DetMap;
}

void 
CFFitter1DHBT::Initialize() 
{
  CFFitter::Initialize();

  Int_t tPCount;
  Pair *tPair;
  Int_t tNorm=0;
  
  for (int tPCount=0; tPCount<mPairManager->GetNormPairCount(); tPCount++) {
    tPair = mPairManager->GetNormPair(tPCount);
    tNorm++;
  }
  ((CalcCF1DHBT *) mCalcCF)->SetNormBins(tNorm);
}

void 
CFFitter1DHBT::Fit()
{
  CFFitter::Fit();
}

void 
CFFitter1DHBT::Write()
{
  TH2D *tSIm;

  TH1D *tROutS;
  TH1D *tRSideS;
  TH1D *tRLongS;
  TH1D *tROut;
  TH1D *tRSide;
  TH1D *tRLong;
  TH1D *tRTime;
  TH1D *tRStar;
  TH1D *tRInv;

  TGraphErrors *tG;

  mChi2Map->WriteHisto();
  mChi2DetMap->WriteHisto();
  
  mSourceModel->SetModelParameters(mBestPar);
  mCalcCF->Generate();
  tG = mCalcCF->GetGraph(mBestPur);
  tG->Write();
  
  // Source imaging histograms
  tSIm = new TH2D("tSIm","tSIm",50,-10.0,10.0,50,-10.0,10.0);
  tROutS  = new TH1D("tRoutS","r*_{out} distribution",100,-50.0,50.0);
  tRSideS = new TH1D("tRsideS","r*_{side} distribution",100,-50.0,50.0);
  tRLongS = new TH1D("tRlongS","r*_{long} distribution",100,-50.0,50.0);
  tROut   = new TH1D("tRout"  ,"r_{out} distribution",100,-50.0,50.0);
  tRSide  = new TH1D("tRside" ,"r_{side} distribution",100,-50.0,50.0);
  tRLong  = new TH1D("tRlong" ,"r_{long} distribution",100,-50.0,50.0);
  tRTime  = new TH1D("tRtime" ,"r_{time} distribution",100,-50.0,50.0);
  tRStar  = new TH1D("tRstar" ,"r* distribution",100,-50.0,50.0);
  tRInv   = new TH1D("tRInv"  ,"r* distribution scaled by 1/r^2",100,-50.0,50.0);

  tSIm = new TH2D("tSIm","tSIm",50,-10.0,10.0,50,-10.0,10.0);
  Pair *tPair;
  for(int ti=0; ti<mPairManager->GetFitPairCount(); ti++){
    tPair = mPairManager->GetFitPair(ti);
    tPair->SetPosition();
    tSIm->Fill(tPair->x2().x,tPair->x2().y);

    double tPx = tPair->p1().x+tPair->p2().x;
    double tPy = tPair->p1().y+tPair->p2().y;
    double tPz = tPair->p1().z+tPair->p2().z;
    double tE  = tPair->p1().t+tPair->p2().t;
    double tPt = tPx*tPx + tPy*tPy;
    double tMt = tE*tE - tPz*tPz;
    double tM  = sqrt(tMt - tPt);
    tMt = sqrt(tMt);
    tPt = sqrt(tPt);
    
    double tDX = tPair->x1().x-tPair->x2().x;
    double tDY = tPair->x1().y-tPair->x2().y;
    double mRLong     = tPair->x1().z-tPair->x2().z;
    double mDTime     = tPair->x1().t-tPair->x2().t;
  
    double mRTrans = tDX>0.? sqrt(tDX*tDX+tDY*tDY) : -1.*sqrt(tDX*tDX+tDY*tDY);
    double mROut = (tDX*tPx + tDY*tPy)/tPt;
    double mRSide = (-tDX*tPy + tDY*tPx)/tPt;
    double mRSidePairCMS = mRSide;
    
    // Lab -> LCMS -> PRF method
    double tBeta = tPz/tE;
    double tGamma = tE/tMt;
    double mRLongPairCMS = tGamma*(mRLong - tBeta* mDTime);
    double mDTimePairLCMS = tGamma*(mDTime - tBeta* mRLong);
    tBeta = tPt/tMt;
    tGamma = tMt/tM;
    double mROutPairCMS = tGamma*(mROut - tBeta* mDTimePairLCMS);
    double mDTimePairCMS = tGamma*(mDTimePairLCMS - tBeta* mROut);
    double mRStar = sqrt(mROutPairCMS*mROutPairCMS + mRSidePairCMS*mRSidePairCMS +
		  mRLongPairCMS*mRLongPairCMS);
    
    tROutS->Fill(mROutPairCMS);
    tRSideS->Fill(mRSidePairCMS);
    tRLongS->Fill(mRLongPairCMS);
    tROut->Fill(mROut);
    tRSide->Fill(mRSide);
    tRLong->Fill(mRLong);
    tRTime->Fill(mDTime);
    if (mRStar != 0.0)
      tRInv->Fill(mROutPairCMS > 0.0 ? mRStar : -1.0 * mRStar, 1.0/(mRStar*mRStar));
    tRStar->Fill(mROutPairCMS > 0.0 ? mRStar : -1.0*mRStar);
    
//     if ((tPair->x1().x > 0.118557) && (tPair->x1().x < 0.118559) &&
// 	(tPair->x2().x > 0.160845) && (tPair->x2().x < 0.160847)) 
//     {
//       PRINT_MESSAGE("Got Px1: "   << tPair->x1().x);
//       PRINT_MESSAGE("Got Px2: "   << tPair->x2().x);
//       PRINT_MESSAGE("Got RStar "  << mRStar);
//       PRINT_MESSAGE("Got ROutS "  << mROutPairCMS);
//       PRINT_MESSAGE("Got RSideS " << mRSidePairCMS);
//       PRINT_MESSAGE("Got RLongS " << mRLongPairCMS);
//       PRINT_MESSAGE("Got RTimeS " << mDTimePairCMS);
//       PRINT_MESSAGE("Got ROut  "  << tPair->p2().x);
//       PRINT_MESSAGE("Got RSide "  << tPair->p2().y);
//       PRINT_MESSAGE("Got RLong "  << tPair->p2().z);
//       PRINT_MESSAGE("Got RTime "  << tPair->p2().t);
//     }
    
}
  tSIm->Write();

  tROutS->Write();
  tRSideS->Write();
  tRLongS->Write();
  tROut->Write();
  tRSide->Write();
  tRLong->Write();
  tRTime->Write();
  tRStar->Write();
  tRInv->Write();
  
  CFFitter::Write();
}

double 
CFFitter1DHBT::GetChi2()
{
  double tChi2=0.;
  const double* tContent = mCalcCF->GetContent();
  const double* tErr2 = mCalcCF->GetError2();
  const double* tContentExp = mExpCF->GetContent();
  const double* tErr2Exp = mExpCF->GetError2();
  double tVal;
  for(int ti=0; ti<((ExpCF1DHBT *) mExpCF)->GetNFitBin(); ti++){
    tVal = tContentExp[ti]-tContent[ti];
    //    tVal = mContent[ti]-((tContent[ti]-1.)/aPurity+1.);
    //    tChi2 += (tVal*tVal/(tErr2Exp[ti]+tErr2[ti]/aPurity/aPurity));   
    tChi2 += (tVal*tVal/(tErr2Exp[ti]+tErr2[ti]));   
    PRINT_DEBUG_3("Chi bin " << 
		  ti << " " <<
		  tContentExp[ti] << " " << 
		  tContent[ti] << " " << 
		  tErr2Exp[ti] << " " << 
		  tErr2[ti] << " " << 
		  tVal << "    " <<
		  (tVal*tVal/(tErr2Exp[ti]+tErr2[ti])) << "    " <<
		  tChi2);
  }
  return tChi2;
}

void     
CFFitter1DHBT::InitializePairManager()
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
  
  //  GeneratePairs(mNFPBuf, mCalcCF->GetNFitBins());

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

  if (mGenerate)
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

// This function is used to generate additional pairs
// in bins where there is only a few of them. 
// This may improve the fit stability
void 
CFFitter1DHBT::GeneratePairs(int *aPairsPerBin, int aNBins)
{
  //  return ;
  Double_t tBinMin, tBinMax, tBinMaxValue, farg, ftest, fval;
  Double_t fout, flong, fside;
  Double_t px1, py1, pz1, e1, px2, py2, pz2, e2;

  Int_t tPairCount = 1;
  Int_t tInitPair = mPairManager->GetFitPairCount();
  Int_t tPairGenerated = 1;
  Int_t tMaxBinCount=mMinPairsPerBin;

  Float_t tPairMom[8];
  Pair *tRealPair;

  TRandom rand;
  Pair *tPair;

  rand.SetSeed(mRandomSeed+3);
  
  Int_t *tGenNFPBuf;
  tGenNFPBuf = (Int_t *) malloc(sizeof(Int_t) * aNBins);
  for (int iter=0; iter<aNBins; iter++) {
    tGenNFPBuf[iter] = 0;
    if (aPairsPerBin[iter] > tMaxBinCount) tMaxBinCount = aPairsPerBin[iter];
  }
  
  TF1 *tFun = ((ExpCF1DHBT *) mExpCF)->mHDenFit;  
  for (int iter=0; iter<4; iter++) {
    PRINT_DEBUG_3("Par " << iter << " " << tFun->GetParameter(iter));
  }
  Double_t tRangeMin = ((ExpCF1DHBT *) mExpCF)->GetLowKStarFit();
  Double_t tRangeMax = ((ExpCF1DHBT *) mExpCF)->GetHighKSstarFit();
  Double_t tRangeStep = (tRangeMax-tRangeMin)/aNBins;
  
  // In each k* bin do the following:
  for (int nbin=0; nbin<aNBins; nbin++){
    tBinMin = tRangeMin + nbin*tRangeStep;
    tBinMax = tRangeMin + (nbin+1)*tRangeStep;
    tBinMaxValue = tFun->Eval(tBinMax);

    PRINT_DEBUG_2("Generating pairs in bin " << nbin << " " << tBinMin << " " << tBinMax << " " << tBinMaxValue << " " << aPairsPerBin[nbin]);
    for (int npair=0; npair<tMaxBinCount-aPairsPerBin[nbin]; npair++) {
      //   Generate k* with Monte-Carlo in the range of the bin
      do{
	farg = rand.Rndm()*tRangeStep+tBinMin;
	ftest = rand.Rndm()*tBinMaxValue;
	fval = tFun->Eval(farg);
      }
      while (ftest > fval);
      
      // Correctly generate k*Out, k*Side and k*Long
      // First generate three components
      double tComp1 = rand.Rndm()*farg;
      double tComp2 = rand.Rndm()*TMath::Sqrt(farg*farg - tComp1*tComp1);
      double tComp3 = TMath::Sqrt(farg*farg - tComp1*tComp1 - tComp2*tComp2);

      double tCase = rand.Rndm()*6.0;
      if (tCase < 1.0) {
	fout = tComp1;
	fside = tComp2;
	flong = tComp3;
      }
      else if (tCase < 2.0) {
	fout = tComp1;
	fside = tComp3;
	flong = tComp2;
      }
      else if (tCase < 3.0) {
	fout = tComp2;
	fside = tComp1;
	flong = tComp3;
      }
      else if (tCase < 4.0) {
	fout = tComp2;
	fside = tComp3;
	flong = tComp1;
      }
      else if (tCase < 5.0) {
	fout = tComp3;
	fside = tComp1;
	flong = tComp2;
      }
      else {
	fout = tComp3;
	fside = tComp2;
	flong = tComp1;
      }
      
      // Set the proper fout sign
      fout = (rand.Rndm() > 0.5) ? fout : -fout;
      
      // Generate k*Long sign
      flong = (rand.Rndm() > 0.5) ? flong : -flong;

      // Generating k*Side sign
      fside = (rand.Rndm() > 0.5) ? fside : -fside;
      
      //   Get px1, py1, and pz1 from an existing pair
      tPair = mPairManager->GetFitPair(tPairCount);
      if (tPairCount++ > tInitPair) tPairCount = 1;

      //   Calculate px2, py2, pz2
      Double_t tPx = tPair->p1().x + tPair->p2().x;
      Double_t tPy = tPair->p1().y + tPair->p2().y;
      Double_t tPz = tPair->p1().z + tPair->p2().z;
      Double_t tPe = tPair->p1().t + tPair->p2().t;
      Double_t tm1sq = tPair->p1().t*tPair->p1().t - tPair->p1().x*tPair->p1().x - tPair->p1().y*tPair->p1().y - tPair->p1().z*tPair->p1().z;

      // Rotating to Py=0
      tPx = TMath::Hypot(tPx,tPy);
      tPy = 0;
      
      // Other variables
      Double_t tMt = tPe*tPe - tPz*tPz;
      Double_t tPt = tPx*tPx;
      Double_t tMass = TMath::Sqrt(tMt - tPt);
      tMt = TMath::Sqrt(tMt);
      tPt = TMath::Sqrt(tPt);
      Double_t beta = tPz/tPe;
      Double_t gamma = tPe/tMt;

      // Now we get the variables
      py1 = fside;
      Double_t tG = flong/gamma;
      Double_t tF = (tPt/tMt)*gamma*(1-beta*beta);
      Double_t tD = (tMass/tMt)*fout-(tPt/tMt)*beta*flong;
      Double_t tP = tm1sq + tD*tD + fside*fside + tG*tG;
      Double_t tK = tF*tF + beta*beta - 1;
      Double_t tL = 2*tD*tF + 2*tG*beta;
      
      e1 = (- tL - TMath::Sqrt(tL*tL - 4*tK*tP))/(2*tK);
      pz1 = tG + beta*e1;
      px1 = tD + tF*e1;

      // And now we get the second particle variables
      px2 = tPx - px1;
      py2 = tPy - py1;
      pz2 = tPz - pz1;
      e2  = tPe - e1;
      
      //      PRINT_DEBUG_3("Generated pair: " << farg << " " << fout << " " << flong << " " << fside << " " << TMath::Sqrt(e2*e2-px2*px2-py2*py2-pz2*pz2));
      
      // Rotating the particles by a random angle
      Double_t tRot = rand.Rndm()*2*TMath::Pi();
      Double_t pt = TMath::Hypot(px1,py1);
      Double_t phi = TMath::ATan2(py1,px1);
      phi += tRot;
      px1 = pt*TMath::Cos(phi);
      py1 = pt*TMath::Sin(phi);

      pt = TMath::Hypot(px2,py2);
      phi = TMath::ATan2(py2,px2);
      phi += tRot;
      px2 = pt*TMath::Cos(phi);
      py2 = pt*TMath::Sin(phi);

      // Store the generated pair
      tPairMom[0] = px1;
      tPairMom[1] = py1;
      tPairMom[2] = pz1;
      tPairMom[3] = e1;
      tPairMom[4] = px2;
      tPairMom[5] = py2;
      tPairMom[6] = pz2;
      tPairMom[7] = e2;
      
      tRealPair = new Pair();
      tRealPair->SetMomentum(tPairMom);
      tRealPair->SetPairNum(10000000+tPairGenerated++);
      
//       PRINT_DEBUG_3("k* gen calc " << farg  << " " << tRealPair->GetKStar());
//       PRINT_DEBUG_3("k* gen calc " << fout  << " " << tRealPair->GetKStarTimeOutSign());
//       PRINT_DEBUG_3("k* gen calc " << fside << " " << tRealPair->GetKStarTimeSideSign());
//       PRINT_DEBUG_3("k* gen calc " << flong << " " << tRealPair->GetKStarTimeLongSign());
      
      if (mCalcCF->InFitRange(tRealPair)) {
	tRealPair->SetBin(mCalcCF->FitRangeBin(tRealPair));
	tGenNFPBuf[mCalcCF->FitRangeBin(tRealPair)]++;
	mPairManager->StoreFitPair(tRealPair);
      }
      
    }
  }
  
  for (int ti=0; ti<mCalcCF->GetNFitBins(); ti++) {
    mCalcCF->SetFitBinCount(ti, aPairsPerBin[ti] + tGenNFPBuf[ti]);
    PRINT_DEBUG_2("Bin number " << ti << " generated additional " << tGenNFPBuf[ti] << " pairs");
  }
}


double  
CFFitter1DHBT::GetChi2(double aPurity)
{
  double tChi2=0.;
  const double* tContent = mCalcCF->GetContent();
  const double* tErr2 = mCalcCF->GetError2();
  const double* tContentExp = mExpCF->GetContent();
  const double* tErr2Exp = mExpCF->GetError2();
  double tVal;
  for(int ti=0; ti<((ExpCF1DHBT *) mExpCF)->GetNFitBin(); ti++){
    // tVal = tContentExp[ti]-tContent[ti];
    tVal = tContentExp[ti]-((tContent[ti]-1.)*aPurity+1.);
    tChi2 += (tVal*tVal/(tErr2Exp[ti]+tErr2[ti]*aPurity*aPurity));   
    // tChi2 += (tVal*tVal/(tErr2Exp[ti]+tErr2[ti]));   
    PRINT_DEBUG_3("Chi bin " << 
		  ti << " " << 
		  tContentExp[ti] << " " << 
		  tContent[ti] << " " << 
		  ((tContent[ti]-1.0)*aPurity + 1.0) << " " << 
		  aPurity << " " << 
		  tErr2Exp[ti] << " " << 
		  tErr2[ti] << " " << 
		  tVal << "    " <<
		  (tVal*tVal/(tErr2Exp[ti]+tErr2[ti])) << "    " <<
		  tChi2);
  }
  return tChi2;
}

void 
CFFitter1DHBT::GenerateParameterStub()
{
  CFFitter::GenerateParameterStub();

  fstream *os;
  os = new fstream(sRPFileName.Data(),ios::ate|ios::out|ios::in);

  (*os) 
    << "# 1DHBT Fitter specific parameters" << endl 
    << "PairType = " << endl << endl;
  (*os) 
    << "# Minimum number of pairs in bin" << endl 
    << "MinPairsPerBin = 10000" << endl << endl;
  (*os) 
    << "# Set to 1 to randomly generate additional pairs (only for testing!)" << endl 
    << "GeneratePairs = 0" << endl << endl;

  os->close();

  CalcCF1DHBT::GenerateParameterStub();
  ExpCF1DHBT::GenerateParameterStub();
  
}

void 
CFFitter1DHBT::ReadParameters()
{
  CFFitter::ReadParameters();
  
  try {
    mPType = atoi((sRPInstance->getPar("PairType")).Data());
    mMinPairsPerBin = atoi((sRPInstance->getPar("MinPairsPerBin")).Data());
  }
  catch (STR e){
    PRINT_MESSAGE("CFFitter1DHBT::ReadParameters - Error reading parameter: " << e.Data());
    exit(0);
  }
  try { 
    mRandomSeed = atoi((sRPInstance->getPar("RandomSeed")).Data());
  }
  catch (STR e){
    PRINT_MESSAGE("CFFitterNonId::ReadParameters - Parameter not found: " << e.Data() << " Using default.");
    mRandomSeed = 21341;
    
  }
  try { 
    mGenerate = atoi((sRPInstance->getPar("GeneratePairs")).Data());
  }
  catch (STR e){
    PRINT_MESSAGE("CFFitterNonId::ReadParameters - Parameter not found: " << e.Data() << " Using default.");
    mGenerate = 0;
    
  }

}

