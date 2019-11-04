#include "CFFitterQuad.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TF1.h"
#include "fstream"
#include "SourceModelGausROut.h"

extern TFile   *sOutFile;
extern ReadPar *sRPInstance;
extern STR      sRPFileName;

using namespace std;

CFFitterQuad::CFFitterQuad() : CFFitter()
{
  int tNType;
  ReadParameters();
  
  mExpCFs  = new ExpCFQuad * [mNPTypes];
  mCalcCFs = new CalcCFQuad * [mNPTypes];
  
  tNType = 0;
  do {
    int tSys = mPTypes[tNType] % 10;
    switch (tSys){
    case 1:
      mExpCFs[tNType] = new ExpCFQuad(EXPCF_SYSTEM_PP, sOutFile);
      break;
    case 2:
      mExpCFs[tNType] = new ExpCFQuad(EXPCF_SYSTEM_MM, sOutFile);
      break;
    case 3:
      mExpCFs[tNType] = new ExpCFQuad(EXPCF_SYSTEM_PM, sOutFile);
      break;
    case 4:
      mExpCFs[tNType] = new ExpCFQuad(EXPCF_SYSTEM_MP, sOutFile);
      break;
    case 5:
      mExpCFs[tNType] = new ExpCFQuad(EXPCF_SYSTEM_ZZ, sOutFile);
      break;
    default:
      PRINT_MESSAGE("Icorrect pair type charge " << mPType << " . Aborting");
      exit(1);
    }
    mCalcCFs[tNType]  = new CalcCFQuad(((ExpCFQuad *) mExpCFs[tNType])->GetLowKStarFit(),
					((ExpCFQuad *) mExpCFs[tNType])->GetHighKSstarFit(),
					((ExpCFQuad *) mExpCFs[tNType])->GetNFitBin(),
					((ExpCFQuad *) mExpCFs[tNType])->GetLowKStarNorm(),
					((ExpCFQuad *) mExpCFs[tNType])->GetHighKStarNorm());
    mCalcCFs[tNType]->SetPairSystem((CorrFitPairSystem) mPTypes[tNType]);
  }
  while (++tNType < mNPTypes);
  mExpCF = mExpCFs[0];
  mCalcCF = mCalcCFs[0];
  if (strstr(mSourceModel->GetModelIdentifier(),"SMGausROut")) 
    {
      ((SourceModelGausROut *) mSourceModel)->mRandom->SetSeed(mRandomSeed);
    }
}

CFFitterQuad::~CFFitterQuad()
{
  if (mPairManager)
    delete mPairManager;
  for (int ti=0; ti<mNPTypes; ti++) {
    delete mCalcCFs[ti];
    delete mExpCFs[ti];
  }
  delete mExpCFs;
  delete mCalcCFs;
  if (mChi2Map)
    delete mChi2Map;
  if (mChi2DetMap)
    delete mChi2DetMap;
}

void     
CFFitterQuad::InitializePairManager()
{
  Int_t   tPCount;
  Pair   *tPair;
  Int_t   tNormN=0, tNormP = 0;
  Int_t  *tNFPBuf;
  Int_t   tPairUsed;
  UInt_t  seed1, seed2;
  mPairManagers[0] = mPairManager = new PairManager(mPairFileName.Data());
  for (int ti=1; ti<mNPTypes; ti++) {
    mPairManagers[ti] = new PairManager(mPFileNames[ti].Data());
  }

  // Save the random seed
  Pair::GetSeed(seed1, seed2);

  for (int tT=0; tT<mNPTypes; tT++) {
    tNormN = 0;
    tNormP = 0;

    tNFPBuf = (Int_t *) malloc(sizeof(Int_t) * mCalcCFs[tT]->GetNFitBins());
    for (int ti=0; ti<mCalcCFs[tT]->GetNFitBins(); ti++)
      tNFPBuf[ti] = 0;
    
    Pair::SetMomRes(mMomRes, mPTypes[tT]); 
    Pair::SetSeed(seed1, seed2);

    for (tPCount=0; tPCount<mPairManagers[tT]->GetPairCount(); tPCount++){
      tPair = mPairManagers[tT]->ReadPair();
      tPairUsed = 0;
      
      if (mCalcCFs[tT]->InFitRange(tPair)) {
	tPair->SetBin(mCalcCFs[tT]->FitRangeBin(tPair));
	if (tNFPBuf[mCalcCFs[tT]->FitRangeBin(tPair)] < maxpairperbin) {
	  tNFPBuf[mCalcCFs[tT]->FitRangeBin(tPair)]++;
	  mPairManagers[tT]->StoreFitPair(tPair);
	  tPairUsed = 1;
	}
      }
      if (mCalcCFs[tT]->InNormRange(tPair)) {
	tPair->SetBin(mCalcCFs[tT]->NormRangeBin(tPair));
      if (tPair->GetKStarTimeOutSign() > 0.0)
	{
	  if (tNormP < maxpairperbin*3) {
	    tNormP++;
	    mPairManagers[tT]->StoreNormPair(tPair);
	    tPairUsed = 1;
	  }
	}
      else
	{
	  if (tNormN < maxpairperbin*3) {
	    tNormN++;
	    mPairManagers[tT]->StoreNormPair(tPair);
	    tPairUsed = 1;
	  }
	}
// 	mPairManagers[tT]->StoreNormPair(tPair);
// 	tPairUsed = 1;
      }
      
      if (!tPairUsed) delete tPair;
    }
    
    mPairManagers[tT]->CloseFile();

    for (int ti=0; ti<mCalcCFs[tT]->GetNFitBins(); ti++) {
      mCalcCFs[tT]->SetFitBinCount(ti, tNFPBuf[ti]);
      PRINT_DEBUG_3("Fit pairs in bin " << ti << " " << tNFPBuf[ti]);
    }
    
//     tNormP = 0;
//     tNormN = 0;
    
//     for (int tPCount=0; tPCount<mPairManagers[tT]->GetNormPairCount(); tPCount++) {
//       tPair = mPairManagers[tT]->GetNormPair(tPCount);
//       if (tPair->GetKStarTimeOutSign() > 0.0)
// 	tNormP++;
//       else
// 	tNormN++;
//     }
    
    ((CalcCFQuad *) mCalcCFs[tT])->SetNormPBins(tNormP);
    ((CalcCFQuad *) mCalcCFs[tT])->SetNormNBins(tNormN);
    ((CalcCFQuad *) mCalcCFs[tT])->SetPairManager(mPairManagers[tT]);
    
    PRINT_MESSAGE("CFFitter::Initialize");
    PRINT_MESSAGE("Number of Fit Pairs read: " << mPairManagers[tT]->GetFitPairCount());
    PRINT_MESSAGE("Number of Norm Pairs read: " << mPairManagers[tT]->GetNormPairCount());
    PRINT_MESSAGE("NormN: " << tNormN);
    PRINT_MESSAGE("NormP: " << tNormP);

    mPairManager = mPairManagers[tT];
    mCalcCF      = mCalcCFs[tT];
    mExpCF       = mExpCFs[tT];
    
    GeneratePairs(tNFPBuf, mCalcCF->GetNFitBins());
    
  }

  mPairManager = mPairManagers[0];
  mCalcCF      = mCalcCFs[0];
  mExpCF       = mExpCFs[0];

}

// This function is used to generate additional pairs
// in bins where there is only a few of them. 
// This may improve the fit stability
void 
CFFitterQuad::GeneratePairs(int *aPairsPerBin, int aNBins)
{
  return;
  
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
  
  TF1 *tFunP = ((ExpCFQuad *) mExpCF)->mHDenOutPFit;  
  TF1 *tFunN = ((ExpCFQuad *) mExpCF)->mHDenOutNFit;
  for (int iter=0; iter<4; iter++) {
    PRINT_DEBUG_3("Par " << iter << " " << tFunP->GetParameter(iter) << " " << tFunN->GetParameter(iter));
  }
  Double_t tRangeMin = ((ExpCFQuad *) mExpCF)->GetLowKStarFit();
  Double_t tRangeMax = ((ExpCFQuad *) mExpCF)->GetHighKSstarFit();
  Double_t tRangeStep = (tRangeMax-tRangeMin)*2/aNBins;
  
  // In each k* bin do the following:
  for (int nbin=0; nbin<aNBins; nbin++){
    if (nbin >= aNBins/2) {
      tBinMin = tRangeMin + (nbin-aNBins/2)*tRangeStep;
      tBinMax = tRangeMin + (nbin-aNBins/2+1)*tRangeStep;
      tBinMaxValue = tFunP->Eval(tBinMax);
    }
    else {
      tBinMin = tRangeMin + (aNBins/2-nbin-1)*tRangeStep;
      tBinMax = tRangeMin + (aNBins/2-nbin)*tRangeStep;
      tBinMaxValue = tFunN->Eval(tBinMax);
    }

    PRINT_DEBUG_2("Generating pairs in bin " << nbin << " " << tBinMin << " " << tBinMax << " " << tBinMaxValue << " " << aPairsPerBin[nbin]);
    for (int npair=0; npair<tMaxBinCount-aPairsPerBin[nbin]; npair++) {
      //   Generate k* with Monte-Carlo in the range of the bin
      do{
	farg = rand.Rndm()*tRangeStep+tBinMin;
	ftest = rand.Rndm()*tBinMaxValue;
	if (nbin >= aNBins/2)
	  fval = tFunP->Eval(farg);
	else
	  fval = tFunN->Eval(farg);
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
      if (!(nbin >= aNBins/2))
	fout = -fout;
      
      // Generate k*Long sign
      flong = (rand.Rndm() > 0.5) ? flong : -flong;

      // Generating k*Side sign
      fside = (rand.Rndm() > 0.5) ? fside : -fside;

//       //   Generate k*Out in range (0,k*)
//       fout = rand.Rndm()*farg;
//       if (!(nbin >= aNBins/2))
// 	fout = -fout;
      
//       //   Generate k*Long in range (0,sqrt(k* **2-K*Out**2))
//       flong = rand.Rndm()*(TMath::Sqrt(farg*farg-fout*fout))*2 - TMath::Sqrt(farg*farg-fout*fout);

//       //   Calculate k*Side
//       fside = TMath::Sqrt(farg*farg - fout*fout - flong*flong);
//       if (rand.Rndm() > 0.5) fside *= -1.0;
      
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
      
      PRINT_DEBUG_3("Generated pair: " << farg << " " << fout << " " << flong << " " << fside << " " << TMath::Sqrt(e2*e2-px2*px2-py2*py2-pz2*pz2));
      
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

void 
CFFitterQuad::Initialize() 
{
  CFFitter::Initialize();
}

void 
CFFitterQuad::Fit()
{
  CFFitter::Fit();
}

void 
CFFitterQuad::Write()
{
  TH2D *tSIm;

  TGraphErrors *tG;
  mChi2Map->WriteHisto();
  mChi2DetMap->WriteHisto();
  
  tSIm = new TH2D("tSIm","tSIm",50,-10.0,10.0,50,-10.0,10.0);
  Pair *tPair;
  for(int ti=0; ti<mPairManager->GetFitPairCount(); ti++){
    tPair = mPairManager->GetFitPair(ti);
    tPair->SetPosition();
    tSIm->Fill(tPair->x2().x,tPair->x2().y);
  }
  tSIm->Write();
  
  for (int ti=0; ti<mNPTypes; ti++) {
    mExpCFs[ti]->Write();
  }
  
  mSourceModel->SetModelParameters(mBestPar);
  for (int ti=0; ti<mNPTypes; ti++) {
    mCalcCFs[ti]->Generate();
    tG = mCalcCFs[ti]->GetGraph(mBestPur);
    tG->Write();
  }

  const double* tContentLike   = mCalcCFs[0]->GetContent();
  const double* tContentUnlike = mCalcCFs[1]->GetContent();
  const double* tErr2Like      = mCalcCFs[0]->GetError2();
  const double* tErr2Unlike    = mCalcCFs[1]->GetError2();

  const double* tContentExpLike   = ((ExpCF *) mExpCFs[0])->GetContent();
  const double* tContentExpUnlike = ((ExpCF *) mExpCFs[1])->GetContent();
  const double* tErr2ExpLike      = mExpCFs[0]->GetError2();
  const double* tErr2ExpUnlike    = mExpCFs[1]->GetError2();

  Double_t *xy;
  Double_t *yk;
  Double_t *xye;
  Double_t *yke;

  xy = (Double_t *) malloc(sizeof(Double_t ) * ((ExpCFQuad *) mExpCFs[0])->GetNFitBin());
  yk = (Double_t *) malloc(sizeof(Double_t ) * ((ExpCFQuad *) mExpCFs[0])->GetNFitBin());
  yke = (Double_t *) malloc(sizeof(Double_t ) * ((ExpCFQuad *) mExpCFs[0])->GetNFitBin());
  
  Int_t nbin = ((ExpCFQuad *) mExpCFs[0])->GetNFitBin();

  double tCalcL, tCalcU;
  for (int iteri = 0; iteri<nbin; iteri++) {
    xy[iteri] = tG->GetX()[iteri];
    tCalcL = (tContentLike[iteri] - 1.0) * ((ExpCFQuad *) mExpCFs[0])->GetPurity(iteri) / mBestPur + 1.0;
    tCalcU = (tContentUnlike[iteri] - 1.0) * ((ExpCFQuad *) mExpCFs[1])->GetPurity(iteri) / mBestPur + 1.0;
    yk[iteri] = tCalcL/tCalcU;
    yke[iteri] = ((tErr2Like[iteri]*tCalcU*tCalcU +
		   tErr2Unlike[iteri]*tCalcL*tCalcL)/
		  (tCalcU*tCalcU*tCalcU*tCalcU));
    yke[iteri] = sqrt(yke[iteri]);
    //    yke[iteri] = yk[iteri]*(sqrt(tErr2Like[iteri]) * (((ExpCFQuad *) mExpCFs[0])->GetPurity(iteri) / mBestPur) / tCalcL +
    //			    sqrt(tErr2Unlike[iteri]) * (((ExpCFQuad *) mExpCFs[1])->GetPurity(iteri) / mBestPur) / tCalcU);
  }

  TGraphErrors *tGCLoU = new TGraphErrors(nbin, xy, yk, NULL, yke);
  tGCLoU->SetName("CalculatedLikeOverUnlike");
  tGCLoU->SetTitle("CalculatedLikeOverUnlike");

  for (int iteri = 0; iteri<nbin; iteri++) {
    yk[iteri] = tContentExpLike[iteri]/tContentExpUnlike[iteri];
    yke[iteri] = ((tErr2ExpLike[iteri]*tContentExpUnlike[iteri]*tContentExpUnlike[iteri] +
		   tErr2ExpUnlike[iteri]*tContentExpLike[iteri]*tContentExpLike[iteri])/
		  (tContentExpUnlike[iteri]*tContentExpUnlike[iteri]*tContentExpUnlike[iteri]*tContentExpUnlike[iteri]));
    yke[iteri] = sqrt(yke[iteri]);
//     yke[iteri] = yk[iteri]*(sqrt(tErr2ExpLike[iteri]) / tContentExpLike[iteri] +
// 			    sqrt(tErr2ExpUnlike[iteri]) / tContentExpUnlike[iteri]);
  }

  TGraphErrors *tGELoU = new TGraphErrors(nbin, xy, yk, NULL, yke);
  tGELoU->SetName("ExperimentalLikeOverUnlike");
  tGELoU->SetTitle("ExperimentalLikeOverUnlike");

  tGCLoU->Write();
  tGELoU->Write();

  int tDim = mSourceModel->GetNParameters();
  TVectorD *tResult = new TVectorD(tDim+1);
  for (int ti=0; ti<tDim; ti++)
    (*tResult)(ti) = mBestPar[ti];
  (*tResult)(tDim) = mBestPur;

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

double 
CFFitterQuad::GetChi2()
{
  // Do the alternate method of Chi2 calculation for multiple functions
  // Weight a contribution from each function by the inverse
  // of chi2/dof for this function
  double tChi2=0.;
  double tChi2Func;
  //  for (int tT=0; tT<mNPTypes; tT++) {
  const double* tContentLike   = mCalcCFs[0]->GetContent();
  const double* tContentUnlike = mCalcCFs[1]->GetContent();
  const double* tErr2Like      = mCalcCFs[0]->GetError2();
  const double* tErr2Unlike    = mCalcCFs[1]->GetError2();

  const double* tContentExpLike   = ((ExpCF *) mExpCFs[0])->GetContent();
  const double* tContentExpUnlike = ((ExpCF *) mExpCFs[1])->GetContent();
  const double* tErr2ExpLike      = mExpCFs[0]->GetError2();
  const double* tErr2ExpUnlike    = mExpCFs[1]->GetError2();

  double tVal, tCalcL, tCalcU, tCalcLErr, tCalcUErr;
  double tExpLErr, tExpUErr, tErrExp, tErrCalc, tErrFull;
  tChi2Func = 0.0;
  for(int ti=0; ti<((ExpCFQuad *) mExpCFs[0])->GetNFitBin(); ti++){
    tCalcL = (tContentLike[ti] - 1.0) * ((ExpCFQuad *) mExpCFs[0])->GetPurity(ti) + 1.0;
    tCalcU = (tContentUnlike[ti] - 1.0) * ((ExpCFQuad *) mExpCFs[1])->GetPurity(ti) + 1.0;
    tCalcLErr = sqrt(tErr2Like[ti]) * ((ExpCFQuad *) mExpCFs[0])->GetPurity(ti);
    tCalcUErr = sqrt(tErr2Unlike[ti]) * ((ExpCFQuad *) mExpCFs[1])->GetPurity(ti);
    tExpLErr = sqrt(tErr2ExpLike[ti]);
    tExpUErr = sqrt(tErr2ExpUnlike[ti]);
    tErrExp = ((tExpLErr*tExpLErr*tContentExpUnlike[ti]*tContentExpUnlike[ti] + 
		tExpUErr*tExpUErr*tContentExpLike[ti]*tContentExpLike[ti])/
	       (tContentExpUnlike[ti]*tContentExpUnlike[ti]*tContentExpUnlike[ti]*tContentExpUnlike[ti]));
    tErrCalc = ((tCalcLErr*tCalcLErr*tCalcU*tCalcU +
		 tCalcUErr*tCalcUErr*tCalcL*tCalcL)/
		(tCalcU*tCalcU*tCalcU*tCalcU));

    // The ratio
    if (1) {
      tVal = (tContentExpLike[ti]/tContentExpUnlike[ti]) - (tCalcL/tCalcU);
      tErrFull = tErrExp + tErrCalc;
      tChi2Func += (tVal*tVal/(tErrFull));   
    }
    // The difference
    else {
      tVal = (tContentExpLike[ti] - tContentExpUnlike[ti]) - (tCalcL - tCalcU);
      tChi2Func += (tVal*tVal/(hypot(tErr2ExpLike[ti], tErr2ExpUnlike[ti]) + hypot(tErr2Like[ti], tErr2Unlike[ti])));   
    }
  }
  tChi2 += tChi2Func;

  return tChi2;
}

double  
CFFitterQuad::GetChi2(double aPurity)
{
  double tChi2=0.;
  double tChi2Func;
  //  for (int tT=0; tT<mNPTypes; tT++) {
  const double* tContentLike   = mCalcCFs[0]->GetContent();
  const double* tContentUnlike = mCalcCFs[1]->GetContent();
  const double* tErr2Like      = mCalcCFs[0]->GetError2();
  const double* tErr2Unlike    = mCalcCFs[1]->GetError2();

  const double* tContentExpLike   = ((ExpCF *) mExpCFs[0])->GetContent();
  const double* tContentExpUnlike = ((ExpCF *) mExpCFs[1])->GetContent();
  const double* tErr2ExpLike      = mExpCFs[0]->GetError2();
  const double* tErr2ExpUnlike    = mExpCFs[1]->GetError2();

  double tVal, tCalcL, tCalcU, tCalcLErr, tCalcUErr;
  double tExpLErr, tExpUErr, tErrExp, tErrCalc, tErrFull;
  tChi2Func = 0.0;
  for(int ti=0; ti<((ExpCFQuad *) mExpCFs[0])->GetNFitBin(); ti++){
    tCalcL = (tContentLike[ti] - 1.0) * ((ExpCFQuad *) mExpCFs[0])->GetPurity(ti) * aPurity + 1.0;
    tCalcU = (tContentUnlike[ti] - 1.0) * ((ExpCFQuad *) mExpCFs[1])->GetPurity(ti) * aPurity + 1.0;
    tCalcLErr = sqrt(tErr2Like[ti]) * ((ExpCFQuad *) mExpCFs[0])->GetPurity(ti) * aPurity;
    tCalcUErr = sqrt(tErr2Unlike[ti]) * ((ExpCFQuad *) mExpCFs[1])->GetPurity(ti) * aPurity;
    tExpLErr = sqrt(tErr2ExpLike[ti]);
    tExpUErr = sqrt(tErr2ExpUnlike[ti]);
    tErrExp = ((tExpLErr*tExpLErr*tContentExpUnlike[ti]*tContentExpUnlike[ti] + 
		tExpUErr*tExpUErr*tContentExpLike[ti]*tContentExpLike[ti])/
	       (tContentExpUnlike[ti]*tContentExpUnlike[ti]*tContentExpUnlike[ti]*tContentExpUnlike[ti]));
    tErrCalc = ((tCalcLErr*tCalcLErr*tCalcU*tCalcU +
		 tCalcUErr*tCalcUErr*tCalcL*tCalcL)/
		(tCalcU*tCalcU*tCalcU*tCalcU));

    PRINT_DEBUG_3("Bin " << ti);
    PRINT_DEBUG_3("CalcL CalcU ExpL ExpU " << tCalcL << " " << tCalcL << " " << tContentExpLike[ti] << " " << tContentExpUnlike[ti]);
    PRINT_DEBUG_3("ErrCL ErrCU ErrEL ErrEU " << tCalcLErr << " " << tCalcUErr << " " << tExpLErr << " " << tExpUErr);

    // The ratio
    if (1) {
      tVal = (tContentExpLike[ti]/tContentExpUnlike[ti]) - (tCalcL/tCalcU);
      tErrFull = tErrExp + tErrCalc;
      tChi2Func += (tVal*tVal/(tErrFull));   
      PRINT_DEBUG_3("Val Err Chi2 " << tVal << " " << tErrFull << " " <<  (tVal*tVal/(tErrFull*tErrFull)));
    }
    // The difference
    else {
      tVal = (tContentExpLike[ti] - tContentExpUnlike[ti]) - (tCalcL - tCalcU);
      tChi2Func += (tVal*tVal/((hypot(tErr2ExpLike[ti], tErr2ExpUnlike[ti]) + hypot(tErr2Like[ti], tErr2Unlike[ti]))*aPurity*aPurity));   
    }
  }
  tChi2 += tChi2Func;
  PRINT_DEBUG_3("Chi2 " << tChi2);

  return tChi2;
}

void    
CFFitterQuad::GenerateCF()
{
  for (int ti=0; ti<mNPTypes; ti++)
    mCalcCFs[ti]->Generate();
}

void     
CFFitterQuad::GeneratePairFileNameStub(fstream *os)
{
  (*os) << endl
	<< "# These are the names of the files containing input pairs" << endl
	<< "# Fill as many names as You need" << endl
	<< "# Leave others blank" << endl
	<< "InPairCalcNameLike = " << endl 
	<< "InPairCalcNameUnlike = " << endl 
	<< endl;

}

void 
CFFitterQuad::GenerateParameterStub()
{
  CFFitter::GenerateParameterStub();

  fstream *os;
  os = new fstream(sRPFileName.Data(),ios::ate|ios::out|ios::in);

  (*os) 
    << "# Quad Fitter specific parameters" << endl 
    << "# Fill as many pair types as You need" << endl 
    << "# Leave others blank" << endl 
    << "PairTypeUnlike = " << endl 
    << "PairTypeLike = " << endl 
  (*os) 
    << "# Minimum number of pairs in bin" << endl 
    << "MinPairsPerBin = 10000" << endl << endl;
  os->close();

  CalcCFQuad::GenerateParameterStub();
  ExpCFQuad::GenerateParameterStub();
  
}

void 
CFFitterQuad::ReadParameters()
{
  CFFitter::ReadParameters();
  STR tType;
  
  try {
    mPType = atoi((sRPInstance->getPar("PairTypeLike")).Data());
    mPTypes[0] = atoi((sRPInstance->getPar("PairTypeLike")).Data());
    mNPTypes = 1;
    tType = sRPInstance->getPar("PairTypeUnlike");
    mPTypes[1] = atoi(tType.Data());
    mNPTypes = 2;

    mPFileNames[0] = sRPInstance->getPar("InPairCalcNameLike");
    mPFileNames[1] = sRPInstance->getPar("InPairCalcNameUnlike");

    mMomRes  = atof((sRPInstance->getPar("MomentumResolutionCorrection")).Data())/100.0;

    mMinPairsPerBin = atoi((sRPInstance->getPar("MinPairsPerBin")).Data());
  }
  catch (STR e){
    PRINT_MESSAGE("CFFitterQuad::ReadParameters - Error reading parameter: " << e.Data());
    exit(0);
  }
  try { 
    mRandomSeed = atoi((sRPInstance->getPar("RandomSeed")).Data());
  }
  catch (STR e){
    PRINT_MESSAGE("CFFitterQuad::ReadParameters - Parameter not found: " << e.Data() << " Using default.");
    mRandomSeed = 21341;
    
  }
}

