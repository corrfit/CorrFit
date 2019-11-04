#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TMinuit.h>
#include <TVectorD.h>

using namespace std;

TH2D    *gHisto2Fit;
TMinuit *gMinuitC2M;

void Chi2Mapfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  static TH2D* tH2;
  static double tSigMin;
  static double tSigMax;
  static double tMuMin;
  static double tMuMax;
  static double tMax;
  static double tSigNBin;
  static double tMuNBin;
  static double tSigStep;
  static double tMuStep;
  if(tH2 != gHisto2Fit){
    if(gHisto2Fit==0) {f=1e6;return;}
    tH2 = gHisto2Fit;
    tSigNBin = tH2->GetNbinsX();
    tSigMin = tH2->GetXaxis()->GetBinCenter(1);
    tSigMax = tH2->GetXaxis()->GetBinCenter((int) tSigNBin);
    tMuNBin = tH2->GetNbinsY();
    tMuMin = tH2->GetYaxis()->GetBinCenter(1);
    tMuMax = tH2->GetYaxis()->GetBinCenter((int) tMuNBin);
    tMax = tH2->GetMaximum();
    tSigStep = (tSigMax-tSigMin)/tSigNBin;
    tMuStep = (tMuMax-tMuMin)/tMuNBin;
  }
  if(par[0]<tSigMin || par[0]>tSigMax || par[1]<tMuMin || par[1]>tMuMax){
    double tMuEdge = par[1]-(tMuMax+tMuMin)/2.;
    double tSigEdge = par[1]-(tSigMax+tSigMin)/2.;
    f = tMax + sqrt(tSigEdge*tSigEdge+tMuEdge*tMuEdge);
  }
  else{ // quad interpolation
    int tSigBin = (int) ((par[0]-tSigMin)/tSigStep);
    double tSigMant = par[0] - tSigMin  - tSigBin*tSigStep;
    int tMuBin = (int) ((par[1]-tMuMin)/tMuStep);
    double tMuMant = par[1] - tMuMin  - tMuBin*tMuStep;
    double tSig0 = tSigMant;//*tSigMant ;
    double tSig1 = (tSigStep-tSigMant);//*(tSigStep-tSigMant);
    double tMu0 = tMuMant;//*tMuMant;
    double tMu1 = (tMuStep-tMuMant);//*(tMuStep-tMuMant);
    f=(tSig1*tMu1*tH2->GetCellContent(tSigBin,tMuBin)+
       tSig0*tMu1*tH2->GetCellContent(tSigBin+1,tMuBin)+
       tSig1*tMu0*tH2->GetCellContent(tSigBin,tMuBin+1)+
       tSig0*tMu0*tH2->GetCellContent(tSigBin+1,tMuBin+1))/tSigStep/tMuStep;
  }
}

void FindMin(TH2D* aH2, 
	     double& aSig, double& aSigErr, double& aMu, double& aMuErr,
	     double& aChi2Min){
  gHisto2Fit = aH2;
  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 1000;
  arglist[1] = 10.;
  gMinuitC2M->mnparm(0, "sig", aSig, 0.01, 0.,50.,ierflg);
  gMinuitC2M->mnparm(1, "mu",  aMu, 0.01, -20.,20.,ierflg);
  gMinuitC2M->mnexcm("MINIMIZE", arglist ,2,ierflg);
  gMinuitC2M->mnexcm("MINOS", arglist ,2,ierflg);
  gMinuitC2M->GetParameter(0, aSig, aSigErr);
  gMinuitC2M->GetParameter(1, aMu, aMuErr); 
  double tPar[2];
  tPar[0] = aSig;
  tPar[1] = aMu;
  int tNPar = 2;
  Chi2Mapfcn(tNPar, 0, aChi2Min, tPar, 0);
}

int main(int argc, char** argv)
{
  double mSig, mSigErr, mMu, mMuErr, mChiMin;
  Double_t *pars;
  
  if (argc < 3) {
    cout << "Usage: GetMin <filename> <histname>" << endl;
    exit(0);
  }

  TFile *infile = new TFile(argv[1]);
  TH2D  *inhist = (TH2D *) infile->Get(argv[2]);
  TVectorD *parstart = (TVectorD *) infile->Get("TVectorD");
  
  gMinuitC2M = new TMinuit(2);  
  gMinuitC2M->SetFCN(Chi2Mapfcn);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 0;
  gMinuitC2M->mnexcm("SET STR", arglist ,1,ierflg);
  arglist[0] = 1.;
  gMinuitC2M->mnexcm("SET ERR", arglist ,1,ierflg);
  
  
  mSig = (*parstart)(0);
  mMu  = (*parstart)(1);
  mSigErr = 0.1;
  mMuErr = 0.1;
  mChiMin = 1000.0;
  
  FindMin(inhist, mSig, mSigErr, mMu, mMuErr, mChiMin);
  
  cout << "Sigma " << mSig << " +/- " << mSigErr << endl;
  cout << "Mean "  << mMu  << " +/- " << mMuErr  << endl;
  cout << "Chi2Min " << mChiMin << endl;

  return(1);
}
