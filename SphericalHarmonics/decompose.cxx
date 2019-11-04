#include "sf.h"
#include <TH3D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <complex>

int main(int argc, char **argv)
{
  TFile *infile = new TFile(argv[1]);

  TH3D *cnum = (TH3D *) infile->Get("NumFake3DLCMSSphModel");
  TH3D *cden = (TH3D *) infile->Get("Den3DLCMSSphModel");
  TH3D *crat = new TH3D(*cnum);
  crat->Divide(cden);
  crat->SetTitle("crat");
  crat->SetName("crat");

  
  TH1D *ylm00 = new TH1D("ylm00", "ylm00", crat->GetNbinsX(), crat->GetXaxis()->GetXmin(), crat->GetXaxis()->GetXmax());
  TH1D *ylm10 = new TH1D("ylm10", "ylm10", crat->GetNbinsX(), crat->GetXaxis()->GetXmin(), crat->GetXaxis()->GetXmax());
  TH1D *ylm11 = new TH1D("ylm11", "ylm11", crat->GetNbinsX(), crat->GetXaxis()->GetXmin(), crat->GetXaxis()->GetXmax());
  TH1D *ylm20 = new TH1D("ylm20", "ylm20", crat->GetNbinsX(), crat->GetXaxis()->GetXmin(), crat->GetXaxis()->GetXmax());
  TH1D *ylm21 = new TH1D("ylm21", "ylm21", crat->GetNbinsX(), crat->GetXaxis()->GetXmin(), crat->GetXaxis()->GetXmax());
  TH1D *ylm22 = new TH1D("ylm22", "ylm22", crat->GetNbinsX(), crat->GetXaxis()->GetXmin(), crat->GetXaxis()->GetXmax());

  ylm00->Sumw2();
  ylm10->Sumw2();
  ylm11->Sumw2();
  ylm20->Sumw2();
  ylm21->Sumw2();
  ylm22->Sumw2();
  
  double npt = crat->GetNbinsY()*crat->GetNbinsZ()/TMath::Sqrt(4*TMath::Pi());
  double phic11 = (TMath::Pi()/6.0)/(TMath::Sin(TMath::Pi()/6.0));
  double phic12 = (TMath::Pi()/3.0)/(TMath::Sin(TMath::Pi()/3.0));

  complex<double> ylmbuffer[9];

  double sumw[6];
  double sumw2[6];

  for (int ik = 1; ik<=crat->GetNbinsX(); ik++) {
    for (int iw=0; iw<6; iw++)
      sumw[iw] = sumw2[iw] = 0;
    for (int ip = 1; ip<=crat->GetNbinsY(); ip++) {
      for (int it = 1; it<=crat->GetNbinsZ(); it++) {
	SpherHarmonics::YlmUpToL(2, crat->GetZaxis()->GetBinCenter(it), crat->GetYaxis()->GetBinCenter(ip), ylmbuffer);
	ylm00->Fill(crat->GetXaxis()->GetBinCenter(ik), crat->GetBinContent(ik, ip, it)*real(ylmbuffer[0])/npt);
	ylm10->Fill(crat->GetXaxis()->GetBinCenter(ik), crat->GetBinContent(ik, ip, it)*real(ylmbuffer[2])/npt);
	ylm11->Fill(crat->GetXaxis()->GetBinCenter(ik), crat->GetBinContent(ik, ip, it)*real(ylmbuffer[3])*phic11/npt);
	ylm20->Fill(crat->GetXaxis()->GetBinCenter(ik), crat->GetBinContent(ik, ip, it)*real(ylmbuffer[6])/npt);
	ylm21->Fill(crat->GetXaxis()->GetBinCenter(ik), crat->GetBinContent(ik, ip, it)*real(ylmbuffer[7])*phic11/npt);
	ylm22->Fill(crat->GetXaxis()->GetBinCenter(ik), crat->GetBinContent(ik, ip, it)*real(ylmbuffer[8])*phic12/npt);
	sumw[0] += 1.0/TMath::Power(crat->GetBinError(ik, ip, it),2);
	sumw[1] += 1.0/TMath::Power(crat->GetBinError(ik, ip, it),2);
	sumw[2] += 1.0/TMath::Power(crat->GetBinError(ik, ip, it),2);
	sumw[3] += 1.0/TMath::Power(crat->GetBinError(ik, ip, it),2);
	sumw[4] += 1.0/TMath::Power(crat->GetBinError(ik, ip, it),2);
	sumw[5] += 1.0/TMath::Power(crat->GetBinError(ik, ip, it),2);
      }
    }
    ylm00->SetBinError(ik, TMath::Sqrt(1.0/sumw[0]));
    ylm10->SetBinError(ik, TMath::Sqrt(1.0/sumw[1]));
    ylm11->SetBinError(ik, TMath::Sqrt(1.0/sumw[2]));
    ylm20->SetBinError(ik, TMath::Sqrt(1.0/sumw[3]));
    ylm21->SetBinError(ik, TMath::Sqrt(1.0/sumw[4]));
    ylm22->SetBinError(ik, TMath::Sqrt(1.0/sumw[5]));
  }

  TFile *ofile = new TFile("decomposed.root","RECREATE");
  ofile->cd();
  ylm00->Write();
  ylm10->Write();
  ylm11->Write();
  ylm20->Write();
  ylm21->Write();
  ylm22->Write();
}
