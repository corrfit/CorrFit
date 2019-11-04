#include <TH1D.h>
#include "sf.h"
#include <complex>
#include <TFile.h>
#include <iostream>
#include <TMath.h>

using namespace std;

int main(int argc, char **argv)
{
  TFile *infile = new TFile(argv[1]);
  
  TH1D *numsr[16];
  TH1D *densr[16];
  TH1D *cfnsr[16];

  TH1D *numsi[16];
  TH1D *densi[16];
  TH1D *cfnsi[16];


  TH1D *cfcs[16];

  char buf[200];
  int nh=0;

  for (int il=0; il<=2; il++) {
    for (int im=0; im<=2*il; im++) {
      sprintf(buf,"NumReYlm%i%iNonIdCYlm",il,im);
      numsr[nh] = (TH1D *) infile->Get(buf);
      sprintf(buf,"DenReYlm%i%iNonIdCYlm",il,im);
      densr[nh] = (TH1D *) infile->Get(buf);
      sprintf(buf,"CfnReYlm%i%iNonIdCYlm",il,im);
      cfnsr[nh] = (TH1D *) infile->Get(buf);

      sprintf(buf,"NumImYlm%i%iNonIdCYlm",il,im);
      numsi[nh] = (TH1D *) infile->Get(buf);
      sprintf(buf,"DenImYlm%i%iNonIdCYlm",il,im);
      densi[nh] = (TH1D *) infile->Get(buf);
      sprintf(buf,"CfnImYlm%i%iNonIdCYlm",il,im);
      cfnsi[nh] = (TH1D *) infile->Get(buf);

      cfcs[nh] = new TH1D(*cfnsr[nh]);
      cfcs[nh]->Reset("ICE");
      nh++;
    }
  }

  double phis[10] = { 0.1, 0.65, 2.33, 0.14, 4.33, 1.77, 2.99, 1.92, 0.99, 1.32 };
  double thet[10] = { 0.11,0.45, 0.66, 0.93, -0.22, -0.33, 0.12, 0.73, 0.12, -0.99 };

  int ibin = 6;

  complex<double> ylmbuffer[16];
  for (int iter=0; iter<10; iter++) {
    double qvec = numsr[0]->GetXaxis()->GetBinCenter(ibin);
    SpherHarmonics::YlmUpToL(3, thet[iter], phis[iter], ylmbuffer);
    double tsumr=0, csumr=0, msumr=0;
    double tsumi=0, csumi=0, msumi=0;

    for (int ilm=0; ilm<9; ilm++) {
      tsumr += numsr[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
      csumr += cfnsr[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
      msumr += densr[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);

      tsumr -= numsi[ilm]->GetBinContent(ibin)*imag(ylmbuffer[ilm]);
      csumr -= cfnsi[ilm]->GetBinContent(ibin)*imag(ylmbuffer[ilm]);
      msumr -= densi[ilm]->GetBinContent(ibin)*imag(ylmbuffer[ilm]);

      tsumi += numsr[ilm]->GetBinContent(ibin)*imag(ylmbuffer[ilm]);
      csumi += cfnsr[ilm]->GetBinContent(ibin)*imag(ylmbuffer[ilm]);
      msumi += densr[ilm]->GetBinContent(ibin)*imag(ylmbuffer[ilm]);

      tsumi += numsi[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
      csumi += cfnsi[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
      msumi += densi[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
    }
    
    cout << "t c*m " << tsumr/TMath::Sqrt(4*TMath::Pi()) << " "<< (csumr*msumr - csumi*msumi) << " " << tsumr/TMath::Sqrt(4*TMath::Pi())/( csumr*msumr - csumi*msumi) << endl;
    cout << "t c*m " << tsumi/TMath::Sqrt(4*TMath::Pi()) << " "<< (csumr*msumi + csumi*msumr) << " " << tsumi/TMath::Sqrt(4*TMath::Pi())/( csumr*msumi + csumi*msumr) << endl;
  }
}
