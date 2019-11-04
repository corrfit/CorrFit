#include <TH1D.h>
#include <sf.h>
#include <complex>

using namespace std;

int main(int argc, char **argv)
{
  TFile *infile = new TFile(argv[0]);
  
  TH1D *nums[16];
  TH1D *dens[16];
  TH1D *cfns[16];
  TH1D *cfcs[16];

  char buf[200];
  int nh=0;

  for (int il=0; il<=3; il++) {
    for (int im=0; im<=2*il; il++) {
      sprintf("NumReYlm%i%iNonIdCYlm",il,im);
      nums[nh] = (TH1D *) infile->Get(buf);
      sprintf("DenReYlm%i%iNonIdCYlm",il,im);
      dens[nh] = (TH1D *) infile->Get(buf);
      sprintf("CfnReYlm%i%iNonIdCYlm",il,im);
      cfns[nh] = (TH1D *) infile->Get(buf);
      cfcs[nh] = new TH1D(*cfns[hn]);
      cfcs[nh]->Reset("ICE");
      nh++;
    }
  }

  double phis[10] = { 0.1, 0.65, 2.33, 0.14, 4.33, 1.77, 2.99, 1.92, 0.99, 1.32 };
  double thet[10] = { 0.11,0.45, 0.66, 0.93, -0.22, -0.33, 0.12, 0.73, 0.12, -0.99 };
  int ibin = 10;

  complex<double> ylmbuffer[16];
  for (int iter=0; iter<10; iter++) {
    double qvec = nums[0]->GetXaxis()->GetBinCenter(10);
    SpherHarmonics::YlmUpToL(3, thet[iter], phis[iter], ylmbuffer);
    double tsum=0, csum=0, msum=0;

    for (int ilm=0; ilm<16; ilm++) {
      tsum += nums[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
      csum += cfns[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
      msum += dens[ilm]->GetBinContent(ibin)*real(ylmbuffer[ilm]);
    }
    
    cout << "t c*m " << tsum << " "<< csum*msum << endl;
  }
}
