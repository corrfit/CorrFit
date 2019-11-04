#include <TH2D.h>
#include <TFile.h>
#include "sf.h"
#include <math.h>

int main(int argc, char **argv)
{
  TH2D *histil = new TH2D("histil","histil",100,-0.2, 0.2, 100, -0.2, 0.2);
  
  for (int ix=1; ix<=histil->GetNbinsX(); ix++) {
    double ko = histil->GetXaxis()->GetBinCenter(ix);
    for (int iy=1; iy<=histil->GetNbinsY(); iy++) {
      double ks = histil->GetYaxis()->GetBinCenter(iy);
      double phi = atan2(ks, ko);
      histil->SetBinContent(ix, iy, 1*SpherHarmonics::ReYlm(0,0,0,1,phi) - 0.006*SpherHarmonics::ReYlm(1,1,1,phi) + 0.006*SpherHarmonics::ReYlm(1,-1,1,phi) + (0.0007)*SpherHarmonics::ImYlm(1,1,1,phi)+(0.0007)*SpherHarmonics::ImYlm(1,-1,1,phi));
    }
  }
  TFile *ofile = new TFile("ofil.root","RECREATE");
  ofile->cd();
  histil->Write();

}
