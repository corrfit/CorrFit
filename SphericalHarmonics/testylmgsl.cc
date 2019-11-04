#include <iostream>
#include "sf.h"
#include "mmatrix.h"
#include <TH1D.h>
#include <TRandom2.h>
#include <math.h>
#include <TFile.h>
#include <complex>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>

using namespace std;

int getbin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag)
{
  return (qbin*MMatrix::GetMaxJM()*MMatrix::GetMaxJM()*4 +
	  (ilmprim*2 + primimag) * MMatrix::GetMaxJM()*2 +
	  ilmzero*2 + zeroimag);
}

int main(int argc, char **argv)
{
  TH1D *shn00r = new TH1D("shn00r","shn00r",50,0.0,0.3);
  TH1D *shn00i = new TH1D("shn00i","shn00i",50,0.0,0.3);
  TH1D *shn12r = new TH1D("shn12r","shn12r",50,0.0,0.3);
  TH1D *shn12i = new TH1D("shn12i","shn12i",50,0.0,0.3);
  TH1D *shn10r = new TH1D("shn10r","shn10r",50,0.0,0.3);
  TH1D *shn10i = new TH1D("shn10i","shn10i",50,0.0,0.3);
  TH1D *shn11r = new TH1D("shn11r","shn11r",50,0.0,0.3);
  TH1D *shn11i = new TH1D("shn11i","shn11i",50,0.0,0.3);
  TH1D *shn24r = new TH1D("shn24r","shn24r",50,0.0,0.3);
  TH1D *shn24i = new TH1D("shn24i","shn24i",50,0.0,0.3);
  TH1D *shn23r = new TH1D("shn23r","shn23r",50,0.0,0.3);
  TH1D *shn23i = new TH1D("shn23i","shn23i",50,0.0,0.3);
  TH1D *shn20r = new TH1D("shn20r","shn20r",50,0.0,0.3);
  TH1D *shn20i = new TH1D("shn20i","shn20i",50,0.0,0.3);
  TH1D *shn21r = new TH1D("shn21r","shn21r",50,0.0,0.3);
  TH1D *shn21i = new TH1D("shn21i","shn21i",50,0.0,0.3);
  TH1D *shn22r = new TH1D("shn22r","shn22r",50,0.0,0.3);
  TH1D *shn22i = new TH1D("shn22i","shn22i",50,0.0,0.3);
  //   TH1D *shn30r = new TH1D("shn30r","shn30r",50,0.0,0.3);
  //   TH1D *shn30i = new TH1D("shn30i","shn30i",50,0.0,0.3);
  //   TH1D *shn31r = new TH1D("shn31r","shn31r",50,0.0,0.3);
  //   TH1D *shn31i = new TH1D("shn31i","shn31i",50,0.0,0.3);
  //   TH1D *shn32r = new TH1D("shn32r","shn32r",50,0.0,0.3);
  //   TH1D *shn32i = new TH1D("shn32i","shn32i",50,0.0,0.3);
  //   TH1D *shn33r = new TH1D("shn33r","shn33r",50,0.0,0.3);
  //   TH1D *shn33i = new TH1D("shn33i","shn33i",50,0.0,0.3);

  TH1D *shd00r = new TH1D("shd00r","shd00r",50,0.0,0.3);
  TH1D *shd00i = new TH1D("shd00i","shd00i",50,0.0,0.3);
  TH1D *shd12r = new TH1D("shd12r","shd12r",50,0.0,0.3);
  TH1D *shd12i = new TH1D("shd12i","shd12i",50,0.0,0.3);
  TH1D *shd10r = new TH1D("shd10r","shd10r",50,0.0,0.3);
  TH1D *shd10i = new TH1D("shd10i","shd10i",50,0.0,0.3);
  TH1D *shd11r = new TH1D("shd11r","shd11r",50,0.0,0.3);
  TH1D *shd11i = new TH1D("shd11i","shd11i",50,0.0,0.3);
  TH1D *shd24r = new TH1D("shd24r","shd24r",50,0.0,0.3);
  TH1D *shd24i = new TH1D("shd24i","shd24i",50,0.0,0.3);
  TH1D *shd23r = new TH1D("shd23r","shd23r",50,0.0,0.3);
  TH1D *shd23i = new TH1D("shd23i","shd23i",50,0.0,0.3);
  TH1D *shd20r = new TH1D("shd20r","shd20r",50,0.0,0.3);
  TH1D *shd20i = new TH1D("shd20i","shd20i",50,0.0,0.3);
  TH1D *shd21r = new TH1D("shd21r","shd21r",50,0.0,0.3);
  TH1D *shd21i = new TH1D("shd21i","shd21i",50,0.0,0.3);
  TH1D *shd22r = new TH1D("shd22r","shd22r",50,0.0,0.3);
  TH1D *shd22i = new TH1D("shd22i","shd22i",50,0.0,0.3);
  //   TH1D *shd30r = new TH1D("shd30r","shd30r",50,0.0,0.3);
  //   TH1D *shd30i = new TH1D("shd30i","shd30i",50,0.0,0.3);
  //   TH1D *shd31r = new TH1D("shd31r","shd31r",50,0.0,0.3);
  //   TH1D *shd31i = new TH1D("shd31i","shd31i",50,0.0,0.3);
  //   TH1D *shd32r = new TH1D("shd32r","shd32r",50,0.0,0.3);
  //   TH1D *shd32i = new TH1D("shd32i","shd32i",50,0.0,0.3);
  //   TH1D *shd33r = new TH1D("shd33r","shd33r",50,0.0,0.3);
  //   TH1D *shd33i = new TH1D("shd33i","shd33i",50,0.0,0.3);

  TH1D *shc00r = new TH1D("shc00r","shc00r",50,0.0,0.3);
  TH1D *shc00i = new TH1D("shc00i","shc00i",50,0.0,0.3);
  TH1D *shc10r = new TH1D("shc10r","shc10r",50,0.0,0.3);
  TH1D *shc10i = new TH1D("shc10i","shc10i",50,0.0,0.3);
  TH1D *shc11r = new TH1D("shc11r","shc11r",50,0.0,0.3);
  TH1D *shc11i = new TH1D("shc11i","shc11i",50,0.0,0.3);
  TH1D *shc12r = new TH1D("shc12r","shc12r",50,0.0,0.3);
  TH1D *shc12i = new TH1D("shc12i","shc12i",50,0.0,0.3);
  TH1D *shc20r = new TH1D("shc20r","shc20r",50,0.0,0.3);
  TH1D *shc20i = new TH1D("shc20i","shc20i",50,0.0,0.3);
  TH1D *shc21r = new TH1D("shc21r","shc21r",50,0.0,0.3);
  TH1D *shc21i = new TH1D("shc21i","shc21i",50,0.0,0.3);
  TH1D *shc22r = new TH1D("shc22r","shc22r",50,0.0,0.3);
  TH1D *shc22i = new TH1D("shc22i","shc22i",50,0.0,0.3);
  TH1D *shc23r = new TH1D("shc23r","shc23r",50,0.0,0.3);
  TH1D *shc23i = new TH1D("shc23i","shc23i",50,0.0,0.3);
  TH1D *shc24r = new TH1D("shc24r","shc24r",50,0.0,0.3);
  TH1D *shc24i = new TH1D("shc24i","shc24i",50,0.0,0.3);

  TH1D *binctn = new TH1D("binctn","binctn",50,0.0,0.3);
  TH1D *binctd = new TH1D("binctd","binctd",50,0.0,0.3);

  shn00r->Sumw2();
  shn00i->Sumw2();
  shn12r->Sumw2();
  shn12i->Sumw2();
  shn10r->Sumw2();
  shn10i->Sumw2();
  shn11r->Sumw2();
  shn11i->Sumw2();
  shn24r->Sumw2();
  shn24i->Sumw2();
  shn23r->Sumw2();
  shn23i->Sumw2();
  shn20r->Sumw2();
  shn20i->Sumw2();
  shn21r->Sumw2();
  shn21i->Sumw2();
  shn22r->Sumw2();
  shn22i->Sumw2();
  //   shn30r->Sumw2();
  //   shn30i->Sumw2();
  //   shn31r->Sumw2();
  //   shn31i->Sumw2();
  //   shn32r->Sumw2();
  //   shn32i->Sumw2();
  //   shn33r->Sumw2();
  //   shn33i->Sumw2();
  binctn->Sumw2();

  // The error matrix 
  // It is a 2Nx2N matrix for each q bin !!!
  double deltlm[20*20*50];
  //  for (int iter=0; iter<50000; iter++) deltlm[iter] = 0.0;
  
  // The error matrix for the correlation function
  double delclm[20*50*50];

  shd00r->Sumw2();
  shd00i->Sumw2();
  shd12r->Sumw2();
  shd12i->Sumw2();
  shd10r->Sumw2();
  shd10i->Sumw2();
  shd11r->Sumw2();
  shd11i->Sumw2();
  shd24r->Sumw2();
  shd24i->Sumw2();
  shd23r->Sumw2();
  shd23i->Sumw2();
  shd20r->Sumw2();
  shd20i->Sumw2();
  shd21r->Sumw2();
  shd21i->Sumw2();
  shd22r->Sumw2();
  shd22i->Sumw2();
  //   shd30r->Sumw2();
  //   shd30i->Sumw2();
  //   shd31r->Sumw2();
  //   shd31i->Sumw2();
  //   shd32r->Sumw2();
  //   shd32i->Sumw2();
  //   shd33r->Sumw2();
  //   shd33i->Sumw2();
  binctd->Sumw2();

  shc00r->Sumw2();
  shc00i->Sumw2();
  shc10r->Sumw2();
  shc10i->Sumw2();
  shc11r->Sumw2();
  shc11i->Sumw2();
  shc12r->Sumw2();
  shc12i->Sumw2();
  shc20r->Sumw2();
  shc20i->Sumw2();
  shc21r->Sumw2();
  shc21i->Sumw2();
  shc22r->Sumw2();
  shc22i->Sumw2();
  shc23r->Sumw2();
  shc23i->Sumw2();
  shc24r->Sumw2();
  shc24i->Sumw2();

  TRandom2 *tRand = new TRandom2();
  tRand->SetSeed(736271);

  Double_t ko, ks, kl;
  Double_t kv, kphi, ktheta;
  Double_t ylm00t;
  Double_t ylm00cr, ylm00ci;
  Double_t ylm10cr, ylm10ci;
  Double_t ylm11cr, ylm11ci;
  Double_t ylm12cr, ylm12ci;
  Double_t ylm24cr, ylm24ci;
  Double_t ylm23cr, ylm23ci;
  Double_t ylm20cr, ylm20ci;
  Double_t ylm21cr, ylm21ci;
  Double_t ylm22cr, ylm22ci;
  Double_t ylm30cr, ylm30ci;
  Double_t ylm31cr, ylm31ci;
  Double_t ylm32cr, ylm32ci;
  Double_t ylm33cr, ylm33ci;

  complex<double> ylm00;
  complex<double> ylm12;
  complex<double> ylm10;
  complex<double> ylm11;
  complex<double> ylm24;
  complex<double> ylm23;
  complex<double> ylm20;
  complex<double> ylm21;
  complex<double> ylm22;
  complex<double> ylmc[10];

  Double_t wgt;
  Double_t radx = 6.0 / 0.197327;
  Double_t rady = 4.0 / 0.197327;
  Double_t radz = 7.0 / 0.197327;

  radx *= radx;
  rady *= rady;
  radz *= radz;

  int flip = 0;

  int nbins = shn00r->GetNbinsX();
  double wbin = (shn00r->GetXaxis()->GetXmax() - shn00r->GetXaxis()->GetXmin())/nbins;
  double xmin = shn00r->GetXaxis()->GetXmin();

  for (int iter=0; iter<20000000; iter++) {
    if (!((iter-1) % 100000))
      cout << "Pair " << iter<< endl;
    
    ko = (tRand->Rndm()*2.0 - 1.0) * 0.32 / sqrt(3.0);
    ks = (tRand->Rndm()*2.0 - 1.0) * 0.32 / sqrt(3.0);
    kl = (tRand->Rndm()*2.0 - 1.0) * 0.32 / sqrt(3.0);

    wgt = 1.0 + exp(-ko*ko*radx
		    -ks*ks*rady
		    -kl*kl*radz);

    //    cout << ko << " " << ks << " " << kl << endl;
    // Calculating omega
    kv = sqrt(ko*ko + ks*ks + kl*kl);
    //     kphi = atan2(ks, ko);
    //     ktheta = atan2(sqrt(ko*ko+ks*ks), kl);
    //     cout << kv << " " << kphi << " " << ktheta << endl;

    // Calculating Ylm
    //    ylm00t = SpherHarmonics::ReYlm(1, 1, ktheta, kphi);
    ylmc[0] = ylm00 = SpherHarmonics::Ylm(0, 0, ko, ks, kl);
    ylmc[1] = ylm12 = SpherHarmonics::Ylm(1, -1, ko, ks, kl);
    ylmc[2] = ylm10 = SpherHarmonics::Ylm(1, 0, ko, ks, kl);
    ylmc[3] = ylm11 = SpherHarmonics::Ylm(1, 1, ko, ks, kl);
    ylmc[4] = ylm24 = SpherHarmonics::Ylm(2, -2, ko, ks, kl);
    ylmc[5] = ylm23 = SpherHarmonics::Ylm(2, -1, ko, ks, kl);
    ylmc[6] = ylm20 = SpherHarmonics::Ylm(2, 0, ko, ks, kl);
    ylmc[7] = ylm21 = SpherHarmonics::Ylm(2, 1, ko, ks, kl);
    ylmc[8] = ylm22 = SpherHarmonics::Ylm(2, 2, ko, ks, kl);
    
    //     ylm00cr = SpherHarmonics::ReYlm(0, 0, ko, ks, kl);
    //     ylm00ci = SpherHarmonics::ImYlm(0, 0, ko, ks, kl);
    //     ylm11cr = SpherHarmonics::ReYlm(1, -1, ko, ks, kl);
    //     ylm11ci = SpherHarmonics::ImYlm(1, -1, ko, ks, kl);
    //     ylm10cr = SpherHarmonics::ReYlm(1, 0, ko, ks, kl);
    //     ylm10ci = SpherHarmonics::ImYlm(1, 0, ko, ks, kl);
    //     ylm11cr = SpherHarmonics::ReYlm(1, 1, ko, ks, kl);
    //     ylm11ci = SpherHarmonics::ImYlm(1, 1, ko, ks, kl);
    //     ylm24cr = SpherHarmonics::ReYlm(2, -2, ko, ks, kl);
    //     ylm24ci = SpherHarmonics::ImYlm(2, -2, ko, ks, kl);
    //     ylm23cr = SpherHarmonics::ReYlm(2, -1, ko, ks, kl);
    //     ylm23ci = SpherHarmonics::ImYlm(2, -1, ko, ks, kl);
    //     ylm20cr = SpherHarmonics::ReYlm(2, 0, ko, ks, kl);
    //     ylm20ci = SpherHarmonics::ImYlm(2, 0, ko, ks, kl);
    //     ylm21cr = SpherHarmonics::ReYlm(2, 1, ko, ks, kl);
    //     ylm21ci = SpherHarmonics::ImYlm(2, 1, ko, ks, kl);
    //     ylm22cr = SpherHarmonics::ReYlm(2, 2, ko, ks, kl);
    //     ylm22ci = SpherHarmonics::ImYlm(2, 2, ko, ks, kl);
    //     ylm30cr = SpherHarmonics::ReYlm(3, 0, ko, ks, kl);
    //     ylm30ci = SpherHarmonics::ImYlm(3, 0, ko, ks, kl);
    //     ylm31cr = SpherHarmonics::ReYlm(3, 1, ko, ks, kl);
    //     ylm31ci = SpherHarmonics::ImYlm(3, 1, ko, ks, kl);
    //     ylm32cr = SpherHarmonics::ReYlm(3, 2, ko, ks, kl);
    //     ylm32ci = SpherHarmonics::ImYlm(3, 2, ko, ks, kl);
    //     ylm33cr = SpherHarmonics::ReYlm(3, 3, ko, ks, kl);
    //     ylm33ci = SpherHarmonics::ImYlm(3, 3, ko, ks, kl);

    //    cout << "Ylm is " << ylm00t << " " << ylm00c << endl << endl;

    // Filling the numerators
    if (flip) {
      if (wgt > tRand->Rndm()*2.0) {
	shn00r->Fill(kv, real(ylm00));
	shn00i->Fill(kv, -imag(ylm00));
	shn12r->Fill(kv, real(ylm12));
	shn12i->Fill(kv, -imag(ylm12));
	shn10r->Fill(kv, real(ylm10));
	shn10i->Fill(kv, -imag(ylm10));
	shn11r->Fill(kv, real(ylm11));
	shn11i->Fill(kv, -imag(ylm11));
	shn24r->Fill(kv, real(ylm24));
	shn24i->Fill(kv, -imag(ylm24));
	shn23r->Fill(kv, real(ylm23));
	shn23i->Fill(kv, -imag(ylm23));
	shn20r->Fill(kv, real(ylm20));
	shn20i->Fill(kv, -imag(ylm20));
	shn21r->Fill(kv, real(ylm21));
	shn21i->Fill(kv, -imag(ylm21));
	shn22r->Fill(kv, real(ylm22));
	shn22i->Fill(kv, -imag(ylm22));
	//       shn30r->Fill(kv, ylm30cr);
	//       shn30i->Fill(kv, ylm30ci);
	//       shn31r->Fill(kv, ylm31cr);
	//       shn31i->Fill(kv, ylm31ci);
	//       shn32r->Fill(kv, ylm32cr);
	//       shn32i->Fill(kv, ylm32ci);
	//       shn33r->Fill(kv, ylm33cr);
	//       shn33i->Fill(kv, ylm33ci);
	binctn->Fill(kv, 1.0);
	
	// Fill in the error matrix
	int nqbin = lrint(floor((kv - xmin)/wbin));
	int tabshift = nqbin*MMatrix::GetMaxJM()*MMatrix::GetMaxJM()*4;
	for (int ilmzero=0; ilmzero<MMatrix::GetMaxJM(); ilmzero++)
	  for (int ilmprim=0; ilmprim<MMatrix::GetMaxJM(); ilmprim++) {
	    // 	    deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2] += real(ylmc[ilmzero])*real(ylmc[ilmprim]);
	    // 	    deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2] += real(ylmc[ilmzero])*-imag(ylmc[ilmprim]);
	    // 	    deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] += -imag(ylmc[ilmzero])*real(ylmc[ilmprim]);
	    // 	    deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] += -imag(ylmc[ilmzero])*-imag(ylmc[ilmprim]);
	    deltlm[getbin(nqbin, ilmzero, 0, ilmprim, 0)] += real(ylmc[ilmzero])*real(ylmc[ilmprim]);
	    deltlm[getbin(nqbin, ilmzero, 0, ilmprim, 1)] += real(ylmc[ilmzero])*-imag(ylmc[ilmprim]);
	    deltlm[getbin(nqbin, ilmzero, 1, ilmprim, 0)] += -imag(ylmc[ilmzero])*real(ylmc[ilmprim]);
	    deltlm[getbin(nqbin, ilmzero, 1, ilmprim, 1)] += -imag(ylmc[ilmzero])*-imag(ylmc[ilmprim]);
	    
	  }
      }
      flip = 0;
    }
    else {
      // Filling denominators
      shd00r->Fill(kv, real(ylm00));
      shd00i->Fill(kv, -imag(ylm00));
      shd12r->Fill(kv, real(ylm12));
      shd12i->Fill(kv, -imag(ylm12));
      shd10r->Fill(kv, real(ylm10));
      shd10i->Fill(kv, -imag(ylm10));
      shd11r->Fill(kv, real(ylm11));
      shd11i->Fill(kv, -imag(ylm11));
      shd24r->Fill(kv, real(ylm24));
      shd24i->Fill(kv, -imag(ylm24));
      shd23r->Fill(kv, real(ylm23));
      shd23i->Fill(kv, -imag(ylm23));
      shd20r->Fill(kv, real(ylm20));
      shd20i->Fill(kv, -imag(ylm20));
      shd21r->Fill(kv, real(ylm21));
      shd21i->Fill(kv, -imag(ylm21));
      shd22r->Fill(kv, real(ylm22));
      shd22i->Fill(kv, -imag(ylm22));
      //    shd30r->Fill(kv, ylm30cr);
      //     shd30i->Fill(kv, ylm30ci);
      //     shd31r->Fill(kv, ylm31cr);
      //     shd31i->Fill(kv, ylm31ci);
      //     shd32r->Fill(kv, ylm32cr);
      //     shd32i->Fill(kv, ylm32ci);
      //     shd33r->Fill(kv, ylm33cr);
      //     shd33i->Fill(kv, ylm33ci);
      binctd->Fill(kv, 1.0);
      flip = 1;
    }
  }

  // Final normalization
  //   shn00r->Divide(binctn);
  //   shn00i->Divide(binctn);
  //   shn11r->Divide(binctn);
  //   shn11i->Divide(binctn);
  //   shn20r->Divide(binctn);
  //   shn20i->Divide(binctn);
  //   shn22r->Divide(binctn);
  //   shn22i->Divide(binctn);
  
  //   shd00r->Divide(binctd);
  //   shd00i->Divide(binctd);
  //   shd11r->Divide(binctd);
  //   shd11i->Divide(binctd);
  //   shd20r->Divide(binctd);
  //   shd20i->Divide(binctd);
  //   shd22r->Divide(binctd);
  //   shd22i->Divide(binctd);
  
  TFile *ofile = new TFile("shmout.root","RECREATE");
  
  ofile->cd();
  shn00r->Write();
  shn00i->Write();
  shn10r->Write();
  shn10i->Write();
  shn11r->Write();
  shn11i->Write();
  shn20r->Write();
  shn20i->Write();
  shn21r->Write();
  shn21i->Write();
  shn22r->Write();
  shn22i->Write();
  //  shn30r->Write();
  //   shn30i->Write();
  //   shn31r->Write();
  //   shn31i->Write();
  //   shn32r->Write();
  //   shn32i->Write();
  //   shn33r->Write();
  //   shn33i->Write();
  binctn->Write();

  shd00r->Write();
  shd00i->Write();
  shd12r->Write();
  shd12i->Write();
  shd10r->Write();
  shd10i->Write();
  shd11r->Write();
  shd11i->Write();
  shd24r->Write();
  shd24i->Write();
  shd23r->Write();
  shd23i->Write();
  shd20r->Write();
  shd20i->Write();
  shd21r->Write();
  shd21i->Write();
  shd22r->Write();
  shd22i->Write();
  //   shd30r->Write();
  //   shd30i->Write();
  //   shd31r->Write();
  //   shd31i->Write();
  //   shd32r->Write();
  //   shd32i->Write();
  //   shd33r->Write();
  //   shd33i->Write();
  binctd->Write();

  complex<double> tMq0[10];
  complex<double> tTq0[10];
  double tMq0packed[20];
  double tTq0packed[20];
  double tMTilde[20*20];
  double tMTildeInv[20*20];
  complex<double> tCq0[10];

  for (int ibin=1; ibin<=shc00r->GetNbinsX(); ibin++) {
    //  for (int ibin=15; ibin<16; ibin++) {

    if (0) {
      tMq0[0] = shd00r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[1] = shd12r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[2] = shd10r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[3] = shd11r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[4] = shd24r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[5] = shd23r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[6] = shd20r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[7] = shd21r->GetBinContent(20)/shd00r->GetEntries();
      tMq0[8] = shd22r->GetBinContent(20)/shd00r->GetEntries();
      //  tMq0[9] = shd33r->GetBinContent(20)/shd00r->GetEntries();
    }

    tMq0[0] = complex<double>(shd00r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd00i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[1] = complex<double>(shd12r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd12i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[2] = complex<double>(shd10r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd10i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[3] = complex<double>(shd11r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd11i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[4] = complex<double>(shd24r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd24i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[5] = complex<double>(shd23r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd23i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[6] = complex<double>(shd20r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd20i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[7] = complex<double>(shd21r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd21i->GetBinContent(ibin)/shd00r->GetEntries());
    tMq0[8] = complex<double>(shd22r->GetBinContent(ibin)/shd00r->GetEntries(),
			      shd22i->GetBinContent(ibin)/shd00r->GetEntries());

    tTq0[0] = complex<double>(shn00r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn00i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[1] = complex<double>(shn12r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn12i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[2] = complex<double>(shn10r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn10i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[3] = complex<double>(shn11r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn11i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[4] = complex<double>(shn24r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn24i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[5] = complex<double>(shn23r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn23i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[6] = complex<double>(shn20r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn20i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[7] = complex<double>(shn21r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn21i->GetBinContent(ibin)/shn00r->GetEntries());
    tTq0[8] = complex<double>(shn22r->GetBinContent(ibin)/shn00r->GetEntries(),
			      shn22i->GetBinContent(ibin)/shn00r->GetEntries());

    // Calculate the proper error matrix for T
    int tabshift = (ibin-1)*MMatrix::GetMaxJM()*MMatrix::GetMaxJM()*4;
    for (int ilmzero=0; ilmzero<MMatrix::GetMaxJM(); ilmzero++)
      for (int ilmprim=0; ilmprim<MMatrix::GetMaxJM(); ilmprim++) {
	// 	deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2] /= shn00r->GetEntries();
	// 	deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2] /= shn00r->GetEntries();
	// 	deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] /= shn00r->GetEntries();
	// 	deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] /= shn00r->GetEntries();

	// 	deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2] -= real(tTq0[ilmzero])*real(tTq0[ilmprim]);
	// 	deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2] -= real(tTq0[ilmzero])*imag(tTq0[ilmprim]);
	// 	deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] -= imag(tTq0[ilmzero])*real(tTq0[ilmprim]);
	// 	deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] -= imag(tTq0[ilmzero])*imag(tTq0[ilmprim]);

	// 	deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2] /= (shn00r->GetEntries() - 1);
	// 	deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2] /= (shn00r->GetEntries() - 1);
	// 	deltlm[tabshift + ((ilmprim*2)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] /= (shn00r->GetEntries() - 1);
	// 	deltlm[tabshift + ((ilmprim*2+1)*MMatrix::GetMaxJM()*2) + ilmzero*2+1] /= (shn00r->GetEntries() - 1);

	deltlm[getbin(ibin-1, ilmzero, 0, ilmprim, 0)] /= shn00r->GetEntries();
	deltlm[getbin(ibin-1, ilmzero, 0, ilmprim, 1)] /= shn00r->GetEntries();
	deltlm[getbin(ibin-1, ilmzero, 1, ilmprim, 0)] /= shn00r->GetEntries();
	deltlm[getbin(ibin-1, ilmzero, 1, ilmprim, 1)] /= shn00r->GetEntries();

	deltlm[getbin(ibin-1, ilmzero, 0, ilmprim, 0)] -= real(tTq0[ilmzero])*real(tTq0[ilmprim]);
	deltlm[getbin(ibin-1, ilmzero, 0, ilmprim, 1)] -= real(tTq0[ilmzero])*imag(tTq0[ilmprim]);
	deltlm[getbin(ibin-1, ilmzero, 1, ilmprim, 0)] -= imag(tTq0[ilmzero])*real(tTq0[ilmprim]);
	deltlm[getbin(ibin-1, ilmzero, 1, ilmprim, 1)] -= imag(tTq0[ilmzero])*imag(tTq0[ilmprim]);

	deltlm[getbin(ibin-1, ilmzero, 0, ilmprim, 0)] /= (shn00r->GetEntries() - 1);
	deltlm[getbin(ibin-1, ilmzero, 0, ilmprim, 1)] /= (shn00r->GetEntries() - 1);
	deltlm[getbin(ibin-1, ilmzero, 1, ilmprim, 0)] /= (shn00r->GetEntries() - 1);
	deltlm[getbin(ibin-1, ilmzero, 1, ilmprim, 1)] /= (shn00r->GetEntries() - 1);
      }

    // Write out error matrices
    cout << "Errors from T " << endl;
    cout << shn00r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn00i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn12r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn12i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn10r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn11i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn11r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn10i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn24r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn24i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn23r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn23i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn20r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn20i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn21r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn21i->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn22r->GetBinError(ibin)/shn00r->GetEntries() << " " 
	 << shn22i->GetBinError(ibin)/shn00r->GetEntries() << endl;
    cout << endl;
    cout << "T error matrix" << endl;
    for (int ilmprim=0; ilmprim<MMatrix::GetMaxJM()*2; ilmprim++) {
      for (int ilmzero=0; ilmzero<MMatrix::GetMaxJM()*2; ilmzero++) {
	int ilmp = ilmprim / 2;
	int ilmz = ilmzero / 2;
	int imgp = ilmprim % 2;
	int imgz = ilmzero % 2;
	cout.width(10);
	//	cout << sqrt(fabs(deltlm[tabshift + (ilmprim*MMatrix::GetMaxJM()*2) + ilmzero]));
	cout << sqrt(fabs(deltlm[getbin(ibin-1, ilmz, imgz, ilmp, imgp)]));
      }
      cout << endl;
    }
    cout << endl;
    
    for (int iter=0; iter<10; iter++) {
      tMq0packed[iter*2] = real(tMq0[iter]);
      tMq0packed[iter*2+1] = imag(tMq0[iter]);

      tTq0packed[iter*2] = real(tTq0[iter]);
      tTq0packed[iter*2+1] = imag(tTq0[iter]);
    }

    cout << "20 bin of M " << endl;
    for (int iter=0; iter<10; iter++)
      cout << tMq0[iter] << "   ";
    cout << endl << endl;

    for (int iter=0; iter<10; iter++)
      cout << tMq0packed[iter*2] << " " << tMq0packed[iter*2+1] << " ";
    cout << endl << endl;
    
    cout << "20 bin of T " << endl;
    for (int iter=0; iter<10; iter++)
      cout << tTq0[iter] << "   ";
    cout << endl << endl;
    
    MMatrix::GetMtilde(tMq0, tMTilde);
  
    cout << "M matrix " << endl;
    for (int row=0; row<MMatrix::GetMaxJM()*2; row++) {
      for (int iter=0; iter<MMatrix::GetMaxJM()*2; iter++) {
	cout.width(10);
	cout.precision(3);
 	cout << tMTilde[row*(2*MMatrix::GetMaxJM())+iter]; 
	// 	cout << real(tMTilde[row*MMatrix::GetMaxJM()+iter]); 
	// 	if (imag(tMTilde[row*MMatrix::GetMaxJM()+iter])>0) cout << ",";
	// 	else cout << ",";
	// 	cout.width(10);
	// 	cout.precision(3);
	// 	cout << (imag(tMTilde[row*MMatrix::GetMaxJM()+iter])) << " | ";
      }
      cout << endl;
    }
    cout << endl;
    
    // Perform the solution for the correlation function itself and the errors

    // Fill in the M matrix
    double mM[20*20];
    for (int iter=0; iter<20*20; iter++)
      mM[iter] = tMTilde[iter];

    gsl_matrix_view matM = gsl_matrix_view_array(mM, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);

    // Prepare halper matrices
    double mU[20*20];
    for (int iter=0; iter<20*20; iter++)
      mU[iter] = tMTilde[iter];

    gsl_matrix_view matU = gsl_matrix_view_array(mU, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);

    double mV[20*20];
    gsl_matrix_view matV = gsl_matrix_view_array(mV, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);

    double vS[20];
    gsl_vector_view vecS = gsl_vector_view_array(vS, MMatrix::GetMaxJM()*2);

    double vW[20];
    gsl_vector_view vecW = gsl_vector_view_array(vW, MMatrix::GetMaxJM()*2);

    // Decomposing the M matrix
    gsl_linalg_SV_decomp(&matU.matrix, &matV.matrix, &vecS.vector, &vecW.vector);
  
    double vB[20];
    for (int iter=0; iter<MMatrix::GetMaxJM(); iter++) {
      vB[iter*2] = real(tTq0[iter]);
      vB[iter*2+1] = imag(tTq0[iter]);
    }

    // Prepare inputs for solving the problem
    gsl_vector_view vecB = gsl_vector_view_array(vB, MMatrix::GetMaxJM()*2);

    double vX[20] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    gsl_vector_view vecX = gsl_vector_view_array(vX, MMatrix::GetMaxJM()*2);

    // Solving the problem
    gsl_linalg_SV_solve(&matU.matrix, &matV.matrix, &vecS.vector, &vecB.vector, &vecX.vector);

    // Double-check the solution
    double vY[20] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    gsl_vector_view vecY = gsl_vector_view_array(vY, MMatrix::GetMaxJM()*2);

    gsl_blas_dgemv(CblasNoTrans, 1.0, &matM.matrix, &vecX.vector, 1.0, &vecY.vector); 
        
    //   double tIdent[10*10];

    //   cout << "M*M^-1 " << endl;
    //   for (int row=0; row<MMatrix::GetMaxJM(); row++) {
    //     for (int iter=0; iter<MMatrix::GetMaxJM(); iter++)
    //       cout << tIdent[row*MMatrix::GetMaxJM()+iter] << "\t";
    //     cout << endl << endl;
    //   }

    //   for (int ilmc=0; ilmc<MMatrix::GetMaxJM(); ilmc++) {
    //     tCq0[ilmc] = 0.0;
    //     for (int ilmt=0; ilmt<MMatrix::GetMaxJM(); ilmt++) {
    //       tCq0[ilmc] += tMTildeInv[ilmc*MMatrix::GetMaxJM() + ilmt] * tTq0[ilmt];
    //     }
    //     //    tCq0[ilmc] /= (shd00r->GetEntries()*shn00r->GetEntries());
    //   }
    
    cout << "20 bin of C " << endl;
    for (int iter=0; iter<MMatrix::GetMaxJM()*2; iter++)
      cout << vX[iter] << "   ";
    cout << endl << endl;
  
    cout << "20 bin of reconstructed M " << endl;
    for (int iter=0; iter<MMatrix::GetMaxJM()*2; iter++)
      cout << vY[iter] << "   ";
    cout << endl << endl;
  
    shc00r->SetBinContent(ibin, vX[0]);
    shc00i->SetBinContent(ibin, vX[1]);
    shc12r->SetBinContent(ibin, vX[2]);
    shc12i->SetBinContent(ibin, vX[3]);
    shc10r->SetBinContent(ibin, vX[4]);
    shc10i->SetBinContent(ibin, vX[5]);
    shc11r->SetBinContent(ibin, vX[6]);
    shc11i->SetBinContent(ibin, vX[7]);
    shc24r->SetBinContent(ibin, vX[8]);
    shc24i->SetBinContent(ibin, vX[9]);
    shc23r->SetBinContent(ibin, vX[10]);
    shc23i->SetBinContent(ibin, vX[11]);
    shc20r->SetBinContent(ibin, vX[12]);
    shc20i->SetBinContent(ibin, vX[13]);
    shc21r->SetBinContent(ibin, vX[14]);
    shc21i->SetBinContent(ibin, vX[15]);
    shc22r->SetBinContent(ibin, vX[16]);
    shc22i->SetBinContent(ibin, vX[17]);
  
    if (0) {
      // Perform the solutions for the error matrix as well
      
      for (int iter=0; iter<20*20; iter++)
	// We must always add errors
	mU[iter] = (tMTilde[iter]);

      // Decomposing the M matrix
      gsl_linalg_SV_decomp(&matU.matrix, &matV.matrix, &vecS.vector, &vecW.vector);

      for (int iterlmz=0; iterlmz<MMatrix::GetMaxJM()*2; iterlmz++) {
	int ilmz = iterlmz/2;
	int imgz = iterlmz%2;
	for (int ilm=0; ilm<MMatrix::GetMaxJM(); ilm++) {
	  vB[ilm*2]   = sqrt(fabs(deltlm[getbin(ibin-1, ilmz, imgz, ilm, 0)]));
	  vB[ilm*2+1] = sqrt(fabs(deltlm[getbin(ibin-1, ilmz, imgz, ilm, 1)]));
	}
	gsl_linalg_SV_solve(&matU.matrix, &matV.matrix, &vecS.vector, &vecB.vector, &vecX.vector);
      
	for (int ilm=0; ilm<MMatrix::GetMaxJM(); ilm++) {
	  if (!(isnan(vX[ilm*2])))
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 0)] = vX[ilm*2];
	  else
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 0)] = 1e-12;
	  if (!(isnan(vX[ilm*2+1])))
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 1)] = vX[ilm*2+1];
	  else
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 1)] = 1e-12;
	}

      }

      // Double-check the solution
      double mY[20*20];
      gsl_matrix_view matY = gsl_matrix_view_array(mY, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);
      gsl_matrix_view matC = gsl_matrix_view_array(delclm, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);
    
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &matM.matrix, &matC.matrix, 1.0, &matY.matrix); 
    
      //     cout << "T error matrix reconstructed" << endl;
      //     for (int ilmprim=0; ilmprim<MMatrix::GetMaxJM()*2; ilmprim++) {
      //       for (int ilmzero=0; ilmzero<MMatrix::GetMaxJM()*2; ilmzero++) {
      // 	int ilmp = ilmprim / 2;
      // 	int ilmz = ilmzero / 2;
      // 	int imgp = ilmprim % 2;
      // 	int imgz = ilmzero % 2;
      // 	cout.width(10);
      // 	//	cout << sqrt(fabs(deltlm[tabshift + (ilmprim*MMatrix::GetMaxJM()*2) + ilmzero]));
      // 	cout << sqrt(fabs(mY[getbin(ibin-1, ilmz, imgz, ilmp, imgp)]));
      //       }
      //       cout << endl;
      //     }
      //     cout << endl;
    }

    if (1) {
      // Perform the solutions for the error matrix as well
      // Second way of doing errors - MM matrix
      
      // Calculating the MM matrix
      double mM1[20*20];
      for (int iter=0; iter<20*20; iter++)
	// We must always add errors
	mM1[iter] = (tMTilde[iter]);
    
      gsl_matrix_view matM1 = gsl_matrix_view_array(mM1, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);
      // Calculating the M^2 = MM matrix
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &matM1.matrix, &matM1.matrix, 1.0, &matU.matrix); 
      
      // Decomposing the U=M^2 matrix
      gsl_linalg_SV_decomp(&matU.matrix, &matV.matrix, &vecS.vector, &vecW.vector);
    
      for (int iterlmz=0; iterlmz<MMatrix::GetMaxJM()*2; iterlmz++) {
	int ilmz = iterlmz/2;
	int imgz = iterlmz%2;
	for (int ilm=0; ilm<MMatrix::GetMaxJM(); ilm++) {
	  // Errors must be positive
	  vB[ilm*2]   = (fabs(deltlm[getbin(ibin-1, ilmz, imgz, ilm, 0)]));
	  vB[ilm*2+1] = (fabs(deltlm[getbin(ibin-1, ilmz, imgz, ilm, 1)]));
	}
	gsl_linalg_SV_solve(&matU.matrix, &matV.matrix, &vecS.vector, &vecB.vector, &vecX.vector);
	
	for (int ilm=0; ilm<MMatrix::GetMaxJM(); ilm++) {
	  if (!(isnan(vX[ilm*2])))
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 0)] = vX[ilm*2];
	  else
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 0)] = 1e-12;
	  if (!(isnan(vX[ilm*2+1])))
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 1)] = vX[ilm*2+1];
	  else
	    delclm[getbin(ibin-1, ilmz, imgz, ilm, 1)] = 1e-12;
	}
      }
      
      // Double-check the solution
      double mY[20*20];
      gsl_matrix_view matY = gsl_matrix_view_array(mY, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);
      gsl_matrix_view matC = gsl_matrix_view_array(delclm, MMatrix::GetMaxJM()*2, MMatrix::GetMaxJM()*2);
      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &matM.matrix, &matC.matrix, 1.0, &matY.matrix); 
      
      //     cout << "T error matrix reconstructed" << endl;
      //     for (int ilmprim=0; ilmprim<MMatrix::GetMaxJM()*2; ilmprim++) {
      //       for (int ilmzero=0; ilmzero<MMatrix::GetMaxJM()*2; ilmzero++) {
      // 	int ilmp = ilmprim / 2;
      // 	int ilmz = ilmzero / 2;
      // 	int imgp = ilmprim % 2;
      // 	int imgz = ilmzero % 2;
      // 	cout.width(10);
      // 	//	cout << sqrt(fabs(deltlm[tabshift + (ilmprim*MMatrix::GetMaxJM()*2) + ilmzero]));
      // 	cout << sqrt(fabs(mY[getbin(ibin-1, ilmz, imgz, ilmp, imgp)]));
      //       }
      //       cout << endl;
      //     }
      //     cout << endl;
      
    }
    
    shc00r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 0, 0, 0, 0)]))));
    shc00i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 0, 1, 0, 1)]))));
    shc12r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 1, 0, 1, 0)]))));
    shc12i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 1, 1, 1, 1)]))));
    shc10r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 2, 0, 2, 0)]))));
    shc10i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 2, 1, 2, 1)]))));
    shc11r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 3, 0, 3, 0)]))));
    shc11i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 3, 1, 3, 1)]))));
    shc24r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 4, 0, 4, 0)]))));
    shc24i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 4, 1, 4, 1)]))));
    shc23r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 5, 0, 5, 0)]))));
    shc23i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 5, 1, 5, 1)]))));
    shc20r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 6, 0, 6, 0)]))));
    shc20i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 6, 1, 6, 1)]))));
    shc21r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 7, 0, 7, 0)]))));
    shc21i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 7, 1, 7, 1)]))));
    shc22r->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 8, 0, 8, 0)]))));
    shc22i->SetBinError(ibin, (sqrt(fabs(delclm[getbin(ibin-1, 8, 1, 8, 1)]))));
    
    cout << "C error matrix" << endl;
    for (int ilmprim=0; ilmprim<MMatrix::GetMaxJM()*2; ilmprim++) {
      for (int ilmzero=0; ilmzero<MMatrix::GetMaxJM()*2; ilmzero++) {
	int ilmp = ilmprim / 2;
	int ilmz = ilmzero / 2;
	int imgp = ilmprim % 2;
	int imgz = ilmzero % 2;
	cout.width(10);
	//	cout << sqrt(fabs(deltlm[tabshift + (ilmprim*MMatrix::GetMaxJM()*2) + ilmzero]));
	cout << (sqrt(fabs(delclm[getbin(ibin-1, ilmz, imgz, ilmp, imgp)])));
      }
      cout << endl;
    }
    cout << endl;
    
  }

  shc00r->Write();
  shc00i->Write();
  shc10r->Write();
  shc10i->Write();
  shc11r->Write();
  shc11i->Write();
  shc12r->Write();
  shc12i->Write();
  shc20r->Write();
  shc20i->Write();
  shc21r->Write();
  shc21i->Write();
  shc22r->Write();
  shc22i->Write();
  shc23r->Write();
  shc23i->Write();
  shc24r->Write();
  shc24i->Write();
  
}

  
