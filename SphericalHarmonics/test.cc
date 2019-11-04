#include <iostream>
#include "sf.h"
#include "mmatrix.h"
#include <TH1D.h>
#include <TRandom2.h>
#include <math.h>
#include <TFile.h>

using namespace std;

// long double getdetone(long double *aMatrix, int size, int rowx, int rowy)
// {
//   long double det = 0;
//   long double psum;
//   for (int iters=0; iters<size-1; iters++) {
//     psum = 1.0;
//     for (int iterk=0; iterk<size-1; iterk++) {
//       psum *= aMatrix[((iters+iterk+(iters+iterk>=rowx)) % 6) * 6 + iterk + (iterk>=rowy)];
//     }
//     det += psum;
//   }
  
//   for (int iters=0; iters<size-1; iters++) {
//     psum = 1.0;
//     for (int iterk=0; iterk<size-1; iterk++) {
//       psum *= aMatrix[((iters-iterk+(iters-iterk>=rowx)) % 6) * 6 + iterk + (iterk>=rowy)];
//     }
//     det += psum;
//   }

//   if ((rowx+rowy) % 2) det *= -1;
  
//   return det;
// }

int main(int argc, char **argv)
{
  TH1D *shn00r = new TH1D("shn00r","shn00r",100,0.0,0.3);
  TH1D *shn00i = new TH1D("shn00i","shn00i",100,0.0,0.3);
  TH1D *shn12r = new TH1D("shn12r","shn12r",100,0.0,0.3);
  TH1D *shn12i = new TH1D("shn12i","shn12i",100,0.0,0.3);
  TH1D *shn10r = new TH1D("shn10r","shn10r",100,0.0,0.3);
  TH1D *shn10i = new TH1D("shn10i","shn10i",100,0.0,0.3);
  TH1D *shn11r = new TH1D("shn11r","shn11r",100,0.0,0.3);
  TH1D *shn11i = new TH1D("shn11i","shn11i",100,0.0,0.3);
  TH1D *shn24r = new TH1D("shn24r","shn24r",100,0.0,0.3);
  TH1D *shn24i = new TH1D("shn24i","shn24i",100,0.0,0.3);
  TH1D *shn23r = new TH1D("shn23r","shn23r",100,0.0,0.3);
  TH1D *shn23i = new TH1D("shn23i","shn23i",100,0.0,0.3);
  TH1D *shn20r = new TH1D("shn20r","shn20r",100,0.0,0.3);
  TH1D *shn20i = new TH1D("shn20i","shn20i",100,0.0,0.3);
  TH1D *shn21r = new TH1D("shn21r","shn21r",100,0.0,0.3);
  TH1D *shn21i = new TH1D("shn21i","shn21i",100,0.0,0.3);
  TH1D *shn22r = new TH1D("shn22r","shn22r",100,0.0,0.3);
  TH1D *shn22i = new TH1D("shn22i","shn22i",100,0.0,0.3);
//   TH1D *shn30r = new TH1D("shn30r","shn30r",100,0.0,0.3);
//   TH1D *shn30i = new TH1D("shn30i","shn30i",100,0.0,0.3);
//   TH1D *shn31r = new TH1D("shn31r","shn31r",100,0.0,0.3);
//   TH1D *shn31i = new TH1D("shn31i","shn31i",100,0.0,0.3);
//   TH1D *shn32r = new TH1D("shn32r","shn32r",100,0.0,0.3);
//   TH1D *shn32i = new TH1D("shn32i","shn32i",100,0.0,0.3);
//   TH1D *shn33r = new TH1D("shn33r","shn33r",100,0.0,0.3);
//   TH1D *shn33i = new TH1D("shn33i","shn33i",100,0.0,0.3);

  TH1D *shd00r = new TH1D("shd00r","shd00r",100,0.0,0.3);
  TH1D *shd00i = new TH1D("shd00i","shd00i",100,0.0,0.3);
  TH1D *shd12r = new TH1D("shd12r","shd12r",100,0.0,0.3);
  TH1D *shd12i = new TH1D("shd12i","shd12i",100,0.0,0.3);
  TH1D *shd10r = new TH1D("shd10r","shd10r",100,0.0,0.3);
  TH1D *shd10i = new TH1D("shd10i","shd10i",100,0.0,0.3);
  TH1D *shd11r = new TH1D("shd11r","shd11r",100,0.0,0.3);
  TH1D *shd11i = new TH1D("shd11i","shd11i",100,0.0,0.3);
  TH1D *shd24r = new TH1D("shd24r","shd24r",100,0.0,0.3);
  TH1D *shd24i = new TH1D("shd24i","shd24i",100,0.0,0.3);
  TH1D *shd23r = new TH1D("shd23r","shd23r",100,0.0,0.3);
  TH1D *shd23i = new TH1D("shd23i","shd23i",100,0.0,0.3);
  TH1D *shd20r = new TH1D("shd20r","shd20r",100,0.0,0.3);
  TH1D *shd20i = new TH1D("shd20i","shd20i",100,0.0,0.3);
  TH1D *shd21r = new TH1D("shd21r","shd21r",100,0.0,0.3);
  TH1D *shd21i = new TH1D("shd21i","shd21i",100,0.0,0.3);
  TH1D *shd22r = new TH1D("shd22r","shd22r",100,0.0,0.3);
  TH1D *shd22i = new TH1D("shd22i","shd22i",100,0.0,0.3);
//   TH1D *shd30r = new TH1D("shd30r","shd30r",100,0.0,0.3);
//   TH1D *shd30i = new TH1D("shd30i","shd30i",100,0.0,0.3);
//   TH1D *shd31r = new TH1D("shd31r","shd31r",100,0.0,0.3);
//   TH1D *shd31i = new TH1D("shd31i","shd31i",100,0.0,0.3);
//   TH1D *shd32r = new TH1D("shd32r","shd32r",100,0.0,0.3);
//   TH1D *shd32i = new TH1D("shd32i","shd32i",100,0.0,0.3);
//   TH1D *shd33r = new TH1D("shd33r","shd33r",100,0.0,0.3);
//   TH1D *shd33i = new TH1D("shd33i","shd33i",100,0.0,0.3);

  TH1D *binctn = new TH1D("binctn","binctn",100,0.0,0.3);
  TH1D *binctd = new TH1D("binctd","binctd",100,0.0,0.3);

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

  TRandom2 *tRand = new TRandom2();
  tRand->SetSeed(33271);

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
  Double_t wgt;
  Double_t radx = 4.0 / 0.197327;
  Double_t rady = 6.0 / 0.197327;
  Double_t radz = 4.0 / 0.197327;

  radx *= radx;
  rady *= rady;
  radz *= radz;

  for (int iter=0; iter<1000000; iter++) {
    ko = (tRand->Rndm()*2.0 - 1.0) * 0.3 / sqrt(3.0);
    ks = (tRand->Rndm()*2.0 - 1.0) * 0.3 / sqrt(3.0);
    kl = (tRand->Rndm()*2.0 - 1.0) * 0.3 / sqrt(3.0);



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
    ylm00cr = SpherHarmonics::ReYlm(0, 0, ko, ks, kl);
    ylm00ci = SpherHarmonics::ImYlm(0, 0, ko, ks, kl);
    ylm11cr = SpherHarmonics::ReYlm(1, -1, ko, ks, kl);
    ylm11ci = SpherHarmonics::ImYlm(1, -1, ko, ks, kl);
    ylm10cr = SpherHarmonics::ReYlm(1, 0, ko, ks, kl);
    ylm10ci = SpherHarmonics::ImYlm(1, 0, ko, ks, kl);
    ylm11cr = SpherHarmonics::ReYlm(1, 1, ko, ks, kl);
    ylm11ci = SpherHarmonics::ImYlm(1, 1, ko, ks, kl);
    ylm24cr = SpherHarmonics::ReYlm(2, -2, ko, ks, kl);
    ylm24ci = SpherHarmonics::ImYlm(2, -2, ko, ks, kl);
    ylm23cr = SpherHarmonics::ReYlm(2, -1, ko, ks, kl);
    ylm23ci = SpherHarmonics::ImYlm(2, -1, ko, ks, kl);
    ylm20cr = SpherHarmonics::ReYlm(2, 0, ko, ks, kl);
    ylm20ci = SpherHarmonics::ImYlm(2, 0, ko, ks, kl);
    ylm21cr = SpherHarmonics::ReYlm(2, 1, ko, ks, kl);
    ylm21ci = SpherHarmonics::ImYlm(2, 1, ko, ks, kl);
    ylm22cr = SpherHarmonics::ReYlm(2, 2, ko, ks, kl);
    ylm22ci = SpherHarmonics::ImYlm(2, 2, ko, ks, kl);
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
    if (wgt > tRand->Rndm()*2.0) {
      shn00r->Fill(kv, ylm00cr);
      shn00i->Fill(kv, ylm00ci);
      shn12r->Fill(kv, ylm12cr);
      shn12i->Fill(kv, ylm12ci);
      shn10r->Fill(kv, ylm10cr);
      shn10i->Fill(kv, ylm10ci);
      shn11r->Fill(kv, ylm11cr);
      shn11i->Fill(kv, ylm11ci);
      shn24r->Fill(kv, ylm24cr);
      shn24i->Fill(kv, ylm24ci);
      shn23r->Fill(kv, ylm23cr);
      shn23i->Fill(kv, ylm23ci);
      shn20r->Fill(kv, ylm20cr);
      shn20i->Fill(kv, ylm20ci);
      shn21r->Fill(kv, ylm21cr);
      shn21i->Fill(kv, ylm21ci);
      shn22r->Fill(kv, ylm22cr);
      shn22i->Fill(kv, ylm22ci);
//       shn30r->Fill(kv, ylm30cr);
//       shn30i->Fill(kv, ylm30ci);
//       shn31r->Fill(kv, ylm31cr);
//       shn31i->Fill(kv, ylm31ci);
//       shn32r->Fill(kv, ylm32cr);
//       shn32i->Fill(kv, ylm32ci);
//       shn33r->Fill(kv, ylm33cr);
//       shn33i->Fill(kv, ylm33ci);
      binctn->Fill(kv, 1.0);
    }
    
    // Filling denominators
    shd00r->Fill(kv, ylm00cr);
    shd00i->Fill(kv, ylm00ci);
    shd12r->Fill(kv, ylm12cr);
    shd12i->Fill(kv, ylm12ci);
    shd10r->Fill(kv, ylm10cr);
    shd10i->Fill(kv, ylm10ci);
    shd11r->Fill(kv, ylm11cr);
    shd11i->Fill(kv, ylm11ci);
    shd24r->Fill(kv, ylm24cr);
    shd24i->Fill(kv, ylm24ci);
    shd23r->Fill(kv, ylm23cr);
    shd23i->Fill(kv, ylm23ci);
    shd20r->Fill(kv, ylm20cr);
    shd20i->Fill(kv, ylm20ci);
    shd21r->Fill(kv, ylm21cr);
    shd21i->Fill(kv, ylm21ci);
    shd22r->Fill(kv, ylm22cr);
    shd22i->Fill(kv, ylm22ci);
//    shd30r->Fill(kv, ylm30cr);
//     shd30i->Fill(kv, ylm30ci);
//     shd31r->Fill(kv, ylm31cr);
//     shd31i->Fill(kv, ylm31ci);
//     shd32r->Fill(kv, ylm32cr);
//     shd32i->Fill(kv, ylm32ci);
//     shd33r->Fill(kv, ylm33cr);
//     shd33i->Fill(kv, ylm33ci);
    binctd->Fill(kv, 1.0);
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

  long double tMq0[10];
  long double tTq0[10];
  long double tMTilde[10*10];
  long double tMTildeInv[10*10];
  long double tCq0[10];

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

  tTq0[0] = shn00r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[1] = shn12r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[2] = shn10r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[3] = shn11r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[4] = shn24r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[5] = shn23r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[6] = shn20r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[7] = shn21r->GetBinContent(20)/shn00r->GetEntries();
  tTq0[8] = shn22r->GetBinContent(20)/shn00r->GetEntries();
  //  tTq0[9] = shn33r->GetBinContent(20)/shn00r->GetEntries();

  cout << "20 bin of M " << endl;
  for (int iter=0; iter<10; iter++)
    cout << tMq0[iter] << "   ";
  cout << endl << endl;

  cout << "20 bin of T " << endl;
  for (int iter=0; iter<10; iter++)
    cout << tTq0[iter] << "   ";
  cout << endl << endl;

  MMatrix::GetMtilde(tMq0, tMTilde);
  
  // Liczymy wyznacznik
  long double det = 0;
//   long double psum;
//   for (int iters=0; iters<6; iters++) {
//     psum = 1.0;
//     for (int iterk=0; iterk<6; iterk++) {
//       psum *= tMTilde[((iters+iterk) % 6) * 6 + iterk];
//     }
//     det += psum;
//   }

//   for (int iters=0; iters<6; iters++) {
//     psum = 1.0;
//     for (int iterk=0; iterk<6; iterk++) {
//       psum *= tMTilde[((iters-iterk) % 6) * 6 + iterk];
//     }
//     det += psum;
//   }


  det = MMatrix::GetDet(tMTilde, MMatrix::GetMaxJM());

  cout << "Wyznacznik wynosi: " << det << endl;
  
  cout << "M matrix " << endl;
  for (int row=0; row<MMatrix::GetMaxJM(); row++) {
    for (int iter=0; iter<MMatrix::GetMaxJM(); iter++)
      cout << tMTilde[row*MMatrix::GetMaxJM()+iter] << " \t";
    cout << endl << endl;
  }

  
  for (int row=0; row<MMatrix::GetMaxJM(); row++) {
    for (int iter=0; iter<MMatrix::GetMaxJM(); iter++)
      tMTildeInv[iter*MMatrix::GetMaxJM()+row] = MMatrix::GetDetOne(tMTilde, MMatrix::GetMaxJM(), row, iter)/det;
  }

  //  MMatrix::InvertMTilde(tMTilde, tMTildeInv);

  cout << "M matrix inverted " << endl;
  for (int row=0; row<MMatrix::GetMaxJM(); row++) {
    for (int iter=0; iter<MMatrix::GetMaxJM(); iter++)
      cout << tMTildeInv[row*MMatrix::GetMaxJM()+iter] << "\t";
    cout << endl << endl;
  }

  long double tIdent[10*10];
  for (int irow=0; irow<MMatrix::GetMaxJM(); irow++) {
    for (int icol=0; icol<MMatrix::GetMaxJM(); icol++) {
      long double sum = 0;
      for (int iter=0; iter<MMatrix::GetMaxJM(); iter++) {
	sum += tMTilde[irow*MMatrix::GetMaxJM() + iter]*tMTildeInv[iter*MMatrix::GetMaxJM()+icol];
      }
      tIdent[irow*MMatrix::GetMaxJM()+icol] = sum;
    }
  }

  cout << "M*M^-1 " << endl;
  for (int row=0; row<MMatrix::GetMaxJM(); row++) {
    for (int iter=0; iter<MMatrix::GetMaxJM(); iter++)
      cout << tIdent[row*MMatrix::GetMaxJM()+iter] << "\t";
    cout << endl << endl;
  }

  for (int ilmc=0; ilmc<MMatrix::GetMaxJM(); ilmc++) {
    tCq0[ilmc] = 0.0;
    for (int ilmt=0; ilmt<MMatrix::GetMaxJM(); ilmt++) {
      tCq0[ilmc] += tMTildeInv[ilmc*MMatrix::GetMaxJM() + ilmt] * tTq0[ilmt];
    }
    //    tCq0[ilmc] /= (shd00r->GetEntries()*shn00r->GetEntries());
  }
    
  cout << "20 bin of C " << endl;
  for (int iter=0; iter<MMatrix::GetMaxJM(); iter++)
    cout << tCq0[iter] << "   ";
  cout << endl << endl;
  
}
  
