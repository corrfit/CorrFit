#include <iostream>
#include "CorrFctnDirectYlm.h"
#include <TH1D.h>
#include <TRandom2.h>
#include <math.h>
#include <TFile.h>
#include <complex>

using namespace std;

int main(int argc, char **argv)
{
  TFile *ofile = new TFile("shmout.root","RECREATE");

  CorrFctnDirectYlm *cylm = new CorrFctnDirectYlm("PionCYlm",2,50,0.0,0.3);

  TRandom2 *tRand = new TRandom2();
  tRand->SetSeed(136571);

  Double_t kvec[3];
  Double_t wgt;
  Double_t radx = 6.0 / 0.197327;
  Double_t rady = 4.0 / 0.197327;
  Double_t radz = 7.0 / 0.197327;

  radx *= radx;
  rady *= rady;
  radz *= radz;

  int flip = 0;

  for (int iter=0; iter<10000000; iter++) {
    if (!((iter-1) % 100000))
      cout << "Pair " << iter<< endl;
    
    kvec[0] = (tRand->Rndm()*2.0 - 1.0) * 0.32 / sqrt(3.0);
    kvec[1] = (tRand->Rndm()*2.0 - 1.0) * 0.32 / sqrt(3.0);
    kvec[2] = (tRand->Rndm()*2.0 - 1.0) * 0.32 / sqrt(3.0);

    wgt = 1.0 + exp(-kvec[0]*kvec[0]*radx
		    -kvec[1]*kvec[1]*rady
		    -kvec[2]*kvec[2]*radz);

    // Filling the numerators
    if (flip) {
//       if (wgt > tRand->Rndm()*2.0) {
//	cylm->AddRealPair(kvec, 1.0);
//       }
      cylm->AddRealPair(kvec, wgt);
      flip = 0;
    }
    else {
      // Filling denominators
      cylm->AddMixedPair(kvec, 1.0);
      flip = 1;
    }
  }
  
  cylm->Finish();

  cylm->Write();
}

  
