#ifndef _MMATRIX_H_
#define _MMATRIX_H_

#include <math.h>
#include <complex>

using namespace std;

namespace MMatrix {
  double ClebschGordan(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm);
  double DeltaJ(double aJot1, double aJot2, double aJot);
  double WignerSymbol(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm);

  void GetMtilde(complex<double>* aMat, double *aMTilde); 
  void InvertMTilde(double *aMTilde, double *mTildeInv);
  
  int  GetMaxJM();
  void GetElEmForIndex(int aIndex, double *aEl, double *aEm);
  void GetElEmForIndex(int aIndex, int *aEl, int *aEm);

  double GetDet(double *aMatrix, int size);
  double GetDetOne(double *aMatrix, int size, int rowx, int rowy);

};

#endif

