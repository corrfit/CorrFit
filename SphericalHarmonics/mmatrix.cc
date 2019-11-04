#include "mmatrix.h"
#include <TMatrix.h>
#include <TMath.h>

const double factorials[11] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0 };
// const int    maxjm = 10; // We cut the decomposition after 2
// const double els[10] = { 0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0 };
// const double ems[10] = { 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 3.0 };
// const int    elsi[10] = { 0, 1, 1, 2, 2, 2, 3, 3, 3, 3 };
// const int    emsi[10] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };
const int    maxjm = 9; // We cut the decomposition after 2
const double els[10] = { 0.0,  1.0, 1.0, 1.0,  2.0,  2.0, 2.0, 2.0, 2.0, 3.0 };
const double ems[10] = { 0.0, -1.0, 0.0, 1.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0 };
const int    elsi[10] = { 0, 1,  1, 1,  2,  2, 2, 2, 2, 3 };
const int    emsi[10] = { 0, -1, 0, 1, -2, -1, 0, 1, 2, 3 };

double MMatrix::ClebschGordan(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm)
{
  int mint, maxt;
  double cgc = 0.0;
  int titer;
  double coef;

  maxt = lrint(aJot1 + aJot2 - aJot);
  mint = 0;
  if (lrint(aJot1 - aEm1) < maxt) maxt = lrint(aJot1 - aEm1);
  if (lrint(aJot2 + aEm2) < maxt) maxt = lrint(aJot2 + aEm2);
  if (lrint(-(aJot-aJot2+aEm1)) > mint) mint = lrint(-(aJot-aJot2+aEm1));
  if (lrint(-(aJot-aJot1-aEm2)) > mint) mint = lrint(-(aJot-aJot1-aEm2));

  for (titer = mint; titer<=maxt; titer ++)
    {
      coef = TMath::Power(-1, titer);
      coef *= TMath::Sqrt((2*aJot+1)*
			  factorials[lrint(aJot1+aEm1)] *
			  factorials[lrint(aJot1-aEm1)] *
			  factorials[lrint(aJot2+aEm2)] *
			  factorials[lrint(aJot2-aEm2)] *
			  factorials[lrint(aJot+aEm)] *
			  factorials[lrint(aJot-aEm)]);
      coef /= (factorials[titer] *
	       factorials[lrint(aJot1+aJot2-aJot-titer)] *
	       factorials[lrint(aJot1-aEm1-titer)] *
	       factorials[lrint(aJot2+aEm2-titer)] *
	       factorials[lrint(aJot-aJot2+aEm1+titer)] *
	       factorials[lrint(aJot-aJot1-aEm2+titer)]);
      
      cgc += coef;
    }

  cgc *= DeltaJ(aJot1, aJot2, aJot);

  return cgc;
}

double MMatrix::DeltaJ(double aJot1, double aJot2, double aJot)
{
  double res = TMath::Sqrt(1.0 * 
			   factorials[lrint(aJot1+aJot2-aJot)] * 
			   factorials[lrint(aJot1-aJot2+aJot)] * 
			   factorials[lrint(-aJot1+aJot2+aJot)] / 
			   factorials[lrint(aJot1+aJot2+aJot+1)]);
  
  return res;
}

double MMatrix::WignerSymbol(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm)
{
  if (lrint(aEm1+aEm2+aEm) != 0.0) 
    return 0.0;
  double cge = ClebschGordan(aJot1, aEm1, aJot2, aEm2, aJot, -aEm);
  if (lrint(abs(aJot1 - aJot2 - aEm)) % 2) 
    cge *= -1.0;
  cge /= sqrt(2*aJot + 1);

  if (cge == -0.0) cge = 0.0;

  return cge;
}


void MMatrix::GetMtilde(complex<double> *aMat, double *aMTilde)
{
  // Create the Mtilde for a given q bin
  double lzero, mzero;
  double lprim, mprim;
  double lbis, mbis;
 
  int lzeroi, mzeroi;
  int lprimi, mprimi;
  int lbisi, mbisi;

  for (int izero = 0; izero<GetMaxJM(); izero++) {
    GetElEmForIndex(izero, &lzero, &mzero);
    GetElEmForIndex(izero, &lzeroi, &mzeroi);
    for (int ibis = 0; ibis<GetMaxJM(); ibis++) {
      GetElEmForIndex(ibis, &lbis, &mbis);
      GetElEmForIndex(ibis, &lbisi, &mbisi);
      complex<double> val = complex<double>(0.0, 0.0);
      complex<double> mcomp[maxjm];
      for (int iprim = 0; iprim<GetMaxJM(); iprim++) {

	GetElEmForIndex(iprim, &lprim, &mprim);
	GetElEmForIndex(iprim, &lprimi, &mprimi);

	if (abs(mzeroi) % 2) mcomp[iprim] = complex<double>(-1.0, 0.0); // (-1)^m
	else mcomp[iprim] = complex<double>(1.0, 0.0);
	
	mcomp[iprim] *= sqrt((2*lzero+1)*(2*lprim+1)*(2*lbis+1));   // P1
	mcomp[iprim] *= WignerSymbol(lzero, 0, lprim, 0, lbis, 0); // W1
	mcomp[iprim] *= WignerSymbol(lzero, -mzero, lprim, mprim, lbis, mbis); // W2
	mcomp[iprim] *= aMat[iprim];
	val += mcomp[iprim];
      }
      aMTilde[(izero*2)*(2*GetMaxJM()) + (ibis*2)]     =  real(val);
      aMTilde[(izero*2+1)*(2*GetMaxJM()) + (ibis*2)]   =  imag(val);
      if (imag(val) != 0.0)
	aMTilde[(izero*2)*(2*GetMaxJM()) + (ibis*2+1)]   = -imag(val);
      else 
	aMTilde[(izero*2)*(2*GetMaxJM()) + (ibis*2+1)]   = 0.0;
      aMTilde[(izero*2+1)*(2*GetMaxJM()) + (ibis*2+1)] =  real(val);
      
    }
  }
}

void MMatrix::InvertMTilde(double *aMTilde, double *mTildeInv)
{
//   TMatrixT<double> tmat(GetMaxJM(), GetMaxJM());
//   for (int iterx=0; iterx<6; iterx++) 
//     for (int itery=0; itery<6; itery++)
//       tmat[iterx][itery] = aMTilde[iterx*6+itery];
//   tmat.Print();
//   TMatrixT<double> tmatinv = tmat.Invert();
//   for (int iter=0; iter<GetMaxJM()*GetMaxJM(); iter++) 
//     {
//       mTildeInv[iter] = tmatinv.GetMatrixArray()[iter];
//     }
  double det = GetDet(aMTilde, GetMaxJM());

  for (int row=0; row<GetMaxJM(); row++) {
    for (int iter=0; iter<GetMaxJM(); iter++)
      mTildeInv[iter+GetMaxJM()*row] = GetDetOne(aMTilde, GetMaxJM(), row, iter)/det;
  }
}

int  MMatrix::GetMaxJM() 
{ return maxjm; }

void MMatrix::GetElEmForIndex(int aIndex, double *aEl, double *aEm)
{
  *aEl = els[aIndex];
  *aEm = ems[aIndex];
}

void MMatrix::GetElEmForIndex(int aIndex, int *aEl, int *aEm)
{
  *aEl = elsi[aIndex];
  *aEm = emsi[aIndex];
}

double MMatrix::GetDet(double *aMatrix, int size)
{
  double det = 0;
  if (size >3) {
    double *mat = (double *)malloc(sizeof(double)*(size-1)*(size-1));
    double detcmp = 0;
    for (int iter=0; iter<size; iter++) {
      //      detcmp = aMatrix[iter];

      // Rewrite the matrix
      int inrow, inncol;
      for (int irow=1, inrow=0; irow<size; irow++, inrow++) {
	for (int icol=0, incol=0; icol<size; icol++, incol++) {
	  if ((icol == iter) && (iter == 0)) icol++;
	  mat[inrow*(size-1) + incol] = aMatrix[irow*size + icol]; 
	  if (icol == iter-1) icol++;
	}
      }
      
      detcmp *= GetDet(mat, size-1);
      if ((detcmp != 0.0) && (iter % 2)) detcmp *= -1; 
      
      det += detcmp;
    }
    free(mat);
  }
  else {
    if (size == 1) 
      //      return aMatrix[0];
    if (size == 2) {
//       det = (aMatrix[0]*aMatrix[3] -
// 	     aMatrix[1]*aMatrix[2]);
    }
    if (size == 3) {
//       det = (aMatrix[0]*aMatrix[4]*aMatrix[8] +
// 	     aMatrix[1]*aMatrix[5]*aMatrix[6] +
// 	     aMatrix[2]*aMatrix[3]*aMatrix[7] -
// 	     aMatrix[0]*aMatrix[5]*aMatrix[7] -
// 	     aMatrix[1]*aMatrix[3]*aMatrix[8] -
// 	     aMatrix[2]*aMatrix[4]*aMatrix[6]);
    }
  }
  return det;
}

double MMatrix::GetDetOne(double *aMatrix, int size, int rowx, int rowy)
{

  double *mat = (double *)malloc(sizeof(double)*(size-1)*(size-1));
  // Rewrite the matrix
  for (int irow=0, inrow=0; irow<size; irow++, inrow++) {
    if ((irow == rowy) && (rowy == 0)) irow++;
    for (int icol=0, incol=0; icol<size; icol++, incol++) {
      if ((icol == rowx) && (rowx == 0)) icol++;
      mat[inrow*(size-1) + incol] = aMatrix[irow*size + icol]; 
      if (icol == rowx-1) icol++;
    }
    if (irow == rowy-1) irow++;
  }
  
  double det = GetDet(mat, size-1);
  
  free(mat);    
  
  if (((rowx+rowy) % 2) && (det !=0)) det *= -1;

  return det;
}
