#include <math.h>
#include <iostream>

#define SIZE 4

using namespace std;

long double getdet(long double *aMatrix, int size)
{
  cout << "Getting det for "<<endl;
  
  for (int row=0; row<size; row++) {
    for (int iter=0; iter<size; iter++) {
      cout << aMatrix[row*size+iter] << " \t";
    }
    cout << endl;
  }


  long double det = 0;
  if (size >3) {
    long double *mat = (long double *)malloc(sizeof(long double)*(size-1)*(size-1));
    long double detcmp = 0;
    for (int iter=0; iter<size; iter++) {
      detcmp = aMatrix[iter];

      // Rewrite the matrix
      int inrow, inncol;
      for (int irow=1, inrow=0; irow<size; irow++, inrow++) {
	for (int icol=0, incol=0; icol<size; icol++, incol++) {
	  if (icol == iter) icol++;
	  mat[inrow*(size-1) + incol] = aMatrix[irow*size + icol]; 
	}
      }
      
      detcmp *= getdet(mat, size-1);
      if ((detcmp != 0.0) && (iter % 2)) detcmp *= -1; 
      
      det += detcmp;
    }
    //    free(mat);
  }
  else {
    if (size == 1) 
      return aMatrix[0];
    if (size == 2) {
      det = (aMatrix[0]*aMatrix[3] -
	     aMatrix[1]*aMatrix[2]);
    }
    if (size == 3) {
      det = (aMatrix[0]*aMatrix[4]*aMatrix[8] +
	     aMatrix[1]*aMatrix[5]*aMatrix[6] +
	     aMatrix[2]*aMatrix[3]*aMatrix[7] -
	     aMatrix[0]*aMatrix[5]*aMatrix[7] -
	     aMatrix[1]*aMatrix[3]*aMatrix[8] -
	     aMatrix[2]*aMatrix[4]*aMatrix[6]);
    }
  }
  return det;
}

long double getdetone(long double *aMatrix, int size, int rowx, int rowy)
{
  cout << "Getting det for " << rowx << " " <<  rowy <<endl;
  
  for (int row=0; row<size; row++) {
    if (row==rowy) continue;
    for (int iter=0; iter<size; iter++) {
      if (iter==rowx) continue;
      cout << aMatrix[row*size+iter] << " \t";
    }
    cout << endl;
  }

  long double *mat = (long double *)malloc(sizeof(long double)*(size-1)*(size-1));
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
  
  long double det = getdet(mat, size-1);
  
  //  free(mat);    
  
  if (((rowx+rowy) % 2) && (det !=0)) det *= -1;

  cout << "Det is " << det << endl;
  
  return det;
}

int main(int argc, char** argv)
{
  
  long double tMat[SIZE*SIZE];
  long double tMatInv[SIZE*SIZE];
  
  tMat[0]  = 1;
  tMat[1]  = 0;
  tMat[2]  = 0;
  tMat[3]  = 5;
  tMat[4]  = 0;
  tMat[5]  = 1;
  tMat[6]  = 0;
  tMat[7]  = 0;
  tMat[8]  = 12;
  tMat[9]  = 0;
  tMat[10] = 1;
  tMat[11] = 0;
  tMat[12] = 16;
  tMat[13] = 3;
  tMat[14] = 2;
  tMat[15] = 1;

  // Liczymy wyznacznik
//   long double det = 0;
//   long double psum;
//   for (int iters=0; iters<SIZE; iters++) {
//     psum = 1.0;
//     for (int iterk=0; iterk<SIZE; iterk++) {
//       psum *= tMat[iterk * SIZE + ((iters + iterk) % 6)];
//     }
//     cout << "psum is " << psum << endl;
//     det += psum;
//   }

//   for (int iters=0; iters<SIZE; iters++) {
//     psum = 1.0;
//     for (int iterk=0; iterk<SIZE; iterk++) {
//       int colx = iters-iterk + SIZE*((iters-iterk)<0);
//       psum *= tMat[(colx % SIZE) + iterk * SIZE];
//     }
//     cout << "psum is " << psum << endl;
//     det -= psum;
//   }

  long double det = getdet(tMat, 4);

  cout << "Wyznacznik wynosi: " << det << endl;

  cout << "M matrix " << endl;
  for (int row=0; row<SIZE; row++) {
    for (int iter=0; iter<SIZE; iter++)
      cout << tMat[row*SIZE+iter] << " \t";
    cout << endl << endl;
  }

  
  for (int row=0; row<SIZE; row++) {
    for (int iter=0; iter<SIZE; iter++)
      tMatInv[iter+SIZE*row] = getdetone(tMat, SIZE, row, iter)/det;
  }

  //  MMatrix::InvertMat(tMat, tMatInv);

  cout << "M matrix inverted " << endl;
  for (int row=0; row<SIZE; row++) {
    for (int iter=0; iter<SIZE; iter++)
      cout << tMatInv[row*SIZE+iter] << "\t";
    cout << endl << endl;
  }

  long double tIdent[SIZE*SIZE];
  for (int irow=0; irow<SIZE; irow++) {
    for (int icol=0; icol<SIZE; icol++) {
      long double sum = 0;
      for (int iter=0; iter<SIZE; iter++) {
	sum += tMat[irow*SIZE + iter]*tMatInv[iter*SIZE+icol];
      }
      tIdent[irow*SIZE+icol] = sum;
    }
  }

  cout << "M*M^-1 " << endl;
  for (int row=0; row<SIZE; row++) {
    for (int iter=0; iter<SIZE; iter++)
      cout << tIdent[row*SIZE+iter] << "\t";
    cout << endl << endl;
  }

  
}
