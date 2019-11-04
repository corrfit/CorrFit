#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
  double am[9] = { 0.155, 0.8331, 1.994,
		   2.558, 1.922,  0.992,
		   4.88,  1.9920, 0.322 };
  
  double mu[9] = { 0.155, 0.8331, 1.994,
		   2.558, 1.922,  0.992,
		   4.88,  1.9920, 0.322 };
  
  double mv[9] = { 0.0, 0.0, 0.0,
		   0.0, 0.0, 0.0,
		   0.0, 0.0, 0.0 };
  
  double vs[3] = { 0.0, 0.0, 0.0 };
  double vw[3] = { 0.0, 0.0, 0.0 };

  double vb[3] = {92.02, 42.01, 63.99};
  double vx[3] = {0.0, 0.0, 0.0};
  double vy[3] = {0.0, 0.0, 0.0};

  gsl_matrix_view matA = gsl_matrix_view_array(am, 3, 3);
  gsl_matrix_view matV = gsl_matrix_view_array(mv, 3, 3);
  gsl_matrix_view matU = gsl_matrix_view_array(mu, 3, 3);

  gsl_vector_view vecS = gsl_vector_view_array(vs, 3);
  gsl_vector_view vecW = gsl_vector_view_array(vw, 3);
  
  gsl_vector_view vecB = gsl_vector_view_array(vb, 3);
  gsl_vector_view vecX = gsl_vector_view_array(vx, 3);

  gsl_vector_view vecY = gsl_vector_view_array(vy, 3);

  cout << "Decomposed " << endl;
  cout << am[0] << " " << am[1] << " " << am[2] << endl;
  cout << am[3] << " " << am[4] << " " << am[5] << endl;
  cout << am[6] << " " << am[7] << " " << am[8] << endl;

  gsl_linalg_SV_decomp(&matU.matrix, &matV.matrix, &vecS.vector, &vecW.vector);
  gsl_linalg_SV_solve( &matU.matrix, &matV.matrix, &vecS.vector, &vecB.vector, &vecX.vector);

  cout << "Into U " << endl;
  cout << mu[0] << " " << mu[1] << " " << mu[2] << endl;
  cout << mu[3] << " " << mu[4] << " " << mu[5] << endl;
  cout << mu[6] << " " << mu[7] << " " << mu[8] << endl;
  
  cout << "V" << endl;
  cout << mv[0] << " " << mv[1] << " " << mv[2] << endl;
  cout << mv[3] << " " << mv[4] << " " << mv[5] << endl;
  cout << mv[6] << " " << mv[7] << " " << mv[8] << endl;

  cout << "S" << endl;
  cout << vs[0] << " " << vs[1] << " " << vs[2] << endl;

  cout << "W" << endl;
  cout << vw[0] << " " << vw[1] << " " << vw[2] << endl;
  
  cout << "Solved Ax=b" << endl;

  cout << "B" << endl;
  cout << vb[0] << " " << vb[1] << " " << vb[2] << endl;
  
  cout << "X" << endl;
  cout << vx[0] << " " << vx[1] << " " << vx[2] << endl;

  gsl_blas_dgemv(CblasNoTrans, 1.0, &matA.matrix, &vecX.vector, 1.0, &vecY.vector);

  cout << "Y" << endl;
  cout << vy[0] << " " << vy[1] << " " << vy[2] << endl;
  
}
