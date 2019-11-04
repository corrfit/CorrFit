#include <math.h>
#include <iostream>
#include "sf.h"
#include <TMath.h>
#include <Math/SpecFunc.h>

using namespace std;

int main(int argc, char **argv)
{
  double x = 0.2;
  double y = 0.5;
  double z = 0.1;

  if (argc > 1) {
    x = atof(argv[1]);
    y = atof(argv[2]);
    z = atof(argv[3]);
  }

  complex<double> ret00 = SpherHarmonics::Ylm(0, 0, x, y, z);
  complex<double> ret10 = SpherHarmonics::Ylm(1, 0, x, y, z);
  complex<double> ret11 = SpherHarmonics::Ylm(1, 1, x, y, z);
  complex<double> ret12 = SpherHarmonics::Ylm(1, -1, x, y, z);
  complex<double> ret20 = SpherHarmonics::Ylm(2, 0, x, y, z);
  complex<double> ret21 = SpherHarmonics::Ylm(2, 1, x, y, z);
  complex<double> ret22 = SpherHarmonics::Ylm(2, 2, x, y, z);
  complex<double> ret23 = SpherHarmonics::Ylm(2, -1, x, y, z);
  complex<double> ret24 = SpherHarmonics::Ylm(2, -2, x, y, z);

  complex<double> ylms[18];
  SpherHarmonics::YlmUpToL(2, x, y, z, ylms);

  cout << "theta phi " << z/sqrt(x*x+y*y+z*z) << " " << atan2(y,x) << endl;

  cout << "leg " << gsl_sf_legendre_sphPlm(2, 2, z/sqrt(x*x+y*y+z*z)) << endl;
  cout << "leg " << gsl_sf_legendre_sphPlm(2, 0, -z/sqrt(x*x+y*y+z*z)) << endl;

  double legpar = TMath::Sqrt((2*2+1)/(4*TMath::Pi()));
  double legsil = TMath::Sqrt(1.0/(4*3*2));

  cout << "leg " << legpar*legsil*ROOT::Math::assoc_legendre(2, 2, z/sqrt(x*x+y*y+z*z)) << endl;
  cout << "leg " << legpar*ROOT::Math::assoc_legendre(2, 0, -z/sqrt(x*x+y*y+z*z)) << endl;

  cout << "Ylm 0  0 is " << real(ret00) << " " << imag(ret00) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 1 -1 is " << real(ret12) << " " << imag(ret12) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 1  0 is " << real(ret10) << " " << imag(ret10) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 1  1 is " << real(ret11) << " " << imag(ret11) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2 -2 is " << real(ret24) << " " << imag(ret24) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2 -1 is " << real(ret23) << " " << imag(ret23) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2  0 is " << real(ret20) << " " << imag(ret20) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2  1 is " << real(ret21) << " " << imag(ret21) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2  2 is " << real(ret22) << " " << imag(ret22) << " for " << x << " " << y << " " << z << endl;

  cout << endl;

  cout << "Ylm 0  0 is " << real(ylms[0]) << " " << imag(ylms[0]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 1 -1 is " << real(ylms[1]) << " " << imag(ylms[1]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 1  0 is " << real(ylms[2]) << " " << imag(ylms[2]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 1  1 is " << real(ylms[3]) << " " << imag(ylms[3]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2 -2 is " << real(ylms[4]) << " " << imag(ylms[4]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2 -1 is " << real(ylms[5]) << " " << imag(ylms[5]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2  0 is " << real(ylms[6]) << " " << imag(ylms[6]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2  1 is " << real(ylms[7]) << " " << imag(ylms[7]) << " for " << x << " " << y << " " << z << endl;
  cout << "Ylm 2  2 is " << real(ylms[8]) << " " << imag(ylms[8]) << " for " << x << " " << y << " " << z << endl;

  
}
