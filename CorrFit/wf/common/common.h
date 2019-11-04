#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstring>
using namespace std;

class Tplanewave{
 public:
  double PI;
  complex<double>ci;
  int q1q2;
  complex<double> couly;
  complex<double> cfact1;
  complex<double> cfact2;
  complex<double> chype[41];
  complex<double> chype1[6];
  complex<double> chype2[6];
  double q,eta; // q is the reduced mom., (p1-p2)/2
  Tplanewave(double eta,int Q1Q2,double qset);
  complex<double> planewave(double r,double ctheta);
  complex<double> hyper(complex<double> a,complex<double> b,
			complex<double> cz);
};

class Tpartwave{
public:
  complex<double>ci;
  double PI;
  int ell,q1q2;
  double q,eta,sigma;
  Tpartwave(double etaset,int q1q2,double qset,int ell);
  complex<double> *phi; /* An incoming Coulomb partial wave
			   which behaves as e^{-i(kr-eta*ln(2kr)+sigma)} */ 
  int nxmax;
  double delx;
  void phi_init();
  complex<double> GetPhiIncoming(double r);
};

class Twavefunction{
 public:
  double PI;
  double ROOT2;
  complex<double> ci;
  double HBARC;
  double MPI,MKAON,MPROTON,MLAMBDA;
  int *ell;
  double **Wepsilon,**delta,**ddeltadq,*eta,*channelweight,*qarray;
  int nqmax,ellmax,nchannels,q1q2;
  double m1,m2,delq,epsilon;
  bool COULOMB,STRONG,generic;
  double mu,muscale,symmweight;
  int q1q2scale;
  Tplanewave **planewave;
  Tpartwave ***partwave;
  // These last two are redefined (overwritten) in inherited classes
  void InitArrays();
  void InitWaves();
  void printCdelta(double Rx,double Ry,double Rz);
  void ParsInit(char *parsfilename);
  double getpsisquared(double *pa,double *xa,double *pb,double *xb);
  double getpsisquared(double q,double r,double ctheta);
  virtual double calcpsisquared(int iq,double r,double ctheta);
  void getqrctheta(double *pa,double *xa,double *pb,double *xb,
		   double *q,double *r,double *ctheta);

  Twavefunction();
};

bool comparestrings(char *s1,char *s2);
double legendre(int ell,double ctheta);
double dgamma(int mm);
double triangle(double m0,double m1,double m2);
complex<double> cgamma(complex<double> c);
complex<double> CWincoming_bigr(int ell,double x,double eta);
complex<double> CWoutgoing_bigr(int ell,double x,double eta);
complex<double> CWoutgoing_smallr(int ell,double x,double eta);
complex<double> CWincoming_smallr(int ell,double x,double eta);
complex<double> hankel(int ell,double x);
complex<double> hankelstar(int ell,double x);
double GetIW(int ell,double epsilon,double q,int q1q2,double eta,
	     double delta);
void phaseshift_CoulombCorrect(int ell,double q,double eta,
			       double &delta,double &ddeltadq);
void getphaseshift_pipi(int I,int ell,double q,double *delta,
			double *ddeltadq);
