#ifndef ALIFEMTOCORRFCTNDIRECTYLM_H
#define ALIFEMTOCORRFCTNDIRECTYLM_H

// Correlation function that is binned in Ylms directly
// Provides a way to store the numerator and denominator
// in Ylms directly and correctly calculate the correlation
// function from them

#include <math.h>
#include <complex>
#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
#include "StHbtCorrFctn.hh"
#include "StHbtYlm.h"

using namespace std;
class TRootIOCtor;

class StHbtCorrFctnDirectYlm: public StHbtCorrFctn {
 public:
  StHbtCorrFctnDirectYlm();
  StHbtCorrFctnDirectYlm(const char *name, int maxl, int ibin, double vmin, double vmax);
  ~StHbtCorrFctnDirectYlm();

  void AddRealPair(double *qvec, double weight=1.0);
  void AddMixedPair(double *qvec, double weight=1.0);

  void AddRealPair(double qout, double qside, double qlong, double weight=1.0);
  void AddMixedPair(double qout, double qside, double qlong, double weight=1.0);

  virtual StHbtString Report();

  virtual void AddRealPair(const StHbtPair* aPair);
  virtual void AddMixedPair(const StHbtPair* aPair);

  void SetR2Factor(bool aFactor);

  virtual void Finish();

  void Write();

  void ReadFromFile(TFile *infile, const char *name, int maxl);

  TH1D *GetNumRealHist(int el, int em);
  TH1D *GetNumImagHist(int el, int em);

  TH1D *GetDenRealHist(int el, int em);
  TH1D *GetDenImagHist(int el, int em);

 private:
  double ClebschGordan(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm);
  double DeltaJ(double aJot1, double aJot2, double aJot);
  double WignerSymbol(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm);
  
  void GetMtilde(complex<double>* aMat, double *aMTilde); 
  
  int  GetMaxJM();
  void GetElEmForIndex(int aIndex, double *aEl, double *aEm);
  void GetElEmForIndex(int aIndex, int *aEl, int *aEm);
  int  GetBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag);

  int  PackYlmVector(double *invec, double *outvec);
  int  PackYlmMatrix(double *inmat, double *outmat);

  int GetIndexForLM(int el, int em);

  void PackCovariances();
  void UnpackCovariances();

  TH1D **numsreal;            //! Real parts of Ylm components of the numerator
  TH1D **numsimag;            //! Imaginary parts of Ylm components of the numerator
  TH1D **densreal;            //! Real parts of Ylm components of the denominator
  TH1D **densimag;            //! Imaginary parts of Ylm components of the denominator

  TH1D **cfctreal;            //! Real parts of Ylm components of the correlation function
  TH1D **cfctimag;	      //! Imaginary parts of Ylm components of the correlation function
  
  TH1D *binctn;               //! Bin occupation for the numerator
  TH1D *binctd;               //! Bin occupation for the denominator

  TH3D *covnum;               //! Numerator covariance matrix packed into TH3D
  TH3D *covden;               //! Denominator covariance matrix packed into TH3D

  double *covmnum;            //! Covariance matrix for the numerator
  double *covmden;            //! Covariance matrix for the denominator
  double *covmcfc;            //! Covariance matrix for the correlation function

  int fMaxL;                  // l cut-off of the decomposition

  int    maxjm;               // number of l-m combinations
  double *els;                //! table of l's
  double *ems;                //! table of m's
  int    *elsi;               //! table of integer l's
  int    *emsi;               //! table of integer m's

  complex<double> *ylmbuffer; //! buffer for ylm calculation
  double *factorials;         //! Helper table of factorials

  int fR2factor;              // If 1 all entries are scaled by 1/r^2
#ifdef __ROOT__ 
  ClassDef(StHbtCorrFctnDirectYlm, 1)
#endif

};

#endif

