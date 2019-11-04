#include "StHbtCorrFctnDirectYlm.h"
#include <TMath.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <iostream>

using namespace std;
#define MAXJM 25

#ifdef __ROOT__ 
ClassImp(StHbtCorrFctnDirectYlm)
#endif

StHbtCorrFctnDirectYlm::StHbtCorrFctnDirectYlm(const char *name, int maxl, int ibin=30, double vmin=0.0, double vmax=0.3):
  numsreal(0),
  numsimag(0),
  densreal(0),
  densimag(0),
  cfctreal(0),
  cfctimag(0),
  binctn(0),
  binctd(0),
  covnum(0),
  covden(0),
  covmnum(0),
  covmden(0),
  covmcfc(0),
  els(0),
  ems(0),
  elsi(0),
  emsi(0),
  ylmbuffer(0),
  factorials(0),
  fR2factor(0)
{
  fMaxL = maxl;
  maxjm = (maxl+1)*(maxl+1);

  cout <<  "Size is " << sizeof(double) << " " << sizeof(complex<double>) << endl;

  // Fill in factorials table
  factorials = (double *) malloc(sizeof(double) * (4 * (maxl + 1)));
  int fac = 1;
  factorials[0] = 1;
  for (int iter=1; iter<4*(maxl+1); iter++)
    {
      fac *= iter;
      factorials[iter] = fac;
    }

  // Fill in els and ems table
  int el = 0;
  int em = 0;
  int il = 0;
  els = (double *) malloc(sizeof(double) * (maxjm));
  ems = (double *) malloc(sizeof(double) * (maxjm));
  elsi = (int *) malloc(sizeof(int) * (maxjm));
  emsi = (int *) malloc(sizeof(int) * (maxjm));
  do {
    els[il] = el;
    ems[il] = em;
    elsi[il] = (int) el;
    emsi[il] = (int) em;

    cout << "il el em " << il << " " << elsi[il] << " " << emsi[il] << endl;
    em++;
    il++;
    if (em > el) {
      el++;
      em = -el;
    }
  }
  while (el <= maxl);
  
  for (il=0; il<maxjm; il++)
    cout << "il el em " << il << " " << elsi[il] << " " << emsi[il] << endl;

  // Create numerator and denominator historgrams
  //  int sthp = sizeof(TH1D *);
  //  numsreal = (TH1D **) malloc(sthp * maxjm);
//   numsreal = new TH1D * [maxjm];
//   numsimag = new TH1D * [maxjm];
//   densreal = new TH1D * [maxjm];
//   densimag = new TH1D * [maxjm];
//   cfctreal = new TH1D * [maxjm];
//   cfctimag = new TH1D * [maxjm];
  numsreal = (TH1D **) malloc(sizeof(TH1D *) * maxjm);
  numsimag = (TH1D **) malloc(sizeof(TH1D *) * maxjm);
  densreal = (TH1D **) malloc(sizeof(TH1D *) * maxjm);
  densimag = (TH1D **) malloc(sizeof(TH1D *) * maxjm);
  cfctreal = (TH1D **) malloc(sizeof(TH1D *) * maxjm);
  cfctimag = (TH1D **) malloc(sizeof(TH1D *) * maxjm);
  
  char bufname[200];
  for (int ihist=0; ihist<maxjm; ihist++) {
    sprintf(bufname, "NumReYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    numsreal[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "NumImYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    numsimag[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "DenReYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    densreal[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "DenImYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    densimag[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "CfnReYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    cfctreal[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);
    sprintf(bufname, "CfnImYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    cfctimag[ihist] = new TH1D(bufname, bufname, ibin, vmin, vmax);

    numsreal[ihist]->Sumw2();
    numsimag[ihist]->Sumw2();
    densreal[ihist]->Sumw2();
    densimag[ihist]->Sumw2();
    cfctreal[ihist]->Sumw2();
    cfctimag[ihist]->Sumw2();
  }

  sprintf(bufname, "BinCountNum%s", name);
  binctn = new TH1D(bufname, bufname, ibin, vmin, vmax);

  sprintf(bufname, "BinCountDen%s", name);
  binctd = new TH1D(bufname, bufname, ibin, vmin, vmax);

  ylmbuffer = (complex<double> *) malloc(sizeof(complex<double>) * maxjm);
  
  // Covariance matrices
  covmnum = (double *) malloc(sizeof(double) * maxjm * maxjm * 4 * ibin);
  covmden = (double *) malloc(sizeof(double) * maxjm * maxjm * 4 * ibin);
  covmcfc = (double *) malloc(sizeof(double) * maxjm * maxjm * 4 * ibin);

  for (int iter=0; iter<maxjm*maxjm*4*ibin; iter++) {
    covmnum[iter] = 0.0;
    covmden[iter] = 0.0;
    covmcfc[iter] = 0.0;
  }

  covnum = 0;
  covden = 0;
}


StHbtCorrFctnDirectYlm::StHbtCorrFctnDirectYlm()
{
  StHbtCorrFctnDirectYlm("StHbtCorrFctnDirectYlm",2);
}

StHbtCorrFctnDirectYlm::~StHbtCorrFctnDirectYlm()
{
  for (int ihist=0; ihist<maxjm; ihist++) {
    delete numsreal[ihist];
    delete numsimag[ihist];
    delete densreal[ihist];
    delete densimag[ihist];
    delete cfctreal[ihist];
    delete cfctimag[ihist];
  }

  delete binctn;
  delete binctd;

  //  delete numsreal;
  //  delete numsimag;
  //  delete densreal;
  //  delete densimag;
  //  delete cfctreal;
  //  delete cfctimag;

  free( numsreal);
  free( numsimag);
  free( densreal);
  free( densimag);
  free( cfctreal);
  free( cfctimag);

  free(factorials);
  free(els);
  free(ems);
  free(elsi);
  free(emsi);
  free(ylmbuffer);

  free(covmnum);
  free(covmden);
  free(covmcfc);

  if (covnum) delete covnum;
  if (covden) delete covden;
}

double StHbtCorrFctnDirectYlm::ClebschGordan(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm)
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

double StHbtCorrFctnDirectYlm::DeltaJ(double aJot1, double aJot2, double aJot)
{
  if ((aJot1+aJot2-aJot) < 0) {
    //    cout << "J1+J2-J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((aJot1-aJot2+aJot) < 0) {
    //    cout << "J1-J2+J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((-aJot1+aJot2+aJot) < 0) {
    //    cout << "-J1+J2+J3 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  if ((aJot1+aJot2+aJot+1) < 0) {
    //    cout << "J1+J2+J3+1 < 0 !!!" << " " << aJot1 << " " << aJot2 << " " << aJot << endl;
    return 0;
  }
  double res = TMath::Sqrt(1.0 * 
			   factorials[lrint(aJot1+aJot2-aJot)] * 
			   factorials[lrint(aJot1-aJot2+aJot)] * 
			   factorials[lrint(-aJot1+aJot2+aJot)] / 
			   factorials[lrint(aJot1+aJot2+aJot+1)]);
  
  return res;
}

double StHbtCorrFctnDirectYlm::WignerSymbol(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm)
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


void StHbtCorrFctnDirectYlm::GetMtilde(complex<double> *aMat, double *aMTilde)
{
  // Create the Mtilde for a given q bin
  double lzero, mzero;
  double lprim, mprim;
  double lbis, mbis;
 
  int lzeroi, mzeroi;
  int lprimi, mprimi;
  int lbisi, mbisi;

  for (int iz=0; iz< GetMaxJM()*2; iz++)
    for (int ip=0; ip<GetMaxJM()*2; ip++)
      aMTilde[iz*GetMaxJM()*2+ip] = 0.0;

  for (int izero = 0; izero<GetMaxJM(); izero++) {
    GetElEmForIndex(izero, &lzero, &mzero);
    GetElEmForIndex(izero, &lzeroi, &mzeroi);
    for (int ibis = 0; ibis<GetMaxJM(); ibis++) {
      GetElEmForIndex(ibis, &lbis, &mbis);
      GetElEmForIndex(ibis, &lbisi, &mbisi);
      complex<double> val = complex<double>(0.0, 0.0);
      complex<double> mcomp[MAXJM];
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

int  StHbtCorrFctnDirectYlm::GetMaxJM() 
{ return maxjm; }

void StHbtCorrFctnDirectYlm::GetElEmForIndex(int aIndex, double *aEl, double *aEm)
{
  *aEl = els[aIndex];
  *aEm = ems[aIndex];
}

void StHbtCorrFctnDirectYlm::GetElEmForIndex(int aIndex, int *aEl, int *aEm)
{
  *aEl = elsi[aIndex];
  *aEm = emsi[aIndex];
}

int StHbtCorrFctnDirectYlm::GetBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag)
{
  return (qbin*GetMaxJM()*GetMaxJM()*4 +
	  (ilmprim*2 + primimag) * GetMaxJM()*2 +
	  ilmzero*2 + zeroimag);
}

void StHbtCorrFctnDirectYlm::AddRealPair(double qout, double qside, double qlong, double weight)
{
  double kv = sqrt(qout*qout + qside*qside + qlong*qlong);
  int nqbin = binctn->GetXaxis()->FindFixBin(kv) - 1;

  if (fR2factor) weight*=(1.0/(kv*kv));
  
  StHbtYlm::YlmUpToL(elsi[GetMaxJM()-1], qout, qside, qlong, ylmbuffer);
  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    //    ylmbuffer[ilm] = StHbtYlm::Ylm(elsi[ilm], emsi[ilm], qout, qside, qlong);

    numsreal[ilm]->Fill(kv, real(ylmbuffer[ilm])*weight);
    numsimag[ilm]->Fill(kv, -imag(ylmbuffer[ilm])*weight);

    binctn->Fill(kv, 1.0);
  }

  // Fill in the error matrix
  //  int tabshift = nqbin*GetMaxJM()*GetMaxJM()*4;
  if (nqbin < binctn->GetNbinsX())
    for (int ilmzero=0; ilmzero<GetMaxJM(); ilmzero++)
      for (int ilmprim=0; ilmprim<GetMaxJM(); ilmprim++) {
	covmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] += real(ylmbuffer[ilmzero])*real(ylmbuffer[ilmprim])*weight*weight;
	covmnum[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] += real(ylmbuffer[ilmzero])*-imag(ylmbuffer[ilmprim])*weight*weight;
	covmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] += -imag(ylmbuffer[ilmzero])*real(ylmbuffer[ilmprim])*weight*weight;
	covmnum[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] += -imag(ylmbuffer[ilmzero])*-imag(ylmbuffer[ilmprim])*weight*weight;
	
      }
  
}

void StHbtCorrFctnDirectYlm::AddMixedPair(double qout, double qside, double qlong, double weight)
{
  double kv = sqrt(qout*qout + qside*qside + qlong*qlong);
  
  if (fR2factor) weight*=(1.0/(kv*kv));

  StHbtYlm::YlmUpToL(elsi[GetMaxJM()-1], qout, qside, qlong, ylmbuffer);
  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    //    ylmbuffer[ilm] = StHbtYlm::Ylm(elsi[ilm], emsi[ilm], qout, qside, qlong);

    densreal[ilm]->Fill(kv, real(ylmbuffer[ilm])*weight);
    densimag[ilm]->Fill(kv, -imag(ylmbuffer[ilm])*weight);

    binctd->Fill(kv, 1.0);
  }

  // Fill in the error matrix
  int nqbin = binctn->GetXaxis()->FindFixBin(kv) - 1;
  //  int tabshift = nqbin*GetMaxJM()*GetMaxJM()*4;
  for (int ilmzero=0; ilmzero<GetMaxJM(); ilmzero++)
    for (int ilmprim=0; ilmprim<GetMaxJM(); ilmprim++) {
      covmden[GetBin(nqbin, ilmzero, 0, ilmprim, 0)] += real(ylmbuffer[ilmzero])*real(ylmbuffer[ilmprim]);
      covmden[GetBin(nqbin, ilmzero, 0, ilmprim, 1)] += real(ylmbuffer[ilmzero])*-imag(ylmbuffer[ilmprim]);
      covmden[GetBin(nqbin, ilmzero, 1, ilmprim, 0)] += -imag(ylmbuffer[ilmzero])*real(ylmbuffer[ilmprim]);
      covmden[GetBin(nqbin, ilmzero, 1, ilmprim, 1)] += -imag(ylmbuffer[ilmzero])*-imag(ylmbuffer[ilmprim]);
      
    }
}

void StHbtCorrFctnDirectYlm::AddRealPair(double *qvec, double weight) {
  AddRealPair(qvec[0], qvec[1], qvec[2], weight);
}

void StHbtCorrFctnDirectYlm::AddMixedPair(double *qvec, double weight) {
  AddMixedPair(qvec[0], qvec[1], qvec[2], weight);
}

void StHbtCorrFctnDirectYlm::Finish()
{
  complex<double> tMq0[MAXJM];
  complex<double> tTq0[MAXJM];
  double tMTilde[MAXJM*MAXJM*4];
  //  complex<double> tCq0[maxjm];

  // Create covariance storage elements
  PackCovariances();

  for (int ibin=1; ibin<=numsreal[0]->GetNbinsX(); ibin++) {
    for (int ilm=0; ilm<maxjm; ilm++) {
      tMq0[ilm] = complex<double>(densreal[ilm]->GetBinContent(ibin)/densreal[0]->GetEntries(),
				  densimag[ilm]->GetBinContent(ibin)/densreal[0]->GetEntries());
      tTq0[ilm] = complex<double>(numsreal[ilm]->GetBinContent(ibin)/numsreal[0]->GetEntries(),
				  numsimag[ilm]->GetBinContent(ibin)/numsreal[0]->GetEntries());
    }

    // Calculate the proper error matrix for T
    // from the temporary covariance matrices
    //    int tabshift = (ibin-1)*GetMaxJM()*GetMaxJM()*4;
    for (int ilmzero=0; ilmzero<GetMaxJM(); ilmzero++)
      for (int ilmprim=0; ilmprim<GetMaxJM(); ilmprim++) {
	covmnum[GetBin(ibin-1, ilmzero, 0, ilmprim, 0)] /= numsreal[0]->GetEntries();
	covmnum[GetBin(ibin-1, ilmzero, 0, ilmprim, 1)] /= numsreal[0]->GetEntries();
	covmnum[GetBin(ibin-1, ilmzero, 1, ilmprim, 0)] /= numsreal[0]->GetEntries();
	covmnum[GetBin(ibin-1, ilmzero, 1, ilmprim, 1)] /= numsreal[0]->GetEntries();

	covmnum[GetBin(ibin-1, ilmzero, 0, ilmprim, 0)] -= real(tTq0[ilmzero])*real(tTq0[ilmprim]);
	covmnum[GetBin(ibin-1, ilmzero, 0, ilmprim, 1)] -= real(tTq0[ilmzero])*imag(tTq0[ilmprim]);
	covmnum[GetBin(ibin-1, ilmzero, 1, ilmprim, 0)] -= imag(tTq0[ilmzero])*real(tTq0[ilmprim]);
	covmnum[GetBin(ibin-1, ilmzero, 1, ilmprim, 1)] -= imag(tTq0[ilmzero])*imag(tTq0[ilmprim]);

	covmnum[GetBin(ibin-1, ilmzero, 0, ilmprim, 0)] /= (numsreal[0]->GetEntries() - 1);
	covmnum[GetBin(ibin-1, ilmzero, 0, ilmprim, 1)] /= (numsreal[0]->GetEntries() - 1);
	covmnum[GetBin(ibin-1, ilmzero, 1, ilmprim, 0)] /= (numsreal[0]->GetEntries() - 1);
	covmnum[GetBin(ibin-1, ilmzero, 1, ilmprim, 1)] /= (numsreal[0]->GetEntries() - 1);
      }

    GetMtilde(tMq0, tMTilde);
  
    // Perform the solution for the correlation function itself and the errors

    if (0) { 
      // This is the "wrong" way 
      // do not use it
      // Fill in the M matrix
      double mM[MAXJM*MAXJM*4];
      for (int iter=0; iter<maxjm*maxjm*4; iter++)
	mM[iter] = tMTilde[iter];
      
      gsl_matrix_view matM = gsl_matrix_view_array(mM, GetMaxJM()*2, GetMaxJM()*2);
      
      // Prepare halper matrices
      double mU[MAXJM*MAXJM*4];
      for (int iter=0; iter<maxjm*maxjm*4; iter++)
	mU[iter] = tMTilde[iter];

      gsl_matrix_view matU = gsl_matrix_view_array(mU, GetMaxJM()*2, GetMaxJM()*2);

      double mV[MAXJM*MAXJM*4];
      gsl_matrix_view matV = gsl_matrix_view_array(mV, GetMaxJM()*2, GetMaxJM()*2);

      double vS[MAXJM*2];
      gsl_vector_view vecS = gsl_vector_view_array(vS, GetMaxJM()*2);

      double vW[MAXJM*2];
      gsl_vector_view vecW = gsl_vector_view_array(vW, GetMaxJM()*2);

      // Decomposing the M matrix
      gsl_linalg_SV_decomp(&matU.matrix, &matV.matrix, &vecS.vector, &vecW.vector);
  
      double vB[MAXJM*2];
      for (int iter=0; iter<GetMaxJM(); iter++) {
	vB[iter*2] = real(tTq0[iter]);
	vB[iter*2+1] = imag(tTq0[iter]);
      }

      // Prepare inputs for solving the problem
      gsl_vector_view vecB = gsl_vector_view_array(vB, GetMaxJM()*2);

      double vX[MAXJM];
      gsl_vector_view vecX = gsl_vector_view_array(vX, GetMaxJM()*2);

      // Solving the problem
      gsl_linalg_SV_solve(&matU.matrix, &matV.matrix, &vecS.vector, &vecB.vector, &vecX.vector);

      for (int ilm=0; ilm<maxjm; ilm++) {
	cfctreal[ilm]->SetBinContent(ibin, vX[ilm*2]);
	cfctimag[ilm]->SetBinContent(ibin, vX[ilm*2+1]);
      }
  
      // Perform the solutions for the error matrix as well
      
      for (int iter=0; iter<maxjm*maxjm*4; iter++)
	// We must always add errors
	mU[iter] = (tMTilde[iter]);
    
      // Decomposing the M matrix
      gsl_linalg_SV_decomp(&matU.matrix, &matV.matrix, &vecS.vector, &vecW.vector);

      for (int iterlmz=0; iterlmz<GetMaxJM()*2; iterlmz++) {
	int ilmz = iterlmz/2;
	int imgz = iterlmz%2;
	for (int ilm=0; ilm<GetMaxJM(); ilm++) {
	  vB[ilm*2]   = sqrt(fabs(covmnum[GetBin(ibin-1, ilmz, imgz, ilm, 0)]));
	  vB[ilm*2+1] = sqrt(fabs(covmnum[GetBin(ibin-1, ilmz, imgz, ilm, 1)]));
	}
	gsl_linalg_SV_solve(&matU.matrix, &matV.matrix, &vecS.vector, &vecB.vector, &vecX.vector);
      
	for (int ilm=0; ilm<GetMaxJM(); ilm++) {
	  if (!(isnan(vX[ilm*2])))
	    covmcfc[GetBin(ibin-1, ilmz, imgz, ilm, 0)] = vX[ilm*2];
	  else
	    covmcfc[GetBin(ibin-1, ilmz, imgz, ilm, 0)] = 1e-12;
	  if (!(isnan(vX[ilm*2+1])))
	    covmcfc[GetBin(ibin-1, ilmz, imgz, ilm, 1)] = vX[ilm*2+1];
	  else
	    covmcfc[GetBin(ibin-1, ilmz, imgz, ilm, 1)] = 1e-12;
	}
      
      }
  
      for (int ilm=0; ilm<maxjm; ilm++) {
	cfctreal[ilm]->SetBinError(ibin, (fabs(covmcfc[GetBin(ibin-1, ilm, 0, ilm, 0)])));
	cfctimag[ilm]->SetBinError(ibin, (fabs(covmcfc[GetBin(ibin-1, ilm, 1, ilm, 1)])));
      }
    }
    
    cout << "C calculation for bin " << (ibin-1) << endl;
    
    if (1) {
      // Here is the "right" way from recent Dave's note
      
      // First we must calculate the Delta^2 C matrix of errors
      
      // First we must have the Delta^2 T matrix
      double mDeltaT[MAXJM*MAXJM*4];
      for (int ilmzero=0; ilmzero<GetMaxJM()*2; ilmzero++)
	for (int ilmprim=0; ilmprim<GetMaxJM()*2; ilmprim++) 
	  mDeltaT[(ilmzero*maxjm*2) + ilmprim] = covmnum[GetBin(ibin-1,ilmzero/2, ilmzero%2, ilmprim/2, ilmprim%2)];
      
      cout << "Delta T matrix " << endl;
      for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++) {
	for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++) {
	  cout.precision(3);
	  cout.width(10);
	  cout << mDeltaT[ilmz*GetMaxJM()*2 + ilmp];
	}
	cout << endl;
      }

      cout << "Delta T matrix packed " << endl;
      double mDeltaTPacked[MAXJM*MAXJM*4];
      int msize = PackYlmMatrix(mDeltaT, mDeltaTPacked);
      cout << "Packed size is " << msize << endl;
      for (int ilmz=0; ilmz<msize; ilmz++) {
	for (int ilmp=0; ilmp<msize; ilmp++) {
	  cout.precision(3);
	  cout.width(10);
	  cout << mDeltaTPacked[ilmz*msize + ilmp];
	}
	cout << endl;
      }

      gsl_matrix_view matDeltaT = gsl_matrix_view_array(mDeltaTPacked, msize, msize);

      // Prepare halper matrices
      double mU[MAXJM*MAXJM*4];
      for (int iter=0; iter<msize*msize; iter++)
	mU[iter] = mDeltaTPacked[iter];

      gsl_matrix_view matU = gsl_matrix_view_array(mU, msize, msize);

      double mM[MAXJM*MAXJM*4];
      double mMPacked[MAXJM*MAXJM*4];
      for (int iter=0; iter<maxjm*maxjm*4; iter++)
 	mM[iter] = tMTilde[iter];
      PackYlmMatrix(mM, mMPacked);
      
      gsl_matrix_view matM = gsl_matrix_view_array(mMPacked, msize, msize);

     cout << "Mtilde matrix " << endl;
       for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++) {
 	for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++) {
 	  cout.precision(3);
 	  cout.width(10);
 	  cout << mM[ilmz*GetMaxJM()*2 + ilmp];
 	}
 	cout << endl;
       }

      cout << "Mtilde matrix packed " << endl;
      for (int ilmz=0; ilmz<msize; ilmz++) {
	for (int ilmp=0; ilmp<msize; ilmp++) {
	  cout.precision(3);
	  cout.width(10);
	  cout << mMPacked[ilmz*msize + ilmp];
	}
	cout << endl;
      }

       // DeltaT * Mtilde = U 

      gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &matU.matrix, &matM.matrix);

       // Mtilde * U = V

      double mX[MAXJM*MAXJM*4];
      for (int iter=0; iter<maxjm*maxjm*4; iter++)
 	mX[iter] = tMTilde[iter];
      
      double mXPacked[MAXJM*MAXJM*4];
      PackYlmMatrix(mX, mXPacked);
      gsl_matrix_view matX = gsl_matrix_view_array(mXPacked, msize, msize);

       double mV[MAXJM*MAXJM*4];
       for (int iter=0; iter<maxjm*maxjm*4; iter++)
	 mV[iter] = 0;
      
       gsl_matrix_view matV = gsl_matrix_view_array(mV, msize, msize);

       gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &matX.matrix, &matM.matrix, 0.0, &matV.matrix);

       // DeltaT * T = X
      
       double vB[MAXJM*2];
       for (int iter=0; iter<GetMaxJM(); iter++) {
	 vB[iter*2] = real(tTq0[iter]);
	 vB[iter*2+1] = imag(tTq0[iter]);
       }
       
       double vBPacked[MAXJM*2];
       PackYlmVector(vB, vBPacked);


       // Prepare inputs for solving the problem
       gsl_vector_view vecB = gsl_vector_view_array(vBPacked, msize);
       
       cout << "B vector " << endl;
       for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++) {
	 cout.precision(3);
	 cout.width(10);
	 cout << vB[ilmz];
       }
       cout << endl;

      cout << "B vector packed " << endl;
      for (int ilmp=0; ilmp<msize; ilmp++) {
	cout.precision(3);
	cout.width(10);
	cout << vBPacked[ilmp];
      }
      cout << endl;

      gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &matU.matrix, &vecB.vector);
      
      cout << "B vector packed solution" << endl;
      for (int ilmp=0; ilmp<msize; ilmp++) {
	cout.precision(3);
	cout.width(10);
	cout << vBPacked[ilmp];
      }
      cout << endl;

      // MtildeT * X = Y

       double vY[MAXJM*2];
       for (int iter=0; iter<GetMaxJM()*2; iter++) {
	 vY[iter] = 0.0;
       }
      
       // Prepare inputs for solving the problem
       gsl_vector_view vecY = gsl_vector_view_array(vY, msize);

       gsl_blas_dgemv (CblasTrans, 1.0, &matX.matrix, &vecB.vector, 0.0, &vecY.vector);

      cout << "Y vector packed solution" << endl;
      for (int ilmp=0; ilmp<msize; ilmp++) {
	cout.precision(3);
	cout.width(10);
	cout << vY[ilmp];
      }
      cout << endl;

      double vC[MAXJM*2];
      for (int iter=0; iter<GetMaxJM()*2; iter++) {
 	vC[iter] = vY[iter];
      }
      
       // Prepare inputs for solving the problem
      gsl_vector_view vecC = gsl_vector_view_array(vC, msize);
      gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, &matV.matrix, &vecC.vector);

      cout << "C vector packed solution" << endl;
      for (int ilmp=0; ilmp<msize; ilmp++) {
	cout.precision(3);
	cout.width(10);
	cout << vC[ilmp];
      }
      cout << endl;
      
      int mpack=0;
      int el, em;
      for (int ilm=0; ilm<maxjm; ilm++) {
 	cfctreal[ilm]->SetBinContent(ibin, vC[mpack++]);
	GetElEmForIndex(ilm, &el, &em);
	if(em==0)
	  cfctimag[ilm]->SetBinContent(ibin, 0);
	else
	  cfctimag[ilm]->SetBinContent(ibin, vC[mpack++]);
      }

      // invert the V matrix to get C errors
      double mS[MAXJM*MAXJM*4];

      for (int iterz=0; iterz<msize; iterz++)
	for (int iterp=0; iterp<msize; iterp++)
	  if (iterp == iterz)
	    mS[iterz*msize + iterp] = 1.0;
	  else
	    mS[iterz*msize + iterp] = 0.0;
      
      gsl_matrix_view matS = gsl_matrix_view_array(mS, msize, msize);
      
      // Invert V 

      gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &matV.matrix, &matS.matrix);

      mpack = 0;
      for (int ilm=0; ilm<maxjm; ilm++) {
 	cfctreal[ilm]->SetBinError(ibin, sqrt(mS[mpack*msize+mpack]));
	mpack++;
	GetElEmForIndex(ilm, &el, &em);
	if(em==0)
	  cfctimag[ilm]->SetBinError(ibin, 0);
	else {
	  cfctimag[ilm]->SetBinError(ibin, sqrt(mS[mpack*msize+mpack]));
	  mpack++;
	}
      }

      

//       *** END of new method ***

//       double vS[MAXJM*2];
//       gsl_vector_view vecS = gsl_vector_view_array(vS, GetMaxJM()*2);


//       double mV[MAXJM*MAXJM*4];
//       gsl_matrix_view matV = gsl_matrix_view_array(mV, GetMaxJM()*2, GetMaxJM()*2);

//       double vS[MAXJM*2];
//       gsl_vector_view vecS = gsl_vector_view_array(vS, GetMaxJM()*2);

//       double vW[MAXJM*2];
//       gsl_vector_view vecW = gsl_vector_view_array(vW, GetMaxJM()*2);
    }
    
  }

}

void StHbtCorrFctnDirectYlm::Write()
{
  for (int ilm=0; ilm<maxjm; ilm++) {
    numsreal[ilm]->Write();
    densreal[ilm]->Write();
    cfctreal[ilm]->Write();
    numsimag[ilm]->Write();
    densimag[ilm]->Write();
    cfctimag[ilm]->Write();
  }
  if (covnum) covnum->Write();
  if (covden) covden->Write();
}

void StHbtCorrFctnDirectYlm::ReadFromFile(TFile *infile, const char *name, int maxl)
{
  if (maxl != fMaxL) {
    cout << "Cannot read function for L " << maxl << " into a container with L "<< fMaxL << endl;
    return;
  }
  cout << "Reading in numerators and denominators" << endl;

  char bufname[200];
  for (int ihist=0; ihist<maxjm; ihist++) {
    sprintf(bufname, "NumReYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    if (numsreal[ihist]) delete numsreal[ihist];
    numsreal[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));

    sprintf(bufname, "NumImYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    if (numsimag[ihist]) delete numsimag[ihist];
    numsimag[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));

    sprintf(bufname, "DenReYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    if (densreal[ihist]) delete densreal[ihist];
    densreal[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));

    sprintf(bufname, "DenImYlm%i%i%s", elsi[ihist], emsi[ihist]<0 ? elsi[ihist]-emsi[ihist] : emsi[ihist], name);
    if (densimag[ihist]) delete densimag[ihist];
    densimag[ihist] = new TH1D(*((TH1D *) infile->Get(bufname)));
  }

  if (covnum) delete covnum;
  sprintf(bufname, "CovNum%s", name);
  covnum = new TH3D (*((TH3D *) infile->Get(bufname)));

  if (covden) delete covden;
  sprintf(bufname, "CovDen%s", name);
  covden = new TH3D (*((TH3D *) infile->Get(bufname)));

  if ((covnum) && (covden)) {
    cout << "Unpacking covariance matrices from file " << endl;
    UnpackCovariances();
  }
  else {

    cout << "Creating fake covariance matrices" << endl;
  
    for (int ibin=1; ibin<=numsreal[0]->GetNbinsX(); ibin++) {
      double nent = numsreal[0]->GetEntries();
      double nentd = densreal[0]->GetEntries();
      for (int ilmx=0; ilmx<GetMaxJM(); ilmx++) {
	for (int ilmy=0; ilmy<GetMaxJM(); ilmy++) {
	  double t1t2rr = numsreal[ilmx]->GetBinContent(ibin)*numsreal[ilmy]->GetBinContent(ibin)/nent/nent;
	  double t1t2ri = numsreal[ilmx]->GetBinContent(ibin)*numsimag[ilmy]->GetBinContent(ibin)/nent/nent;
	  double t1t2ir = numsimag[ilmx]->GetBinContent(ibin)*numsreal[ilmy]->GetBinContent(ibin)/nent/nent;
	  double t1t2ii = numsimag[ilmx]->GetBinContent(ibin)*numsimag[ilmy]->GetBinContent(ibin)/nent/nent;
	  if (ilmx == ilmy) {
	    covmnum[GetBin(ibin-1, ilmx, 0, ilmy, 0)] = nent*(TMath::Power(numsreal[ilmx]->GetBinError(ibin)/nent,2)*(nent-1) + t1t2rr);
	    covmnum[GetBin(ibin-1, ilmx, 0, ilmy, 1)] = nent*t1t2ri;
	    covmnum[GetBin(ibin-1, ilmx, 1, ilmy, 0)] = nent*t1t2ir;
	    covmnum[GetBin(ibin-1, ilmx, 1, ilmy, 1)] = nent*(TMath::Power(numsimag[ilmx]->GetBinError(ibin)/nent,2)*(nent-1) + t1t2rr);
	  }
	  else {
	    covmnum[GetBin(ibin-1, ilmx, 0, ilmy, 0)] = nent*t1t2rr;
	    covmnum[GetBin(ibin-1, ilmx, 0, ilmy, 1)] = nent*t1t2ri;
	    covmnum[GetBin(ibin-1, ilmx, 1, ilmy, 0)] = nent*t1t2ir;
	    covmnum[GetBin(ibin-1, ilmx, 1, ilmy, 1)] = nent*t1t2ii;
	  }
	  t1t2rr = densreal[ilmx]->GetBinContent(ibin)*densreal[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  t1t2ri = densreal[ilmx]->GetBinContent(ibin)*densimag[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  t1t2ir = densimag[ilmx]->GetBinContent(ibin)*densreal[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  t1t2ii = densimag[ilmx]->GetBinContent(ibin)*densimag[ilmy]->GetBinContent(ibin)/nentd/nentd;
	  
	  covmden[GetBin(ibin-1, ilmx, 0, ilmy, 0)] = nentd*t1t2rr;
	  covmden[GetBin(ibin-1, ilmx, 0, ilmy, 1)] = nentd*t1t2ri;
	  covmden[GetBin(ibin-1, ilmx, 1, ilmy, 0)] = nentd*t1t2ir;
	  covmden[GetBin(ibin-1, ilmx, 1, ilmy, 1)] = nentd*t1t2ii;
	}
      }
    }
  }

  // Recalculating the correlation functions
  Finish();
}

int StHbtCorrFctnDirectYlm::PackYlmVector(double *invec, double *outvec)
{
  int ioutcount = 0;
  int em, el;
  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    GetElEmForIndex(ilm, &el, &em);
    outvec[ioutcount++] = invec[ilm*2];
    if (em == 0)
      continue;
    outvec[ioutcount++] = invec[ilm*2 + 1];
  }
  
  return ioutcount;
}

int StHbtCorrFctnDirectYlm::PackYlmMatrix(double *inmat, double *outmat)
{
  int ioutcountz = 0;
  int ioutcountp = 0;
  int emz, elz;
  int emp, elp;
  int finalsize = 0;

  for (int ilm=0; ilm<GetMaxJM(); ilm++) {
    GetElEmForIndex(ilm, &elz, &emz);
    finalsize++;
    if (emz == 0) continue;
    finalsize++;
  }

  for (int ilmz=0; ilmz<GetMaxJM(); ilmz++) {
    GetElEmForIndex(ilmz, &elz, &emz);
    ioutcountp = 0;
    for (int ilmp=0; ilmp<GetMaxJM(); ilmp++) {
      GetElEmForIndex(ilmp, &elp, &emp);
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 0, ilmp, 0)];
      ioutcountp++;
      if (emp == 0) continue;
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 0, ilmp, 1)];
      ioutcountp++;
    }
    ioutcountz++;

    if (emz == 0) continue;
    ioutcountp = 0;
    for (int ilmp=0; ilmp<GetMaxJM(); ilmp++) {
      GetElEmForIndex(ilmp, &elp, &emp);
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 1, ilmp, 0)];
      ioutcountp++;
      if (emp == 0) continue;
      outmat[ioutcountz*finalsize + ioutcountp] = inmat[GetBin(0, ilmz, 1, ilmp, 1)];
      ioutcountp++;
    }
    ioutcountz++;    
  }	
  
  return ioutcountz;  
}

void StHbtCorrFctnDirectYlm::PackCovariances()
{
  char bufname[200];
  sprintf(bufname, "CovNum%s", numsreal[0]->GetName()+10);

  if (covnum) delete covnum;
  covnum = new TH3D(bufname,bufname, 
		    numsreal[0]->GetNbinsX(), numsreal[0]->GetXaxis()->GetXmin(), numsreal[0]->GetXaxis()->GetXmax(),
		    GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5,
		    GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5);
  
  for (int ibin=1; ibin<=covnum->GetNbinsX(); ibin++)
    for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
      for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	covnum->SetBinContent(ibin, ilmz+1, ilmp+1, covmnum[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)]);

  if (covden) delete covden;
  sprintf(bufname, "CovDen%s", numsreal[0]->GetName()+10);
  covden  = new TH3D(bufname,bufname, 
		     densreal[0]->GetNbinsX(), densreal[0]->GetXaxis()->GetXmin(), densreal[0]->GetXaxis()->GetXmax(),
		     GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5,
		     GetMaxJM()*2, -0.5, GetMaxJM()*2 - 0.5);
		     
  for (int ibin=1; ibin<=covden->GetNbinsX(); ibin++)
    for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
      for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	covden->SetBinContent(ibin, ilmz+1, ilmp+1, covmden[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)]);

}

void StHbtCorrFctnDirectYlm::UnpackCovariances()
{
  if (covnum) {
    for (int ibin=1; ibin<=covnum->GetNbinsX(); ibin++)
      for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
	for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	  covmnum[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)] = covnum->GetBinContent(ibin, ilmz+1, ilmp+1);
    
  }
  if (covden) {
    for (int ibin=1; ibin<=covden->GetNbinsX(); ibin++)
      for (int ilmz=0; ilmz<GetMaxJM()*2; ilmz++)
	for (int ilmp=0; ilmp<GetMaxJM()*2; ilmp++)
	  covmden[GetBin(ibin-1, ilmz/2, ilmz%2, ilmp/2, ilmp%2)] = covden->GetBinContent(ibin, ilmz+1, ilmp+1);
  }
}

int StHbtCorrFctnDirectYlm::GetIndexForLM(int el, int em)
{
  for (int iter=0; iter<maxjm; iter++)
    if ((el == elsi[iter]) && (em == emsi[iter]))
      return iter;
  return -1;
}

TH1D *StHbtCorrFctnDirectYlm::GetNumRealHist(int el, int em)
{
  if (GetIndexForLM(el, em)>=0)
    return numsreal[GetIndexForLM(el, em)];
  else 
    return 0;
}
TH1D *StHbtCorrFctnDirectYlm::GetNumImagHist(int el, int em)
{
  if (GetIndexForLM(el, em)>=0)
    return numsimag[GetIndexForLM(el, em)];
  else 
    return 0;
}

TH1D *StHbtCorrFctnDirectYlm::GetDenRealHist(int el, int em)
{
  if (GetIndexForLM(el, em)>=0)
    return densreal[GetIndexForLM(el, em)];
  else 
    return 0;
}
TH1D *StHbtCorrFctnDirectYlm::GetDenImagHist(int el, int em)
{
  if (GetIndexForLM(el, em)>=0)
    return densimag[GetIndexForLM(el, em)];
  else 
    return 0;
}

StHbtString StHbtCorrFctnDirectYlm::Report()
{
  return "StHbtCorrFctnDirectYlm::Finish";
}

void StHbtCorrFctnDirectYlm::AddRealPair(const StHbtPair* aPair)
{
  AddRealPair(aPair->dKOut(), aPair->dKSide(), aPair->dKLong(), 1.0);
}
void StHbtCorrFctnDirectYlm::AddMixedPair(const StHbtPair* aPair)
{
  AddMixedPair(aPair->dKOut(), aPair->dKSide(), aPair->dKLong(), 1.0);
}

void StHbtCorrFctnDirectYlm::SetR2Factor(bool aFactor)
{
  fR2factor = aFactor;
}
