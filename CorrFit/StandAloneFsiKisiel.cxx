#include "StandAloneFsiKisiel.h"

#include "fstream"

#define COULOMBSTEPS 120

extern ReadPar *sRPInstance;
extern STR      sRPFileName;

const double mPi = 0.13957;
const double mK  = 0.493677;
const double mK0 = 0.497648;
//const double mP  = 0.938272;
const double mP  = 0.93827231;
const double mLa = 1.115683;
const double mFmToGeV = 0.197326968;
//const double mFmToGeV = 0.1973;
const double mAlfa = 137.036;
long double mEta = 0.0;

dcomplex mult(dcomplex arga, dcomplex argb)
{
  dcomplex res;
  
  res.re = arga.re * argb.re - arga.im * argb.im;
  res.im = arga.re * argb.im + argb.re * arga.im;

  return res;
}

dcomplex conj(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re;
  res.im = -arg.im;

  return res;
}

dcomplex mult(dcomplex arga, long double argb)
{
  dcomplex res;
  
  res.re = arga.re * argb;
  res.im = arga.im * argb;

  return res;
}

long double modl2(dcomplex arg)
{
  return arg.re*arg.re + arg.im*arg.im;
}

long double modl(dcomplex arg)
{
  return hypot(arg.re, arg.im);
}

dcomplex invr(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re/modl2(arg);
  res.im = - arg.im/modl2(arg);

  return res;
}

StandAloneFsiKisiel::StandAloneFsiKisiel()
{
  mLL = 0;
  mPrint = 0;
}

double StandAloneFsiKisiel::getWeight(Pair& aPair)
{
  PairKinematics(aPair);
  if (mIqs && !mIch && !mIsi) return GetQuantum(aPair);
  if (mIch && !mIqs && !mIsi) return GetCoulomb(aPair);
  if (mIch &&  mIqs && !mIsi) return GetQuantumCoulomb(aPair);
  if (mIch &&  mIqs &&  mIsi) return GetQuantumCoulombStrong(aPair);
  else return 1.0;
}

void StandAloneFsiKisiel::InitializeGamow()
{
  long double twopioverac = 2.0*TMath::Pi()/mBohrRadius;
  long double tpaoverk;
  
  for (int iter=0; iter<4000; iter++) {
    tpaoverk = twopioverac/(iter*0.0002 + 0.0001);
    mGamovMesh[iter] = tpaoverk * 1.0 / (exp(tpaoverk) - 1);
  }
  mEta = twopioverac;
}

long double StandAloneFsiKisiel::Gamow(double arg)
{
  return ((mEta/arg) * 1.0/(exp(mEta/ arg) - 1.0));
}

// Calculates the confluent hypergeometric function F
// from single orientation of cos(theta*)
// For non-symmetrized wave-function (non-identical particles)
void StandAloneFsiKisiel::GetFFsingle(dcomplex *ffp, int sign, Pair &aPair)
{
  long double comprep[COULOMBSTEPS];
  long double compimp[COULOMBSTEPS];
  long double eta, ksip;
  dcomplex alfa, zetp;

  int nsteps;
  
  long double kstar = aPair.GetKStar();
  long double tKstRst = aPair.GetKStarOut()*mROS + aPair.GetKStarSide()*mRSS + aPair.GetKStarLong()*mRLS;
  long double coskr = sign * tKstRst/(aPair.GetKStar() * mRSt);

  if (mPrint) {
    cout << "kr veckvecr ";
    cout.precision(10);
    cout << (kstar*mRSt);
    cout << " ";
    cout.precision(10);
    cout << tKstRst << " " << (kstar*mRSt + tKstRst) << endl;
    cout << "1/Ka " << (1.0/(kstar * mBohrRadius)) << endl;
  }

  if (kstar*mRSt*(1+coskr) > 10.0)
    nsteps = 110;
  else if (kstar*mRSt*(1+coskr) > 5.0)
    nsteps = 45;
  else
    nsteps = 25;

  eta = 1.0/(kstar * mBohrRadius);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp;
  long double tcomp;
  dcomplex sump;
  dcomplex fcmult;

  long double rad = mRSt;

  ksip = kstar*rad*(1+coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);

    sump = mult(sump, 1.0/(tcomp*tcomp));
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    tcomp *= (istep+1);
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1];

  if (mPrint) {
    PRINT_MESSAGE(ffp->re << " " << ffp->im);
  }
}

// Calculates the confluent hypergeometric function
// For two orientations of cos(theta*) 
// For symmetrized wave-function (identical particles)
void StandAloneFsiKisiel::GetFFdouble(dcomplex *ffp, dcomplex *ffm, Pair &aPair)
{
  long double comprep[COULOMBSTEPS];
  long double compimp[COULOMBSTEPS];
  long double comprem[COULOMBSTEPS];
  long double compimm[COULOMBSTEPS];
  long double eta, ksip, ksim;
  dcomplex alfa, zetp, zetm;
  
  long double kstar = aPair.GetKStar();
  long double tKstRst = aPair.GetKStarOut()*mROS + aPair.GetKStarSide()*mRSS + aPair.GetKStarLong()*mRLS;
  long double coskr = tKstRst/(aPair.GetKStar() * mRSt);

  if (mPrint) {
    cout << "kr veckvecr +";
    cout.precision(10);
    cout << (kstar*mRSt);
    cout << " ";
    cout.precision(10);
    cout << tKstRst << " " << (kstar*mRSt + tKstRst) << endl;
    //    PRINT_MESSAGE("kr veckvecr " << (kstar*mRSt) << " " << tKstRst);
  }
  int nsteps;

  if ((kstar*mRSt*(1+coskr) < 5.0) &&
      (kstar*mRSt*(1-coskr) < 5.0))
    nsteps = 25;
  else if ((kstar*mRSt*(1+coskr) < 10.0) &&
	   (kstar*mRSt*(1-coskr) < 10.0))
    nsteps = 45;
  else if ((kstar*mRSt*(1+coskr) < 15.0) &&
	   (kstar*mRSt*(1-coskr) < 15.0))
    nsteps = 110;
  else
    nsteps = 150;
  eta = 1.0/(kstar * mBohrRadius);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp, scompm;
  long double tcomp;
  dcomplex sump, summ;
  dcomplex fcmult;

  long double rad = mRSt;

  ksip = kstar*rad*(1+coskr);
  ksim = kstar*rad*(1-coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  zetm.re = 0.0;
  zetm.im = ksim;

  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  scompm.re = 1.0; 
  scompm.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);
    summ = mult(fcomp, scompm);

    sump = mult(sump, 1.0/(tcomp*tcomp));
    summ = mult(summ, 1.0/(tcomp*tcomp));
	
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
      
      comprem[istep] = summ.re;
      compimm[istep] = summ.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
      
      comprem[istep] = comprem[istep-1] + summ.re;
      compimm[istep] = compimm[istep-1] + summ.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    scompm = mult(scompm, zetm);
    tcomp *= (istep+1);
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1];

  ffm->re = comprem[nsteps-1];
  ffm->im = compimm[nsteps-1];

  if (mPrint) {
    PRINT_MESSAGE(ffp->re << " " << ffp->im);
    PRINT_MESSAGE(ffm->re << " " << ffm->im);
  }
}

// Calculates the h function for the strong interaction
long double StandAloneFsiKisiel::geth(long double eta)
{
  long double etasum = log(1.0/eta) - euler;
  long double series = 0.0;
  long double x2inv = (eta*eta);
  long double element;
  long double save;
  for (int iter=1; iter<1000000; iter++) {
    element = ((1.0*iter)*(1.0*iter) + x2inv)*(1.0*iter);
    //    cout << "Element " << iter << " is " << element << endl;
    element = 1.0/element;
    if (iter == 1) save = 1.0e-10 * element;
    //    cout << "Element " << iter << " is " << element << endl;

    series += element;
    if (element < save) break;
  }
  series *= x2inv;
  etasum += series;

  return etasum;
}

long double StandAloneFsiKisiel::chiim(long double eta)
{
  return Gamow(1.0/(eta*mBohrRadius))/(2.0*eta);
}

void StandAloneFsiKisiel::bfunpfun(long double eta, long double rho, long double &bret, long double &pret)
{
  long double b0 = 1;
  long double b1 = eta*rho;
  long double bsum = b0+b1;
  long double bnpu, bn, bnmu;
  long double p0 = 1.0;
  long double p1 = 0.0;
  long double psum = p0+p1;
  long double pnpu, pn, pnmu;

  bn = b1; bnmu = b0;
  pn = p1; pnmu = p0;
  for (int iter=1; iter<100; iter++) {
    bnpu = 2 * eta*rho *bn - rho*rho*bnmu;
    bnpu /= (1.0*iter+1.0)*(1.0*iter+2.0);
    bsum += bnpu;
    //    cout << "B E " << iter << " " << bnpu << endl;

    pnpu = 2 * eta*rho *pn - rho*rho*pnmu - (2.0*iter+1.0)*2.0*eta*rho*bn;
    pnpu /= (1.0*iter)*(1.0*iter+1.0);
    psum += pnpu;
    //    cout << "P E " << iter << " " << pnpu << endl;

    bnmu = bn;
    bn   = bnpu;

    pnmu = pn;
    pn   = pnpu;
    if ((fabs(bnpu) + fabs(pnpu)) < 1.0e-20) break;
  }

  bret = bsum;
  pret = psum;
}

// Calculate the Gtilde function
dcomplex StandAloneFsiKisiel::GetG(long double eta, long double rho, long double hfun)
{
  dcomplex gtemp;
  long double bres, pres;
  long double bmult;

  bfunpfun(eta, rho, bres, pres);
  bmult = 2.0*eta*rho*bres;

  gtemp.re = pres + bmult*(log(fabs(2.0*eta*rho))+2.0*euler-1.0+hfun);
  gtemp.im = bmult * chiim(eta);

  if (mPrint) {
    cout << "Gre bmult " << gtemp.re << " " << bmult << endl;
    cout << "2er r " << 2.0*eta*rho << " " << rho << endl;
    cout << "p b " << pres << " " << bres << endl;
  }

  return gtemp;
}

void StandAloneFsiKisiel::Getfc(long double kstar, long double eta, long double hfun, dcomplex &fcs, dcomplex &fct)
{
  dcomplex ci;
  dcomplex cia;

  dcomplex fis;
  dcomplex dis;
  dcomplex fcinvs;

  dcomplex fit;
  dcomplex dit;
  dcomplex fcinvt;

  ci.re = hfun;
  ci.im = chiim(eta);
  cia = mult(ci, 2.0/mBohrRadius); 

  fis = invr(f0s);
  dis = mult(d0s, 0.5 * kstar * kstar);
  fcinvs.re = fis.re + dis.re - cia.re;
  fcinvs.im = fis.im + dis.im - cia.im;
  fcs = invr(fcinvs);
  
  fit = invr(f0t);
  dit = mult(d0t, 0.5 * kstar * kstar);
  fcinvt.re = fit.re + dit.re - cia.re;
  fcinvt.im = fit.im + dit.im - cia.im;
  fct = invr(fcinvt);
  
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized) or fermions (antisymmetrized)
// with Coulomb interaction and Strong interaction included
double StandAloneFsiKisiel::GetQuantumCoulombStrong(Pair &aPair)
{
  if (mRSt < 0.0000000001)
    return 1.0;

  if (mRSt < 1.0/0.197327)
    return GetQuantumCoulomb(aPair);

  double tKstRst = aPair.GetKStarOut()*mROS + aPair.GetKStarSide()*mRSS + aPair.GetKStarLong()*mRLS;
  long double kstar = aPair.GetKStar();
  long double rho = mRSt * kstar;
  
  int ccase = 0;
  static int pcount = 0;
  int wavesign = 1;

  if (mPrint == -1)
    mPrint = 1;
  else
    mPrint = 0;

  if (mPrint) {
    PRINT_MESSAGE("KstRst " << tKstRst);
    PRINT_MESSAGE("ROut RSide RLong DTime " << mROS << " " << mRSS << " " << mRLS << " " << mRTS);
  }
  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  //  if (fabs(mRSt) > fabs(mBohrRadius)) return (1.0 + wavesign*cos(2*tKstRst));

  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  long double testp = rho * (1.0 + tKstRst/(rho));
  long double testm = rho * (1.0 - tKstRst/(rho));

  dcomplex ffplus, ffminus;
  if ((testp> 15.0) && (testm> 15.0))
    {
      double asym;
      asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/rho)*mBohrRadius*kstar*kstar))/Gamow(kstar);
      //      cout << "as1 " << asym << endl;
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffminus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffminus.re = 1.0 + (asym -1.0) /2.0;
      ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);

      asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/rho)*mBohrRadius*kstar*kstar))/Gamow(kstar);
      //      cout << "as2 " << asym << endl;
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffplus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffplus.re = 1.0 + (asym -1.0) /2.0;
      ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      
    }
  
  // Check for the classical limit in both functions separately
  else if (((testp< 15.0) && (testm< 15.0))) // ||
    {
      // Calculate the F function
      GetFFdouble(&ffplus, &ffminus, aPair);
      ccase = 1;
    }
  else if (testp< 15.0)
    {
      double asym;
      GetFFsingle(&ffplus,1,aPair);
      GetFFsingle(&ffminus, -1, aPair);
      if ((fabs(ffminus.re) > 10.0) || fabs(ffminus.im) > 10.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/(rho)*mBohrRadius*kstar*kstar)))/Gamow(kstar);
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffminus.re = 1.0 + (asym -1.0) *2.0;
	else
	  ffminus.re = 1.0 + (asym -1.0) /2.0;
	ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);
      }
      ccase = 2;
    }
  else 
    {
      double asym;
      GetFFsingle(&ffminus, -1, aPair);
      GetFFsingle(&ffplus,   1, aPair);
      if ((fabs(ffplus.re) > 10.0) || fabs(ffplus.im) > 10.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/(rho)*mBohrRadius*kstar*kstar)))/Gamow(kstar);
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffplus.re = 1.0 + (asym -1.0) *2.0;
	else
	  ffplus.re = 1.0 + (asym -1.0) /2.0;
	ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      }
      ccase = 3;
    }

  long double eta  = 1.0/(kstar*mBohrRadius);
  long double hfun = geth(eta);
  dcomplex gtilde  = GetG(eta, rho, hfun);
  dcomplex gtilor  = mult(gtilde, 1.0/mRSt);
  if (mPrint) {
    cout << "Gfun " << gtilor.re << " " << gtilor.im << endl;
    cout << "hfun " << hfun + 2*euler - 1.0 << endl;
  }

  dcomplex fcouls, fcoult;
  Getfc(kstar, eta, hfun, fcouls, fcoult);
  
  dcomplex fgs       = mult(gtilor, fcouls);
  long double fgmods = modl2(fgs);

  dcomplex fgt       = mult(gtilor, fcoult);
  long double fgmodt = modl2(fgt);

  dcomplex expikr;
  expikr.re = cos(tKstRst);
  expikr.im = -sin(tKstRst);
  dcomplex expikrc = conj(expikr);
  dcomplex ffplusc = conj(ffplus);
  dcomplex ffminusc = conj(ffminus);

  dcomplex expikr2 = mult(expikr, expikr);
  dcomplex expikrc2 = conj(expikr2);
  dcomplex sterm = mult(expikr2,  mult(ffplus, ffminusc));
  dcomplex tterm = mult(expikrc2, mult(ffminus, ffplusc));

  dcomplex epfpc = mult(expikr, ffplus);
  dcomplex emfmc = mult(expikrc, ffminus);

  long double fcgefhs = (fgs.re*emfmc.re + fgs.im*emfmc.im);
  long double fcgefgs = (fgs.re*epfpc.re + fgs.im*epfpc.im);

  long double fcgefht = (fgt.re*emfmc.re + fgt.im*emfmc.im);
  long double fcgefgt = (fgt.re*epfpc.re + fgt.im*epfpc.im);

  long double smult = 1+wavesign;

  if (!finite(ffminus.re))
    cout << "FFMinus Re not a number ! " << testp << " " << testm << " " << ccase<< endl;
  
  if (!finite(ffminus.im))
    cout << "FFMinus Im not a number !" << testp << " " << testm << " " << ccase<< endl;

  if ((ffplus.re > 2.0) || (ffplus.re < -2.0))
    cout << "FFplus Re wild !" << ffplus.re << endl;

  if ((ffplus.im > 2.0) || (ffplus.im < -2.0))
    cout << "FFplus Im wild !" << ffplus.im << endl;

  if ((ffminus.re > 2.0) || (ffminus.re < -2.0))
    cout << "FFminus Re wild !" << ffminus.re << " " << ccase << endl;
  
  if ((ffminus.im > 2.0) || (ffminus.im < -2.0))
    cout << "FFminus Im wild !" << ffminus.im << " " << ccase <<endl;

  //  coulqscpart = 0.5 * Gamow(kstar) * (modl2(ffplus) + modl2(ffminus));

  //   return (0.5 * Gamow(kstar) * 
  // 	  (modl2(ffplus) + wavesign*sterm.re + wavesign*tterm.re + modl2(ffminus)));
  if (twospin == 1) {
    wavesign = 1;
    smult = 2;
    long double singlet = (0.5 * Gamow(kstar) *
			   (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
			    wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefhs + fcgefgs)));
    wavesign = -1;
    smult = 0;
    long double triplet = (0.5 * Gamow(kstar) *
			   (2.0 * fgmodt * smult + modl2(ffplus) + modl2(ffminus) + 
			    wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefht + fcgefgt)));
    //    scmp = singlet;
    //    tcmp = triplet;

    //    ccmp = 0.5 * Gamow(kstar) * (modl2(ffplus) + modl2(ffminus));
    //    gcmp = fgmod;

    return (0.25 * singlet + 0.75 * triplet);
    //    return triplet;
  }
  else
    return (0.5 * Gamow(kstar) *
	    (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
	     wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefhs + fcgefgs)));
		      
}


void StandAloneFsiKisiel::PairKinematics(Pair &aPair)
{
  // Calculate pair variables
  double tPx = aPair.p1().x+aPair.p2().x;
  double tPy = aPair.p1().y+aPair.p2().y;
  double tPz = aPair.p1().z+aPair.p2().z;
  double tE  = aPair.p1().t+aPair.p2().t;
  double tPt = tPx*tPx + tPy*tPy;
  double tMt = tE*tE - tPz*tPz;//mCVK;
  double tM  = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);

  // Boost to LCMS
  double tBeta = tPz/tE;
  double tGamma = tE/tMt;	    

  double tDX = aPair.x1().x-aPair.x2().x;
  double tDY = aPair.x1().y-aPair.x2().y;
  double tRLong = aPair.x1().z-aPair.x2().z;
  double tDTime = aPair.x1().t-aPair.x2().t;
  
  double tROut = (tDX*tPx + tDY*tPy)/tPt;
  double tRSide = (-tDX*tPy + tDY*tPx)/tPt;
  double tRSidePairCMS = tRSide;
  mRSS = tRSidePairCMS/mFmToGeV;

  double tRLongPairCMS = tGamma*(tRLong - tBeta* tDTime);
  mRLS = tRLongPairCMS/mFmToGeV;
  double tDTimePairLCMS = tGamma*(tDTime - tBeta* tRLong);

  tBeta = tPt/tMt;
  tGamma = tMt/tM;
  double tROutPairCMS = tGamma*(tROut - tBeta* tDTimePairLCMS);
  mROS = tROutPairCMS/mFmToGeV;
  double tDTimePairCMS = tGamma*(tDTimePairLCMS - tBeta* tROut);
  mRTS = tDTimePairCMS/mFmToGeV;
  double tRStar = ::sqrt(tROutPairCMS*tROutPairCMS + tRSidePairCMS*tRSidePairCMS +
			 tRLongPairCMS*tRLongPairCMS);
  mRSt = tRStar/mFmToGeV;
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized)
// with no Coulomb interaction
double StandAloneFsiKisiel::GetQuantum(Pair &aPair)
{
  double tKstRst = aPair.GetKStarOut()*mROS + aPair.GetKStarSide()*mRSS + aPair.GetKStarLong()*mRLS;

  return (1.0 + cos(2.0*tKstRst));
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized)
// with Coulomb interaction included
double StandAloneFsiKisiel::GetQuantumCoulomb(Pair &aPair)
{
  double tKstRst = aPair.GetKStarOut()*mROS + aPair.GetKStarSide()*mRSS + aPair.GetKStarLong()*mRLS;

  if (mPrint == -1)
    mPrint = 1;
  else
    mPrint = 0;

  if (mPrint) {
    PRINT_MESSAGE("KstRst " << tKstRst);
    PRINT_MESSAGE("ROut RSide RLong DTime " << mROS << " " << mRSS << " " << mRLS << " " << mRTS);
  }
  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  if (fabs(mRSt) > fabs(mBohrRadius)) {
    if (twospin)
      return (1.0 - 0.5* cos(2*tKstRst));
    else 
      return (1.0 + cos(2*tKstRst));
  }
   
  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  if ((aPair.GetKStar() * mRSt * (1.0 + tKstRst/(mRSt*aPair.GetKStar()))> 15.0) &&
      (aPair.GetKStar() * mRSt * (1.0 - tKstRst/(mRSt*aPair.GetKStar()))> 15.0))
    {
      double fasymplus  = (1.0 - 1.0/((mRSt+tKstRst)*mBohrRadius*aPair.GetKStar()*aPair.GetKStar()));
      double fasymminus = (1.0 - 1.0/((mRSt-tKstRst)*mBohrRadius*aPair.GetKStar()*aPair.GetKStar()));
      if (mPrint)
	PRINT_MESSAGE("Double approximation");
      if (twospin) 
	return 0.5 * ((fasymplus + fasymminus) - (sqrt(fasymplus*fasymminus)*cos(2*tKstRst)));
      else
	return 0.5 * ((fasymplus + fasymminus) + (2.0*sqrt(fasymplus*fasymminus)*cos(2*tKstRst)));

    }
  //    return (1.0 - 1.0/(mRSt*mBohrRadius*aPair.GetKStar()*aPair.GetKStar()))*(1.0+cos(2*tKstRst));
  
  dcomplex ffplus, ffminus;
  // If RStar is zero, F is 1
  if (mRSt == 0) {
    if (twospin)
      return Gamow(aPair.GetKStar()) * (1.0 - 0.5*cos(2*tKstRst));
    else
      return Gamow(aPair.GetKStar()) * (1.0 + cos(2*tKstRst));
  }
  
  // Check for the classical limit in both functions separately
  if ((aPair.GetKStar() * mRSt * (1.0 + tKstRst/(mRSt*aPair.GetKStar()))< 15.0) &&
      (aPair.GetKStar() * mRSt * (1.0 - tKstRst/(mRSt*aPair.GetKStar()))< 15.0))
    {
      // Calculate the F function
      GetFFdouble(&ffplus, &ffminus, aPair);
      if (mPrint)
	PRINT_MESSAGE("Full calculation");
    }
  else if (aPair.GetKStar() * mRSt * (1.0 + tKstRst/(mRSt*aPair.GetKStar()))< 15.0)
    {
      double asym;
      GetFFsingle(&ffplus, 1, aPair);
      GetFFsingle(&ffminus, -1, aPair);
      if ((fabs(ffminus.re) > 10.0) || fabs(ffminus.im) > 10.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/(mRSt*aPair.GetKStar())*mBohrRadius*aPair.GetKStar()*aPair.GetKStar())))/Gamow(aPair.GetKStar());
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffminus.re = 1.0 + (asym -1.0) *2.0;
	else
	  ffminus.re = 1.0 + (asym -1.0) /2.0;
	ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);
      }
      if (mPrint)
	PRINT_MESSAGE("FFPlus full, FFMinus interpolated");
    }
  else 
    {
      double asym;
      GetFFsingle(&ffminus, -1, aPair);
      GetFFsingle(&ffplus, 1, aPair);
      if ((fabs(ffplus.re) > 10.0) || fabs(ffplus.im) > 10.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/(mRSt*aPair.GetKStar())*mBohrRadius*aPair.GetKStar()*aPair.GetKStar())))/Gamow(aPair.GetKStar());
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffplus.re = 1.0 + (asym - 1.0) * 2.0;
	else
	  ffplus.re = 1.0 + (asym - 1.0) / 2.0;
	ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      }
      if (mPrint) {
	PRINT_MESSAGE("FFMinus full, FFPlus interpolated");
	PRINT_MESSAGE("asym ffplus.re  ffplus.im  " << asym << " " << ffplus.re  << " " << ffplus.im);
	PRINT_MESSAGE("asym ffminus.re ffminus.im " << asym << " " << ffminus.re << " " << ffminus.im);
	// 	PRINT_MESSAGE("Rstar ROut RSide Rlong " << mRSt << " " << mROS << " " << mRSS << " " << mRLS);
	// 	PRINT_MESSAGE("r(1+cos(theta)) " << tKstRst/(mRSt*aPair.GetKStar()));
	// 	PRINT_MESSAGE("down " << mRSt*(1.0 + tKstRst/(mRSt*aPair.GetKStar())*mBohrRadius*aPair.GetKStar()*aPair.GetKStar()));
      }
    }

  dcomplex expikr;
  expikr.re = cos(tKstRst);
  expikr.im = -sin(tKstRst);
  dcomplex expikrc = conj(expikr);
  dcomplex ffplusc = conj(ffplus);
  dcomplex ffminusc = conj(ffminus);

  dcomplex expikr2 = mult(expikr, expikr);
  dcomplex expikrc2 = conj(expikr2);
  dcomplex sterm = mult(expikr2,  mult(ffplus, ffminusc));
  dcomplex tterm = mult(expikrc2, mult(ffminus, ffplusc));

  if (mPrint)
    {
      PRINT_MESSAGE("Gamow sqrt(Gamow) " << Gamow(aPair.GetKStar()) << " " << TMath::Sqrt(Gamow(aPair.GetKStar())));
      dcomplex f12 = mult(ffplus, expikr);
      dcomplex f21 = mult(ffminus, expikrc);
      f12 = mult(f12,  TMath::Sqrt(Gamow(aPair.GetKStar())));
      f21 = mult(f21,  TMath::Sqrt(Gamow(aPair.GetKStar())));

      PRINT_MESSAGE("cosh sinh " << expikr.re << " " << expikr.im);

      PRINT_MESSAGE("f12 " << f12.re << " " << f12.im);
      PRINT_MESSAGE("f21 " << f21.re << " " << f21.im);
      
      dcomplex psi;
      psi.re = f12.re + f21.re;
      psi.im = f12.im + f21.im;
      
      PRINT_MESSAGE("psi " << psi.re << " " << psi.im);
      PRINT_MESSAGE("weight " << (0.5 * modl2(psi)));
    }
  if (twospin) {
    return (0.5 * Gamow(aPair.GetKStar()) * 
	    (modl2(ffplus) + modl2(ffminus) - 0.5*(sterm.re + tterm.re)));
  }
  else {
    return (0.5 * Gamow(aPair.GetKStar()) * 
	    (modl2(ffplus) + sterm.re + tterm.re + modl2(ffminus)));;
  }
}

// Calculate the wave-function modulus sqaured
// for non-identical particles
// that is the Coulomb interaction
double StandAloneFsiKisiel::GetCoulomb(Pair &aPair)
{
  double tKstRst = aPair.GetKStarOut()*mROS + aPair.GetKStarSide()*mRSS + aPair.GetKStarLong()*mRLS;

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  if (fabs(mRSt) > fabs(mBohrRadius)) return 1.0;

  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  if (aPair.GetKStar() * mRSt > 15.0)
    return (1.0 - 1.0/(mRSt*mBohrRadius*aPair.GetKStar()*aPair.GetKStar()));
  
  // Calculate the F function
  dcomplex ffplus;
  GetFFsingle(&ffplus, 1, aPair);

  return (Gamow(aPair.GetKStar()) * modl2(ffplus));
}

void StandAloneFsiKisiel::setPairType(const int aPairType){ 
  mLL = aPairType;
  long double tReducedMass;
  int tSign;

  d0s.re = 0.0;
  d0s.im = 0.0;
  f0s.re = 0.0;
  f0s.im = 0.0;
  d0t.re = 0.0;
  d0t.im = 0.0;
  f0t.re = 0.0;
  f0t.im = 0.0;
  twospin = 0;
  euler = 0.577215665;

  switch (aPairType) {
  case PiPlusPiPlus:
    tReducedMass = mPi/2.0;
    tSign = 1;
    break;
  case PiPlusPiMinus:
    tReducedMass = mPi/2.0;
    tSign = -1;
    break;
  case KaonPlusKaonPlus:
    tReducedMass = mK/2.0;
    tSign = 1;
    break;
  case KaonPlusKaonMinus:
    tReducedMass = mK/2.0;
    tSign = -1;
    break;
  case ProtonProton:
    tReducedMass = mP/2.0;
    tSign = 1;
    d0s.re = 2.77 / mFmToGeV;
    f0s.re = 7.77 / mFmToGeV;
    d0t.re = 1.7  / mFmToGeV;
    f0t.re = -5.4 / mFmToGeV;
    twospin = 1;
    break;
  case PiPlusProton:
    tReducedMass = mP*mPi/(mP+mPi);
    tSign = 1;
    break;
  case PiMinusProton:
    tReducedMass = mP*mPi/(mP+mPi);
    tSign = -1;
    break;
  case PiPlusKaonMinus:
    tReducedMass = mK*mPi/(mK+mPi);
    tSign = -1;
    break;
  case PiPlusKaonPlus:
    tReducedMass = mK*mPi/(mK+mPi);
    tSign = 1;
    break;
  case KaonPlusProton:
    tReducedMass = mK*mP/(mK+mP);
    tSign = 1;
    break;
  case KaonMinusProton:
    tReducedMass = mK*mP/(mK+mP);
    tSign = -1;
    break;
  default:
    PRINT_MESSAGE("Adam Kisiel's code does not support that pair type: " << aPairType);
    break;
  }
  
  // Calculate system Bohr radius in 1/GeV
  mBohrRadius = tSign * mAlfa/ tReducedMass;

  // Initialize Gamow factor table
  InitializeGamow();
}

void StandAloneFsiKisiel::setDefaultCalcPar(){ 
  mIch=1;
  mIqs=0;
  mIsi=0;
  mI3c=0;
}

void StandAloneFsiKisiel::setCoulOn(){ 
  mIch = 1;
}

void StandAloneFsiKisiel::setCoulOff(){ 
  mIch = 0;
}

void StandAloneFsiKisiel::setQuantumOn(){
  mIqs = 1;
}

void StandAloneFsiKisiel::setQuantumOff(){
  mIqs = 0;
}

void StandAloneFsiKisiel::setStrongOn(){
  //  PRINT_MESSAGE("The weight calculator of Adam Kisiel code does not support strong interaction calculation yet.");
  mIsi = 1;
}

void StandAloneFsiKisiel::setStrongOff(){
  mIsi = 0;
}

void StandAloneFsiKisiel::set3BodyOn(){
  mI3c = 1;
}

void StandAloneFsiKisiel::set3BodyOff(){
  mI3c = 0;
}

void StandAloneFsiKisiel::ReadParameters()
{
  /* no-op */
}

void StandAloneFsiKisiel::GenerateParameterStub()
{
  /* no-op */
}

void StandAloneFsiKisiel::PrintNext()
{
  mPrint = -1;
}
