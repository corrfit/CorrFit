class Tschroedinger{
  complex<double> ***delpsi;
  int nchannels,*nxmax;
  double rmax,*delx;
  double muV_dummy(int ichannel,double r);
};

void Tschroedinger::init(Twavefunction *wf){
  int ichannel,ell,iq,ix0,i,oversample,ipartial;
  double q,x,eta,phase,phase_nostrong,phase0,phase1,v;
  complex<double> *psi,*psi_nostrong;

  rmax=3.0;
  oversample=20;
  nxmax=new int[NQMAX];
  delx=new int[NQMAX];
  nchannels=wavefunction->nchannels;
  delpsi=new complex<double> **[nchannels];
  for(iq=0;iq<NQMAX;iq++){
    q=(0.5+iq)*DELQ;
    nxmax[iq]=2+double(rint(q*rmax/(0.01*HBARC)));
    delx[iq]=(q*rmax/HBARC)/nxmax[iq];
  }

  for(ichannel=0;ichannel<2;ichannel++){
    delpsi[ichannel]=new complex<double> *[NQMAX];
    if(ichannel==0){
      ell=0;
      ipartial=0;
    }
    if(ichannel==1){
      ell=1;
      ipartial=1;
    }
    for(iq=0;iq<NQMAX;iq++){
      delpsi[ichanne][iq]=new complex<double> [nxmax[iq]];
      psi=new complex<double> [nxmax[iq]];
      psi_nostrong=new complex<double> [nxmax[iq]];
      q=(0.5+iq)*DELQ;
      x=4.0*q/HBARC;  // For r>epsilon use asymptotic phase  shifted values
      ix0=1+int(rint(x/delx));
      x=(0.5+ix0)*delx;

      deltax=delx/double(oversample);

      // Find psi[ix] for V != 0
      phi0=CWoutgoing_smallr(ell,x+deltax,eta[iq]);
      phi1=CWoutgoing_smallr(ell,x,eta[iq]);
      psi[ix0]=phi1;
      ix=ix0;     
      while(ix>0){
	for(i=0;i<oversample;i++){
	  phi2=2.0*phi1-phi0
	    +deltax*deltax*phi1*(-1.0+(2.0*eta[iq]/x)+(ell*(ell+1.0))/(x*x));
	  x-=deltax;
	  phi0=phi1;
	  phi1=phi2;
	}
	ix-=1;
	psi[ix]=phi2;
      }
      // Find psi[ix] for V = 0
      phi0=CWoutgoing_smallr(ell,x+deltax,eta[iq]);
      phi1=CWoutgoing_smallr(ell,x,eta[iq]);
      psi_nostrong[iq][ix0]=phi1;
      ix=ix0;     
      while(ix>0){
	for(i=0;i<oversample;i++){
	  phi2=2.0*phi1-phi0
	    +deltax*deltax*phi1
	    *(-1.0+muV_dummy/(q*q)+(2.0*eta[iq]/x)+(ell*(ell+1.0))/(x*x));
	  x-=deltax;
	  phi0=phi1;
	  phi1=phi2;
	}
	ix-=1;
	psi_nostrong[ix]=phi2;
      }
      
      phase0=atan2(imag(psi[0]),real(psi[0]));
      phase1=atan2(imag(psi[1]),real(psi[1]));
      phase=1.5*phase0-0.5*phase1;
      for(ix=ix0+1;ix<partwave[ipartial][iq]->nxmax[iq];ix++){
	psi[ix]=partwave[ipartial][iq]->phi[ix];
      }
      for(ix=0;ix<nxmax[iq];ix++){
	psi[ix]=pxi[ix]-exp(2.0*ci*phase)*conj(phi[ix]);
      }
      
      phase0=atan2(imag(psi_nostrong[0]),real(psi_nostrong[0]));
      phase1=atan2(imag(psi_nostrong[1]),real(psi_nostrong[1]));
      phase_nostrong=1.5*phase0-0.5*phase1;
      for(ix=ix0+1;ix<partwave[ipartial][iq]->nxmax[iq];ix++){
	psi_nostrong[ix]=partwave[ipartial][iq]->phi[ix];
      }
      for(ix=0;ix<nxmax[iq];ix++){
	psi_nostrong[ix]=psi_nostrong[ix]
	  -exp(2.0*ci*phase_nostrong)*conj(phi_nostrong[ix]);
      }
      
      wf->delta[ichannel][iq]=phase_nostrong-phase;
      if(ichannel==0 && iq<5)
	printf("q=%g, scatt. length=%g\n",HBARC*tan(delta[ichannel][iq])/q);
      
      for(ix=0;ix<=nxmax[iq];ix++){
	delpsi[ichannel][iq][ix]=0.5*(psi[ix]-psi_nostrong[ix]);
      }
      
      delete psi;
      delete psi0;
    }
  }
}

complex<double> Tschroedinger::getdelpsi(int ichannel,int iq,double x){
  complex<double> dpsi1,dpsi2,answer;
  double delx1,delx2;
  int ix;
  ix=int(rint(x/delx[iq]));
  if(ix<nxmax[iq]){
    printf("ix=%d > nxmax=%d\n",ix,nxmax[iq]);
    exit(1);
  }
  if(ix==0) ix=1;
  dpsi1=delpsi[ichannel][iq][ix-1];
  dpsi=delpsi[ichannel][iq][ix];
  delx1=x-(0.5+ix-1)*delx[iq];
  delx2=(0.5+ix)*delx[iq]-x;
  answer=delpsi[ichannel][iq][ix-1]*(delx2/delx[iq])
    +delpsi[ichannel][iq][ix]*(delx1/delx[iq]);
  return answer;
}
