Tpartwave::Tpartwave(double etaset,int q1q2set,double qset,int ellset){
  complex<double> cg;
  PI=4.0*atan(1.0);
  ci=complex<double>(0.0,1.0);
  q=qset;
  q1q2=q1q2set;
  eta=etaset;
  ell=ellset;
  sigma=0.0;
  if(q1q2!=0){
    cg=cgamma(ell+1.0+ci*eta);
    sigma=atan2(imag(cg),real(cg));
  }
  phi_init();
}

void Tpartwave::phi_init(){
  complex<double> phi0,phi1,phi2,phitest;
  int i,ix,ix0,oversample;
  double deltax,x;

  delx=0.01;
  nxmax=int(floor(12.0/delx));
  phi=new complex<double> [nxmax+1];

  if(q1q2!=0){
    oversample=20;
    deltax=delx/double(oversample);
    ix0=int(floor(0.5/delx));
    for(ix=0;ix<=ix0;ix++){
      x=(0.5+ix)*delx;
      phi[ix]=CWincoming_smallr(ell,x,eta);
    }
    
    x=(0.5+ix0)*delx;
    ix=ix0;
    phi0=CWincoming_smallr(ell,x-deltax,eta);
    phi1=CWincoming_smallr(ell,x,eta);
    
    while(ix<nxmax){
      for(i=0;i<oversample;i++){
	phi2=2.0*phi1-phi0
	  +deltax*deltax*phi1*(-1.0+(2.0*eta/x)+(ell*(ell+1.0))/(x*x));
	x+=deltax;
	phi0=phi1;
	phi1=phi2;
      }
      ix+=1;
      phi[ix]=phi2;
    }
  }
  else{
    for(ix=0;ix<nxmax;ix++){
      x=(0.5+ix)*delx;
      phi[ix]=hankelstar(ell,x);
    }
  }

}

complex<double> Tpartwave::GetPhiIncoming(double r){
  const double HBARC=197.3269602;
  double x,a;
  complex<double> answer;
  int ix;
  x=q*r/HBARC;
  ix=int(floor(x/delx));
  if(x>0.01 && ix<nxmax){
    a=(x-delx*(ix+0.5))/delx;
    if(a>0.0 || ix==0){
      answer=(1.0-a)*phi[ix]+a*phi[ix+1];
    }
    else{
      answer=(1.0+a)*phi[ix]-a*phi[ix-1];
    }
  }
  else if(ix>=nxmax){
    if(q1q2!=0)
      answer=CWincoming_bigr(ell,x,eta);
    else
      answer=hankelstar(ell,x);
  }
  else{
    if(q1q2!=0)
      answer=CWincoming_smallr(ell,x,eta);
    else
      answer=hankelstar(ell,x);
  }
  answer=answer*exp(-ci*sigma);
  return answer;
}
