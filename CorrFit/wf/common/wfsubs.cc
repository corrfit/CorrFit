complex<double> cgamma(complex<double> c){
  /* This calc.s gamma functions which are in the form gamma(n+i*y)
     where n is an int and y is real. */
  complex<double> cg,cphase,ci(0.0,1.0);
  const double PI=4.0*atan(1.0);
  double rz[8]={0.0,0.0,1.64493406684822641,1.20205690315031610,
		1.08232323371113792,1.03692775514333801,1.01734306198444879,
		1.00834927738192071};
  double rzleftover[8];
  int mm,j,ir;
  double x,y,yj,phase,delp,cgmag,EULER=0.5772156649015328606;
  x=real(c);
  y=imag(c);
  phase=-EULER*y;
  for(ir=3;ir<=7;ir++) rzleftover[ir]=rz[ir];
  j=0;
  do{
    j=j+1;
    for(ir=3;ir<=7;ir++) rzleftover[ir]-=1.0/pow(double(j),ir);
    yj=y/double(j);
    delp=(yj)-atan(yj);
    phase=phase+delp;
  }while(j<4 || fabs(delp)>1.0E-4);
  for(ir=3;ir<=7;ir+=2) 
  delp=pow(y,3)*rzleftover[3]/3.0;
  phase+=delp;
  delp=-pow(y,5)*rzleftover[5]/5.0;
  phase+=delp;
  delp=pow(y,7)*rzleftover[7]/7.0;
  phase+=delp;
CGAMMA_ESCAPE:
  phase=phase-2.0*PI*floor(phase/(2.0*PI));
  cphase=exp(ci*phase);
  cgmag=sqrt(PI*y/sinh(PI*y));
  mm=(int)floor(x+0.5);
  cg=cgmag*cphase;
  if(mm<1){
    for(j=1;j<=-mm+1;j++){
      cg=cg/(1.0+double(-j)+ci*y);
    }
  }
  if(mm>1) {
    for(j=1;j<=mm-1;j++){
      cg=cg*(double(j)+ci*y);
    }
  }
  return cg;
}

complex<double> CWincoming_bigr(int ell,double x,double eta){
  // See Abramowitz and Stegun, page 540, Eq. 14.5.8
  const double PI=4.0*atan(1.0);
  double  fk,gk,ak,bk,f,g,ftemp;
  double arg,sigma;
  complex<double> cg,answer,ci(0.0,1.0);
  int k,kmax;
  cg=cgamma(ell+1.0+ci*eta);
  sigma=atan2(imag(cg),real(cg));
  f=fk=1.0;
  g=gk=0.0;
  kmax=1+int(floor(x));
  if(kmax>7) kmax=7;
  for(k=0;k<kmax;k++){
    ak=(2.0*k+1.0)*eta/(2.0*x*(k+1.0));
    bk=(double(ell*(ell+1)-k*(k+1))+eta*eta)/(2.0*x*(k+1.0));
    ftemp=fk;
    fk=ak*fk-bk*gk;
    gk=ak*gk+bk*ftemp;
    f=f+fk;
    g=g+gk;
  }
  
  arg=x-eta*log(2.0*x)-0.5*double(ell+1)*PI+sigma;
  arg=arg-2.0*PI*floor(arg/(2.0*PI));
  answer=(f+ci*g)*exp(ci*arg);
  return conj(answer);
  } 
/* ************************************* */

complex<double> CWoutgoing_bigr(int ell,double x,double eta){
  // See Abramowitz and Stegun, page 540, Eq. 14.5.8
  const double PI=4.0*atan(1.0);
  double  fk,gk,ak,bk,f,g,ftemp;
  double arg,sigma;
  complex<double> cg,answer,ci(0.0,1.0);
  int k,kmax;
  cg=cgamma(ell+1.0+ci*eta);
  sigma=atan2(imag(cg),real(cg));
  f=fk=1.0;
  g=gk=0.0;
  kmax=1+int(floor(x));
  if(kmax>7) kmax=7;
  for(k=0;k<kmax;k++){
    ak=(2.0*k+1.0)*eta/(2.0*x*(k+1.0));
    bk=(double(ell*(ell+1)-k*(k+1))+eta*eta)/(2.0*x*(k+1.0));
    ftemp=fk;
    fk=ak*fk-bk*gk;
    gk=ak*gk+bk*ftemp;
    f=f+fk;
    g=g+gk;
  }
  
  arg=x-eta*log(2.0*x)-0.5*double(ell+1)*PI+sigma;
  arg=arg-2.0*PI*floor(arg/(2.0*PI));
  answer=(f+ci*g)*exp(ci*arg);
  return answer;
  } 
/* ************************************* */



complex<double> CWoutgoing_smallr(int ell,double x,double eta){
  /* The notation is like Gr. + R, page 1063.
     The Coulomb wave function is the same as W(i*eta,l+1/2,2*i*rho) */
  const double PI=4.0*atan(1.0);
  complex<double> ci(0.0,1.0),factor,lterm,fact1,fact2;
  complex<double> psi1,psi2,psi3,cx,sum1,sum2,delsum1;
  complex<double> delp,cdcon1,cdcon2,cg,answer;
  double sigma,EULER=0.5772156649015328606;
  int k;

  cg=cgamma(ell+1.0+ci*eta);
  sigma=atan2(imag(cg),real(cg));

  cdcon1=cgamma(-double(ell)-ci*eta);
  cdcon2=cgamma(double(ell+1)-ci*eta);
  factor=pow(-1.0,double(2*ell+1))*pow(2.0*ci*x,double(ell+1))
    *exp(-ci*x)/(cdcon1*cdcon2);
  psi1=-EULER;
  psi2=-EULER;
  for(k=1;k<=2*ell+1;k++){
    psi2=psi2+1.0/double(k);
  }
  cx=double(ell+1)-ci*eta;
  psi3=-EULER-(1.0/cx)+cx*(PI*PI/6.0);
  for(k=1;k<100000;k++){
    delp=-cx*cx/(double(k*k)*(cx+double(k)));
    psi3=psi3+delp;
    if(abs(delp)<1.0E-12) goto CONVERGE1;
  }
  printf("never escaped loop1 in CWoutgoing_smallr!\n");
CONVERGE1:
  lterm=log(2.0*ci*x);
  fact1=cdcon2/dgamma(2*ell+2);
  sum1=fact1*(psi1+psi2-psi3-lterm);
  for(k=1;k<=10000;k++){
    fact1=fact1*(2.0*ci*x)
      *(double(ell+k)-ci*eta)/double(k*(2*ell+1+k));
    psi1=psi1+1.0/double(k);
    psi2=psi2+1.0/double(2*ell+1+k);
    psi3=psi3+1.0/(double(k-1)+cx);
    delsum1=fact1*(psi1+psi2-psi3-lterm);
    sum1=sum1+delsum1;
    if(abs(delsum1)<1.0E-15) goto CONVERGE2;
  }
  printf("never escaped loop2 in CWoutgoing_smallr!\n");
CONVERGE2:
  fact2=dgamma(2*ell+1)*cdcon1/pow(-2.0*ci*x,2*ell+1);
  sum2=fact2;
  for(k=1;k<=2*ell;k++){
    fact2=fact2*(double(k-ell-1)-ci*eta)*
      (-2.0*ci*x)/(double(k)*double(2*ell-k+1));
    sum2=sum2+fact2;
  }
  sum1=factor*sum1;
  sum2=factor*sum2;
  answer=(sum1+sum2)*exp(PI*eta/2.0);
  answer=exp(-0.5*ci*double(ell+1)*PI+ci*sigma)*conj(answer);
  return answer;
}

/* ***************************************** */

complex<double> CWincoming_smallr(int ell,double x,double eta){
  /* The notation is like Gr. + R, page 1063.
     The Coulomb wave function is the same as W(i*eta,l+1/2,2*i*rho) */
  const double PI=4.0*atan(1.0);
  complex<double> ci(0.0,1.0),factor,lterm,fact1,fact2;
  complex<double> psi1,psi2,psi3,cx,sum1,sum2,delsum1;
  complex<double> delp,cdcon1,cdcon2,cg,answer;
  double sigma,EULER=0.5772156649015328606;
  int k;

  cg=cgamma(ell+1.0+ci*eta);
  sigma=atan2(imag(cg),real(cg));

  cdcon1=cgamma(-double(ell)-ci*eta);
  cdcon2=cgamma(double(ell+1)-ci*eta);
  factor=pow(-1.0,double(2*ell+1))*pow(2.0*ci*x,double(ell+1))
    *exp(-ci*x)/(cdcon1*cdcon2);
  psi1=-EULER;
  psi2=-EULER;
  for(k=1;k<=2*ell+1;k++){
    psi2=psi2+1.0/double(k);
  }
  cx=double(ell+1)-ci*eta;
  psi3=-EULER-(1.0/cx)+cx*(PI*PI/6.0);
  for(k=1;k<100000;k++){
    delp=-cx*cx/(double(k*k)*(cx+double(k)));
    psi3=psi3+delp;
    if(abs(delp)<1.0E-12) goto CONVERGE1;
  }
  printf("never escaped loop1 in CWoutgoing_smallr!\n");
CONVERGE1:
  lterm=log(2.0*ci*x);
  fact1=cdcon2/dgamma(2*ell+2);
  sum1=fact1*(psi1+psi2-psi3-lterm);
  for(k=1;k<=10000;k++){
    fact1=fact1*(2.0*ci*x)
      *(double(ell+k)-ci*eta)/double(k*(2*ell+1+k));
    psi1=psi1+1.0/double(k);
    psi2=psi2+1.0/double(2*ell+1+k);
    psi3=psi3+1.0/(double(k-1)+cx);
    delsum1=fact1*(psi1+psi2-psi3-lterm);
    sum1=sum1+delsum1;
    if(abs(delsum1)<1.0E-15) goto CONVERGE2;
  }
  printf("never escaped loop2 in CWoutgoing_smallr!\n");
CONVERGE2:
  fact2=dgamma(2*ell+1)*cdcon1/pow(-2.0*ci*x,2*ell+1);
  sum2=fact2;
  for(k=1;k<=2*ell;k++){
    fact2=fact2*(double(k-ell-1)-ci*eta)*
      (-2.0*ci*x)/(double(k)*double(2*ell-k+1));
    sum2=sum2+fact2;
  }
  sum1=factor*sum1;
  sum2=factor*sum2;
  answer=(sum1+sum2)*exp(PI*eta/2.0);
  answer=exp(-0.5*ci*double(ell+1)*PI+ci*sigma)*conj(answer);
  return conj(answer);
}

/* ***************************************** */

complex<double> hankel(int ell,double x){
  complex<double> ci(0.0,1.0),answer;
  double j0,j1,j2,n0,n1,n2;
  int i;
  if(ell==0){
    j0=sin(x);
    n0=-cos(x);
    answer=(j0+ci*n0);
  }
  else if(ell==1){
    j1=(sin(x)/x)-cos(x);
    n1=-(cos(x)/x)-sin(x);
    answer=(j1+ci*n1);
  }
  else{
    j0=sin(x);
    n0=-cos(x);
    j1=(sin(x)/x)-cos(x);
    n1=-(cos(x)/x)-sin(x);
    for(i=2;i<=ell;i++){
      j2=-j0+(2*double(i)-1.0)*j1/x;
      n2=-n0+(2*double(i)-1.0)*n1/x;
      j0=j1; j1=j2;
      n0=n1; n1=n2;
    }
    answer=(j2+ci*n2);
  }
  return answer;
}

complex<double> hankelstar(int ell,double x){
  return conj(hankel(ell,x));
}

double GetIW(int ell,double epsilon,double q,int q1q2,double eta,
	     double delta){
  complex<double> ci(0.0,1.0),psi0,psi,psi2,psiminus0,psiminus,psiminus2;
  complex<double> psiplus0,psiplus,psiplus2;
  complex<double> ddeta_psi,ddeta_psiminus,ddeta_psiplus;
  complex<double> I,e2idelta;
  double x,deleta,a2,root;
  const double HBARC=197.3269602;

  if(q1q2!=0){
    deleta=0.001*eta;
    x=q*epsilon/HBARC;
    e2idelta=exp(2.0*ci*delta);
    if(ell==0){
      a2=1.0+eta*eta;
      root=sqrt(a2);

      psi0=CWoutgoing_smallr(ell,x,eta-0.5*deleta);
      psi0=psi0*e2idelta+conj(psi0);

      psiplus0=CWoutgoing_smallr(ell+1,x,eta-0.5*deleta);
      psiplus0=psiplus0*e2idelta+conj(psiplus0);

      psi=CWoutgoing_smallr(ell,x,eta);
      psi=psi*e2idelta+conj(psi);

      psiplus=CWoutgoing_smallr(ell+1,x,eta);
      psiplus=psiplus*e2idelta+conj(psiplus);

      psi2=CWoutgoing_smallr(ell,x,eta+0.5*deleta);
      psi2=psi2*e2idelta+conj(psi2);

      psiplus2=CWoutgoing_smallr(ell+1,x,eta+0.5*deleta);
      psiplus2=psiplus2*e2idelta+conj(psiplus2);

      ddeta_psi=(psi2-psi0)/deleta;
      ddeta_psiplus=(psiplus2-psiplus0)/deleta;

      I=(conj(psi)*psi+conj(psiplus)*psiplus)*x*a2;
      I-=conj(psi)*psiplus*(a2*(1.0+2.0*eta*x)+eta*eta)/root;
      I+=(conj(psiplus)*ddeta_psi-conj(psi)*ddeta_psiplus)*eta*root;
    }
    else{
      a2=double(ell*ell)+eta*eta;
      root=sqrt(a2);

      psi0=CWoutgoing_smallr(ell,x,eta-0.5*deleta);
      psi0=psi0*e2idelta+conj(psi0);

      psiminus0=CWoutgoing_smallr(ell-1,x,eta-0.5*deleta);
      psiminus0=psiminus0*e2idelta+conj(psiminus0);

      psi=CWoutgoing_smallr(ell,x,eta);
      psi=psi*e2idelta+conj(psi);

      psiminus=CWoutgoing_smallr(ell-1,x,eta);
      psiminus=psiminus*e2idelta+conj(psiminus);

      psi2=CWoutgoing_smallr(ell,x,eta+0.5*deleta);
      psi2=psi2*e2idelta+conj(psi2);

      psiminus2=CWoutgoing_smallr(ell-1,x,eta+0.5*deleta);
      psiminus2=psiminus2*e2idelta+conj(psiminus2);

      ddeta_psi=(psi2-psi0)/deleta;
      ddeta_psiminus=(psiminus2-psiminus0)/deleta;

      ddeta_psi=(psi2-psi0)/deleta;
      ddeta_psiminus=(psiminus2-psiminus0)/deleta;

      I=(conj(psi)*psi+conj(psiminus)*psiminus)*a2*x/double(ell*ell);
      I-=conj(psiminus)*psi
	*((2.0*ell+1.0)*a2*double(ell)+2.0*eta*x*a2
	  -eta*eta*double(ell))/(double(ell*ell)*root);
      I+=(conj(ddeta_psiminus)*psi-conj(ddeta_psi)*psiminus)
	*eta*root/double(ell);

    }
    I=0.25*I;
  }
  else{
    x=q*epsilon/HBARC;
    if(ell==0){
      psi=sin(x+delta);
      psiplus=-cos(x+delta)+sin(x+delta)/x;
    
      I=(conj(psi)*psi+conj(psiplus)*psiplus)*x;
      I-=conj(psi)*psiplus;
    }
    else{
      psi=hankel(ell,x);
      psi=(exp(2.0*ci*delta)*psi+conj(psi))/(2.0);
      psiminus=hankel(ell-1,x);
      psiminus=(exp(2.0*ci*delta)*psiminus+conj(psiminus))/(2.0);

      I=(conj(psi)*psi+conj(psiminus)*psiminus)*x;
      I-=(2.0*double(ell)+1)*conj(psiminus)*psi;
    
    }
  }
  return -real(I)/q;
}

void phaseshift_CoulombCorrect(int ell,double q,double eta,
			       double &delta,double &ddeltadq){
  const double PI=4.0*atan(1.0);
  double *gamow,*dgamowdq;
  double tandelta0,tandelta,dtandelta0dq,dtandeltadq,x,y;
  int i;

  if(fabs(eta)>1.0E-10){
    gamow=new double[ell+1];
    dgamowdq=new double[ell+1];
    gamow[0]=2.0*PI*eta/(exp(2.0*PI*eta)-1.0);
    dgamowdq[0]=(gamow[0]/q)*(gamow[0]*exp(2.0*PI*eta)-1.0);
    for(i=0;i<ell;i++){
      gamow[i+1]=gamow[i]*(double((i+1)*(i+1))+eta*eta);
      dgamowdq[i+1]=dgamowdq[i]*(double((i+1)*(i+1))+eta*eta)
	-(2.0*eta*eta/q)*gamow[i];
    }
    
    tandelta0=tan(delta);
    dtandelta0dq=ddeltadq*(1.0+tandelta0*tandelta0);
    
    tandelta=tandelta0*gamow[ell];
    if(cos(delta)>0.0) delta=atan(tandelta);
    else delta=atan(tandelta)+PI;
    if(delta>PI) delta=delta-2.0*PI;
    
    dtandeltadq=gamow[ell]*dtandelta0dq+tandelta0*dgamowdq[ell];
    ddeltadq=dtandeltadq/(1.0+tandelta*tandelta);
    delete gamow;
    delete dgamowdq;
  }
}
