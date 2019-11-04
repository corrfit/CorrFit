Tplanewave::Tplanewave(double etaset,int Q1Q2,double qset){
  complex<double> a,b,a0,b0;
  complex<double> c1top1,c2top2,bot;
  complex<double> c1top2,c2top1;
  int j;
  PI=4.0*atan(1.0);
  ci=complex<double>(0.0,1.0);

  eta=etaset;
  q1q2=Q1Q2;
  q=qset;
  if(q1q2!=0){
    couly=cgamma(1.0+ci*eta);
    couly=couly*exp(-0.5*eta*PI);
    a=-ci*eta;
    b=1.0;
    a0=a;
    b0=b;
    for(j=1;j<=40;j++){
      chype[j]=a/(b*double(j));
      a=a+1.0;
      b=b+1.0;
    }
    c1top1=1.0;
    c2top1=1.0;
    c1top2=1.0;
    c2top2=1.0;
    bot=1.0;
    for(j=1;j<=5;j++){
      c1top1=c1top1*(double(j)+a0-1.0);
      c2top1=c2top1*(double(j)-a0);
      c1top2=c1top2*(double(j)-b0+a0);
      c2top2=c2top2*(double(j)+b0-a0-1.0);
      bot=bot*double(j);
      chype1[j]=(c1top1*c1top2)/(bot);
      chype2[j]=(c2top1*c2top2)/(bot);
    }
    cfact1=1.0/(cgamma(b0-a0));
    cfact2=1.0/(cgamma(a0));
  }
}

complex<double> Tplanewave::planewave(double r,double ctheta){
  const double HBARC=197.3269602;
  complex<double> bb,answer;
  double zq,arg;
  /* See appendix of Messiah.  "cgamma" is the gamma function. "hyper" is
     appropriate hyperbolic function.  This is explained in Messiah's
     section on coul wave func.s. Notation should be explanatory
     when reading the book. 
     This has been modified so that the outgoing waves correspond to a plane
     wave rather than the incoming waves */
  if(q1q2!=0){
    bb=1.0;
    zq=-r*ctheta; // (By flipping ctheta, then taking c.c. below, we get
    // wave that leavew with momentum q, rather than enters with q
    answer=couly*hyper(-ci*eta,bb,ci*q*(r-zq)/HBARC);
    arg=zq*q/HBARC;
    arg=arg-2.0*PI*floor(arg/(2.0*PI));
    answer=answer*exp(ci*arg);
    answer=conj(answer);
  }
  else{
    answer=exp(ci*q*r*ctheta/HBARC);
  }

  return answer;
}

/* **************************************************** */


complex<double> Tplanewave::hyper(complex<double> a,complex<double> b,complex<double> cz){
  complex<double> cw1,cw2,cf1,cf2,czarg,czstarj,answer;
#define rcrit 10.0
  double dmag;
  int j;
  dmag=abs(cz);
  if(dmag<rcrit){
    cf1=1.0;
    czstarj=1.0;
    for(j=1;j<=40;j++){
      czstarj=czstarj*cz*chype[j];
      cf1=cf1+czstarj;
      if(abs(czstarj)<0.0001) goto GOOD_ENOUGH;
    }
    printf("hyper not coverging!.\n");
  GOOD_ENOUGH:
    answer=cf1;
  }
  else{
    cf1=1.0;
    cf2=1.0;
    for(j=1;j<=5;j++){
      cf1=cf1+chype1[j]/(pow(-cz,j));
      cf2=cf2+chype2[j]/(pow(cz,j));
    }
    cw1=cf1*cfact1*(pow(-cz,-a));
    czarg=cz-ci*2.0*PI*floor(real(-ci*cz/(2.0*PI)));
    cw2=cf2*cfact2*pow(cz,a-b)*exp(czarg);
    answer=cw1+cw2;
  }
  return answer;
}
