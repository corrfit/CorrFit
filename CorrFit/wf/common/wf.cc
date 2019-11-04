Twavefunction::Twavefunction(){
  generic=0;
  PI=3.14159265358979323844;
  ROOT2=sqrt(2.0);
  ci=complex<double>(0.0,1.0);
  HBARC=197.3269602;
  MPI=139.58;
  MKAON=493.677;
  MPROTON=938.271;
  MLAMBDA=1115.7;
}

void Twavefunction::ParsInit(char *parsfilename){
  FILE *fptr,*qarrayfile;
  int iq;
  bool filetest=0;
  char value[100],variablename[100],dummy[100];
  char teststring[100],qarrayfilename[100];
  // DEFAULT VALUES
  nqmax=10;
  delq=10.0;
  epsilon=1.0;
  ellmax=2;
  COULOMB=1;
  STRONG=1;
  strcpy(qarrayfilename,"\0");
  // READ IN OTHER VALUES, file should have lines which have "value   NAME"
  // NAME needs to be uppercase, e.g. DELQ, NQMAX, EPSILON or ELLMAX
//   fptr=fopen(parsfilename,"r");
//   while(feof(fptr)==0){
//     strcpy(variablename,"\0");
//     fscanf(fptr,"%s %s",&value,&variablename);
//     fgets(dummy,100,fptr);
//     strcpy(teststring,"DELQ\0");
//     if(comparestrings(variablename,teststring)==1){
//       delq=atof(value);
//       printf("delq set to %g\n",delq);
//     }
//     strcpy(teststring,"NQMAX\0");
//     if(comparestrings(variablename,teststring)==1){
//       nqmax=atoi(value);
//       printf("nqmax set to %d\n",nqmax);
//     }
//     strcpy(teststring,"EPSILON\0");
//     if(comparestrings(variablename,teststring)==1){
//       epsilon=atof(value);
//       printf("epsilon set to %g\n",epsilon);
//     }
//     strcpy(teststring,"ELLMAX\0");
//     if(comparestrings(variablename,teststring)==1){
//       ellmax=atoi(value);
//       printf("ellmax set to %d\n",ellmax);
//     }
//     strcpy(teststring,"STRONG\0");
//     if(comparestrings(variablename,teststring)==1){
//       STRONG=atoi(value);
//       printf("STRONG set to %d\n",STRONG);
//     }
//     strcpy(teststring,"COULOMB\0");
//     if(comparestrings(variablename,teststring)==1){
//       COULOMB=atoi(value);
//       printf("COULOMB set to %d\n",COULOMB);
//     }
//     strcpy(teststring,"QARRAYFILENAME\0");
//     if(comparestrings(variablename,teststring)==1){
//       filetest=1;
//       strcpy(qarrayfilename,value);
//     }
//   }
//   fclose(fptr);
  nqmax=gPrattNQmax;
  delq=gPrattDelQ;
  epsilon=gPrattEpsilon;
  ellmax=gPrattEllMax;
  COULOMB=gPrattCoulomb;
  STRONG=gPrattStrong;
  if (gPrattQArrayFilename!="") {
    strcpy(qarrayfilename,gPrattQArrayFilename.Data());
    filetest=1;
  }
  else
    filetest = 0;
  
  if(filetest==1){
    delq=-1.0;
    printf("will read qarray from %s\n",qarrayfilename);
    qarrayfile=fopen(qarrayfilename,"r");
    fscanf(qarrayfile,"%d",&nqmax);
    qarray=new double[nqmax];
    for(iq=0;iq<nqmax;iq++){
      fscanf(qarrayfile,"%lf",&qarray[iq]);
    }
    fclose(qarrayfile);
  }
  else{
    qarray=new double[nqmax];
    for(iq=0;iq<nqmax;iq++){
      qarray[iq]=(iq+0.5)*delq;
    }
  }
}

void Twavefunction::InitArrays(){
  int iq,l,ichannel;
  
  eta=new double[nqmax];
  planewave=new Tplanewave*[nqmax];

  partwave=new Tpartwave **[ellmax+1];
  for(l=0;l<=ellmax;l++){
    partwave[l]=new Tpartwave *[nqmax];
    for(iq=0;iq<nqmax;iq++) partwave[l][iq]=NULL;
  }

  delta=new double *[nchannels];
  ddeltadq=new double *[nchannels];
  Wepsilon=new double *[nchannels];
  for(ichannel=0;ichannel<nchannels;ichannel++){
    delta[ichannel]=new double[nqmax];
    ddeltadq[ichannel]=new double[nqmax];
    Wepsilon[ichannel]=new double[nqmax];
  }

  ell=new int[nchannels];
  channelweight=new double[nchannels];

}

double Twavefunction::getpsisquared(double q,double r,double ctheta){
  int iq,iqlow,iqhigh;
  double wlow,whigh,interpolate,qscaled,rscaled;

  if(generic==1){
    qscaled=q*(muscale/mu)*q1q2scale/double(q1q2);
    rscaled=q*r/qscaled;
  }
  else{
    qscaled=q;
    rscaled=r;
  }

  if(delq<0){
    iq=0;
    while(qscaled>qarray[iq]+1.0E-5 && iq<nqmax){
      iq+=1;
    }
    iqlow=iq-1;
    iqhigh=iq;
  }
  else{
    iqlow=int(floor(qscaled/delq)-0.5);
    iqhigh=iqlow+1;
  }
  if(iqhigh==nqmax && nqmax>1 
     && (qscaled-qarray[nqmax-1])<0.5*(qarray[nqmax-1]-qarray[nqmax-2])){
    iqhigh=nqmax-1;
    iqlow=iqhigh-1;
  }
    
  if(iqhigh==0) return calcpsisquared(0,rscaled,ctheta);
  else if(iqhigh>=nqmax && qscaled-qarray[nqmax-1]>1.0E-5) return 1.0;
  else{
    wlow=(qarray[iqhigh]-qscaled)/(qarray[iqhigh]-qarray[iqlow]);
    whigh=1.0-wlow;
    if(fabs(wlow)<1.0E-5) interpolate=calcpsisquared(iqhigh,rscaled,ctheta);
    else if(fabs(whigh)<1.0E-5) 
      interpolate=calcpsisquared(iqlow,rscaled,ctheta);
    else{
      interpolate=wlow*calcpsisquared(iqlow,rscaled,ctheta)
	+whigh*calcpsisquared(iqhigh,rscaled,ctheta);
    }
    return interpolate;
  }
}

double Twavefunction::getpsisquared(double *pa,double *xa,
				    double *pb,double *xb){
  double q,r,ctheta;
  getqrctheta(pa,xa,pb,xb,&q,&r,&ctheta);
  return getpsisquared(q,r,ctheta);
}

double Twavefunction::calcpsisquared(int iq,double r,double ctheta){
  return 1.0;  // This is a dummy function to be overwritten by inherited class
}

void Twavefunction::InitWaves(){
  int iq,ichannel;
  double q,e1,e2;
  
  for(iq=0;iq<nqmax;iq++){
    q=qarray[iq];
    e1=sqrt(m1*m1+q*q);
    e2=sqrt(m2*m2+q*q);
    eta[iq]=double(q1q2)*e1*e2/((e1+e2)*137.036*q);
    planewave[iq]=new Tplanewave(eta[iq],q1q2,q);
  }
  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      if(partwave[ell[ichannel]][iq]==NULL){
	partwave[ell[ichannel]][iq]
	  =new Tpartwave(eta[iq],q1q2,q,ell[ichannel]);
      }
    }
  }
}

void Twavefunction::getqrctheta(double *p1,double *r1,double *p2,double *r2,
				double *q,double *r,double *ctheta){
  int alpha;
  const double g[4]={1.0,-1.0,-1.0,-1.0};
  double n[4],qvec[4],rvec[4],nnorm,ndotq,ndotr,qdotr;

  nnorm=0.0;
  ndotq=0.0;
  ndotr=0.0;
  for(alpha=0;alpha<4;alpha++){
    n[alpha]=p1[alpha]+p2[alpha];
    qvec[alpha]=0.5*(p1[alpha]-p2[alpha]);
    rvec[alpha]=r1[alpha]-r2[alpha];
    nnorm+=g[alpha]*n[alpha]*n[alpha];
    ndotq+=n[alpha]*qvec[alpha]*g[alpha];
    ndotr+=n[alpha]*rvec[alpha]*g[alpha];
  }
  nnorm=sqrt(nnorm);
  ndotq=ndotq/nnorm;
  ndotr=ndotr/nnorm;

  *ctheta=0.0;
  *r=0.0;
  *q=0.0;
  for(alpha=0;alpha<4;alpha++){
    n[alpha]=n[alpha]/nnorm;
    rvec[alpha]=rvec[alpha]-ndotr*n[alpha];
    qvec[alpha]=qvec[alpha]-ndotq*n[alpha];
    *r-=g[alpha]*rvec[alpha]*rvec[alpha];
    *q-=g[alpha]*qvec[alpha]*qvec[alpha];
    *ctheta+=g[alpha]*rvec[alpha]*qvec[alpha];
  }
  if(*r<0.0 || *q<0.0 || fabs(*ctheta)>sqrt(*q**r)){
    printf("Disaster, r^2=%g, q^2=%g, ctheta=%g\n",*r,*q,*ctheta/sqrt(*r**q));
    exit(1);
  }
  *r=sqrt(*r);
  *q=sqrt(*q);
  *ctheta=*ctheta/(*r**q);
}

void Twavefunction::printCdelta(double Rx,double Ry,double Rz){
  double q,clocal;
  int ichannel,iq;
  printf("! Qinv  C(Q)_estimated ~ ddelta/dq\n");
  for(iq=0;iq<nqmax;iq++){
    q=qarray[iq];
    clocal=1.0;
    for(ichannel=0;ichannel<nchannels;ichannel++){
      clocal+=channelweight[ichannel]*(2.0*ell[ichannel]+1.0)
	*(2.0*PI)*pow(HBARC,3)
	/(q*q*Rx*Ry*Rz*pow(4.0*PI,1.5))
	*ddeltadq[ichannel][iq];
    }
    printf("%6.2f  %8.4f  %g\n",q,clocal,4.0*q*q*(clocal-1.0));    
  }
}
