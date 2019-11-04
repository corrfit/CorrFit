void Twavefunction_pk::read_phaseshifts(){
  int iq,iread,iqsmooth;
  double qsmooth,q,a,elab,plab,roots,deltaread[61],qread[61];
  char dummy[200];
  FILE *fptr;


  iqsmooth=1;
  qsmooth=25.0;
  fptr=fopen("wf/pk_phaseshifts/s11.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MPI*MPI+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(triangle(roots,MPROTON,MPI));
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth) a=tan(deltaread[iread])/qread[iread];
  }
  fclose(fptr);
  
  iread=0;
  for(iq=0;iq<nqmax;iq++){
    q=qarray[iq];
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[0][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[0][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[0][iq]=atan(fabs(a*q))*a/fabs(a);
      ddeltadq[0][iq]=a*pow(cos(delta[0][iq]),2);
    }
  }

  // ________________________________________________________

  iqsmooth=5;
  qsmooth=75.0;
  fptr=fopen("wf/pk_phaseshifts/p11.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MPI*MPI+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(triangle(roots,MPROTON,MPI));
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth){
      a=tan(deltaread[iread])/pow(qread[iread],3);
    }
  }
  fclose(fptr);

  iread=0;
  for(iq=0;iq<nqmax;iq++){
    q=qarray[iq];
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[1][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[1][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[1][iq]=atan(fabs(a*q*q*q))*a/fabs(a);
      ddeltadq[1][iq]=3*a*q*q*pow(cos(delta[1][iq]),2);
    }
  }

  // ________________________________________________________

  iqsmooth=5;
  qsmooth=75.0;
  fptr=fopen("wf/pk_phaseshifts/p13.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MPI*MPI+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(triangle(roots,MPROTON,MPI));
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth) a=tan(deltaread[iread])/pow(qread[iread],3);
  }
  fclose(fptr);

  iread=0;
  for(iq=0;iq<nqmax;iq++){
    q=qarray[iq];
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[2][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[2][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[2][iq]=atan(fabs(a*q*q*q))*a/fabs(a);
      ddeltadq[2][iq]=3*a*q*q*pow(cos(delta[2][iq]),2);
    }
  }

}

