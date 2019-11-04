// See wavefunction.cc for parent class

Twavefunction_kpi::Twavefunction_kpi(char *parsfilename){
  int iq,ichannel;
  double q;

  ParsInit(parsfilename);

  m1=MKAON;
  m2=MPI;
  q1q2=1;
  nchannels=0;
  InitArrays();
  InitWaves();

}

double Twavefunction_kpi::calcpsisquared(int iq,double r,double ctheta){
  double psisquared,q;
  complex<double> psi,psi0;

  if(iq>=nqmax) psisquared=1.0;
  else{
    q=qarray[iq];
    psi0=planewave[iq]->planewave(r,ctheta);
    if(iq>=nqmax){
      printf("iq too large!\n");
      exit(1);
    }
    
    psisquared=real(psi0*conj(psi0));
  }
  return psisquared;

}


