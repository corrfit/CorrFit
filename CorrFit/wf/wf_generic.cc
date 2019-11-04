// See wavefunction.cc for parent class

Twavefunction_generic::Twavefunction_generic(
					     char *parsfilename,
					     int q1q2set,double m1set,
					     double m2set,
					     double symmweightset): Twavefunction(){
  generic=1;
  ParsInit(parsfilename);
  m1=m1set;
  m2=m2set;
  muscale=m1*m2/(m1+m2);
  mu=muscale;
  symmweight=symmweightset;
  q1q2scale=q1q2set;
  q1q2=q1q2scale;
  nchannels=0;
  InitArrays();
  InitWaves();
  printf("initialization finished\n");
}

void Twavefunction_generic::reset(int q1q2set,double m1set,double m2set,
				  double symmweightset){
  m1=m1set;
  m2=m2set;
  mu=m1*m2/(m1+m2);
  if(q1q2*q1q2set<0){
    printf("Illegal: Trying to reset q1q2 to opposite charge\n");
    exit(1);
  }
  q1q2=q1q2set;
  symmweight=symmweightset;
}

double Twavefunction_generic::calcpsisquared(int iq,double r,
					     double ctheta){
  double psisquared,e1,e2,qscaled,rscaled,vscaled,asymmweight,q;
  complex<double> psi1,psi2,psisymm,psiasymm;

  if(iq>=nqmax){
    printf("iq too large!\n");
    psisquared=1.0;
  }
  else{
    psi1=planewave[iq]->planewave(r,ctheta);
    psi2=conj(psi1);
    
    asymmweight=1.0-symmweight;
    psisymm=(psi1+psi2)/ROOT2;
    psiasymm=(psi1-psi2)/ROOT2;
    
    psisquared=symmweight*real(psisymm*conj(psisymm))
      +asymmweight*real(psiasymm*conj(psiasymm));
  }

  return psisquared;

}




