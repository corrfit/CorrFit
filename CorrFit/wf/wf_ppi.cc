Twavefunction_ppi::Twavefunction_ppi(char *parsfilename) : Twavefunction() {
  int iq,ichannel;
  double q;
  ParsInit(parsfilename);

  m1=MPI;
  m2=MPROTON;
  q1q2=1;
  if(COULOMB==0) q1q2=0;
  nchannels=3;

  InitArrays();
  printf("Arrays Initialized\n");

  ell[0]=0;
  ell[1]=1;
  ell[2]=1;

  InitWaves();
  printf("Partial Waves Initialized\n");

  channelweight[0]=1.0;
  channelweight[1]=1.0/3.0;
  channelweight[2]=2.0/3.0;

  read_phaseshifts();

  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      if(COULOMB==1)
	phaseshift_CoulombCorrect(ell[ichannel],q,eta[iq],
				  delta[ichannel][iq],ddeltadq[ichannel][iq]);
      printf("ichannel=%d, q=%g, delta=%g\n",
	     ichannel,q,delta[ichannel][iq]*180.0/PI);
      Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
	-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
	+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
      Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
	/(4.0*PI*pow(epsilon,3));
    }
  }
  printf("Initialization finished\n");
}

double Twavefunction_ppi::calcpsisquared(int iq,double r,double ctheta){
  double psisquared,x,dpsi2,q;
  double delta_s31,delta_p31,delta_p33;
  complex<double> psi,hstar,psi0;
  int ichannel;

  q=qarray[iq];
  if(iq>=nqmax){
    printf("iq too large!\n");
    exit(1);
  }
  psi0=planewave[iq]->planewave(r,ctheta);

  if(STRONG==1){
    if(r<epsilon){
      psisquared=real(psi0*conj(psi0));
      for(ichannel=0;ichannel<nchannels;ichannel++){
	dpsi2=(2.0*ell[ichannel]+1.0)*2.0*PI*Wepsilon[ichannel][iq]
	  *pow(HBARC,3)/(q*q);
	psisquared+=channelweight[ichannel]*dpsi2;
      }
    }
    else{
      x=q*r/HBARC;
      delta_s31=delta[0][iq];
      delta_p31=delta[1][iq];
      delta_p33=delta[2][iq];
 
      psi=psi0;
      // this refers to m_s=+1/2 s wave
      ichannel=0;
      hstar=partwave[ell[ichannel]][iq]->GetPhiIncoming(r)/x;
      psi+=0.5*hstar*(exp(-2.0*ci*delta_s31)-1.0);
       
      // m_s still equals 1/2, but now ell=1
      ichannel=1;
      hstar=partwave[ell[ichannel]][iq]->GetPhiIncoming(r)/x;
      psi+=0.5*hstar*((2.0/3.0)*exp(-2.0*ci*delta_p33)
      		   +(1.0/3.0)*exp(-2.0*ci*delta_p31)-1.0)*ci*(3.0)*ctheta; 
      psisquared=real(psi*conj(psi));
       
      //Now consider m_s=-1/2, m_ell=1, Note: doesn't interfere with m_s=+1/2
      ichannel=2;
      psi=0.5*hstar*(exp(-2.0*ci*delta_p33)-exp(-2.0*ci*delta_p31))
	*ci*(3.0)*sqrt(1.0-ctheta*ctheta)/3.0;
      psisquared+=real(psi*conj(psi));

    }
  }
  else psisquared=real(psi0*conj(psi0));
  return psisquared;

}

#include "ppi_phaseshifts/ppi_phaseshift_reader.cc"
