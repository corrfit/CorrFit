Twavefunction_pk::Twavefunction_pk(char *parsfilename) : Twavefunction(){
  int iq,ichannel;
  double q;

  ParsInit(parsfilename);

  m1=MPROTON;
  m2=MKAON;
  q1q2=1;
  if(COULOMB==0) q1q2=0;
  nchannels=3;

  InitArrays();

  ell[0]=0;
  ell[1]=1;
  ell[2]=1;

  InitWaves();

  channelweight[0]=1.0;
  channelweight[1]=1.0/3.0;
  channelweight[2]=2.0/3.0;

  read_phaseshifts();

  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      if(q1q2!=0)
      phaseshift_CoulombCorrect(ell[ichannel],q,eta[iq],
      				delta[ichannel][iq],ddeltadq[ichannel][iq]);
      //printf("ichannel=%d, q=%g, delta=%g\n",
      //     ichannel,q,delta[ichannel][iq]*180.0/PI);
      Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
	-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
	+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
      Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
	/(4.0*PI*pow(epsilon,3));
    }
  }
}

double Twavefunction_pk::calcpsisquared(int iq,double r,double ctheta){
  double q,psisquared,x,dpsi2;
  double delta_s11,delta_p11,delta_p13;
  complex<double> psi,hstar0,hstar1,psi0;
  int ipartial,ichannel;

  q=qarray[iq];
  x=q*r/HBARC;
  hstar0=partwave[0][iq]->GetPhiIncoming(r)/x;
  hstar1=partwave[1][iq]->GetPhiIncoming(r)/x;
  psi0=planewave[iq]->planewave(r,ctheta);
  if(iq>=nqmax){
    printf("iq too large!\n");
    exit(1);
  }

  if(r<epsilon){
    psisquared=real(psi0*conj(psi0));
    if(STRONG==1){
      for(ichannel=0;ichannel<nchannels;ichannel++){
	dpsi2=(2.0*ell[ichannel]+1.0)*2.0*PI*Wepsilon[ichannel][iq]
	  *pow(HBARC,3)/(q*q);
	psisquared+=channelweight[ichannel]*dpsi2;
      }
    }
  }
  else{
    if(STRONG==1){
      delta_s11=delta[0][iq];
      delta_p11=delta[1][iq];
      delta_p13=delta[2][iq];
      
      // s-wave m_s=+1/2
      psi=psi0+0.5*hstar0*(exp(-2.0*ci*delta_s11)-1.0);
      
      // m_s still equals 1/2, but now ell=1
      psi+=0.5*hstar1*((2.0/3.0)*exp(-2.0*ci*delta_p13)
		       +(1.0/3.0)*exp(-2.0*ci*delta_p11)-1.0)*ci*(3.0)*ctheta;
      
      psisquared=real(psi*conj(psi));
      
      //Now consider m_s=-1/2, m_ell=1, Note: doesn't interfere with m_s=+1/2
      psi=0.5*hstar1*(exp(-2.0*ci*delta_p13)-exp(-2.0*ci*delta_p11))
	*ci*(3.0)*sqrt(1.0-ctheta*ctheta)/3.0;
    }
    else{
      psi=psi0;
    }
    psisquared+=real(psi*conj(psi));
  }
  return psisquared;

}

#include "pk_phaseshifts/pk_phaseshift_reader.cc"
