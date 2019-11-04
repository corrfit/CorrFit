// See wavefunction.cc for parent class

#ifndef PIPI_PHASESHIFTS_INCLUDED
#include "pipi_phaseshifts/pipi_phaseshifts.cc"
#endif

Twavefunction_pipluspiplus::Twavefunction_pipluspiplus(char *parsfilename){
  int iq,ichannel;
  double q;
  int *I;

  ParsInit(parsfilename);

  m1=MPI;
  m2=MPI;;
  q1q2=1;
  nchannels=2;

  InitArrays();
  ell[0]=0;
  ell[1]=2;
  InitWaves();
  channelweight[0]=channelweight[1]=2.0;

  I=new int[nchannels];
  I[0]=I[1]=2;
  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      getphaseshift_pipi(I[ichannel],ell[ichannel],q,&delta[ichannel][iq],
			   &ddeltadq[ichannel][iq]);
      if(q1q2!=0)
	phaseshift_CoulombCorrect(ell[ichannel],q,eta[iq],
				  delta[ichannel][iq],ddeltadq[ichannel][iq]);
      Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
	-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
	+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
      Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
	/(4.0*PI*pow(epsilon,3));
    }
  }
}

double Twavefunction_pipluspiplus::calcpsisquared(int iq,double r,
						 double ctheta){
  double psisquared,x,dpsi2,root2=sqrt(2.0);
  double delta_s,delta_p,delta_d;
  complex<double> psi,hstar,psi0,psisymm,psia,psib;
  double q;
  int ichannel;

  if(iq>=nqmax){
    psisquared=1.0;
  }
  else{
    q=qarray[iq];
    psia=planewave[iq]->planewave(r,ctheta);
    psib=planewave[iq]->planewave(r,-ctheta);
    psi=(psia+psib)/root2;
    if(r<epsilon){
      psisquared=real(psi*conj(psi));
      if(STRONG==1){
	for(ichannel=0;ichannel<nchannels;ichannel++){
	  dpsi2=(2.0*ell[ichannel]+1.0)*2.0*PI*Wepsilon[ichannel][iq]
	    *pow(HBARC,3)/(q*q);
	  psisquared+=2.0*dpsi2; // factor of two due to symmetrization
	}
      }
    }
    else{
      if(STRONG==1){
	x=q*r/HBARC;
	delta_s=delta[0][iq];
	delta_d=delta[1][iq];
	//delta_s=delta_d=0.0;
	
	ichannel=0;  // This will be the s wave
	hstar=partwave[ell[ichannel]][iq]->GetPhiIncoming(r)/x;
	psi+=root2*0.5*hstar*(exp(-2.0*ci*delta_s)-1.0);

	ichannel=1; // This will be the d wave
	hstar=partwave[ell[ichannel]][iq]->GetPhiIncoming(r)/x;
	psi+=root2*0.5*ci*ci*5.0*legendre(2,ctheta)
	  *hstar*(exp(-2.0*ci*delta_d)-1.0);
      }
      psisquared=real(psi*conj(psi));
    }
  }

  //printf("psisquared=%g\n",psisquared);
  return psisquared;

}
