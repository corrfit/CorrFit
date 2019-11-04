// See wavefunction.cc for parent class

#ifndef PIPI_PHASESHIFTS_INCLUDED
#include "pipi_phaseshifts/pipi_phaseshifts.cc"
#endif

class Twavefunction_pipluspiminus : public Twavefunction {
public:
  double getpsisquared(double q,double r,double ctheta);
  void init(char *parsfilename);
};

void Twavefunction_pipluspiminus::init(char *parsfilename){
  int iq,ichannel,*I;
  double q;

  ParsInit(parsfilename);

  m1=MKAON;
  m2=MKAON;
  q1q2=-1;
  nchannels=5;

  InitArrays();
  ell[0]=0;
  ell[1]=2;
  ell[2]=1;
  ell[3]=0;
  ell[4]=2;
  InitWaves();
  channelweight[0]=channelweight[1]=2.0/3.0;
  channelweight[2]=1.0;
  channelweight[3]=channelweight[4]=1.0/3.0;
  I=new int[nchannels];
  I[0]=I[1]=0;
  I[2]=1;
  I[3]=I[4]=2;

  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=delq*(0.5+iq);
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

double Twavefunction_pipluspiminus::getpsisquared(double q,double r,
						  double ctheta){
  double psisquared,x,dpsi2,root2=sqrt(2.0);
  double delta_s,delta_p,delta_d;
  complex<double> psi,hstar0,hstar1,hstar2,psi0,psisymm,psia,psib;
  int iq,ipartial,ichannel;

  iq=int(floor(q/delq));
  if(iq>=nqmax){
    psisquared=1.0;
  }
  else{
    q=(0.5+iq)*delq;
    x=q*r/HBARC;
    hstar0=partwave[0][iq]->GetPhiIncoming(r)/x;
    hstar1=partwave[1][iq]->GetPhiIncoming(r)/x;
    hstar2=partwave[2][iq]->GetPhiIncoming(r)/x;



    if(r<epsilon){
      psi0=conj(planewave[iq]->planewave(r,-ctheta));
      psisquared=real(psi0*conj(psi0));
      for(ichannel=0;ichannel<nchannels;ichannel++){
	dpsi2=(2.0*ell[ichannel]+1.0)*2.0*PI*Wepsilon[ichannel][iq]
	  *pow(HBARC,3)/(q*q);
	//dpsi2=0.0;
	if(ichannel==0 || ichannel==1)
	  psisquared+=2.0*dpsi2/3.0;
	if(ichannel==2)
	  psisquared+=dpsi2;
	if(ichannel==3 || ichannel==4)
	  psisquared+=(1.0/3.0)*dpsi2;
      }
    }
    else{
      psia=planewave[iq]->planewave(r,ctheta);
      psib=planewave[iq]->planewave(r,-ctheta);

      // For I=0
      psi=(psia+psib)/root2;
      delta_s=delta[0][iq];
      delta_d=delta[1][iq];
      ichannel=0;  // This will be the s wave
      psi+=root2*0.5*hstar0*(exp(-2.0*ci*delta_s)-1.0);
      ichannel=1; // This will be the d wave
      psi+=root2*0.5*ci*ci*5.0*legendre(2,ctheta)*hstar2*(exp(-2.0*ci*delta_d)-1.0);
      psisquared=(1.0/3.0)*real(psi*conj(psi));

      // For I=1
      psi=(psia-psib)/root2;
      delta_p=delta[2][iq];
      ichannel=2;  // This will be the p wave
      psi+=root2*0.5*ci*3.0*legendre(1,ctheta)*hstar1*(exp(-2.0*ci*delta_p)-1.0);
      psisquared+=0.5*real(psi*conj(psi));

      // For I=2
      psi=(psia+psib)/root2;
      delta_s=delta[3][iq];
      delta_d=delta[4][iq];
      ichannel=3;  // This will be the s wave
      psi+=root2*0.5*hstar0*(exp(-2.0*ci*delta_s)-1.0);
      ichannel=4; // This will be the d wave
      psi+=root2*0.5*ci*ci*5.0*legendre(2,ctheta)
	*hstar2*(exp(-2.0*ci*delta_d)-1.0);
      psisquared+=(1.0/6.0)*real(psi*conj(psi));
    }
  }
  return psisquared;

}
