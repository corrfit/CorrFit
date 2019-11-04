// See wavefunction.cc for parent class

Twavefunction_Xipi::Twavefunction_Xipi(char *parsfilename){
  int iq,ichannel;
  double q;

  ParsInit(parsfilename);

  m1=1321.3;
  m2=139.58;
  q1q2=-1;
  if(COULOMB==0) q1q2=0;
  nchannels=1;

  InitArrays();

  ell[0]=1;

  InitWaves();

  channelweight[0]=1.0;

  ichannel=0;
  for(iq=0;iq<nqmax;iq++){
    q=qarray[iq];
    get_phaseshifts_Xistar(q,delta[0][iq],ddeltadq[0][iq]);
    if (q1q2!=0)
      phaseshift_CoulombCorrect(ell[ichannel],q,eta[iq],
				delta[ichannel][iq],ddeltadq[ichannel][iq]);
    printf("ichannel=%d, q=%g, E=%g, delta=%g, ddelta/dq=%g\n",
	   ichannel,q,sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2),
	   delta[ichannel][iq]*180.0/PI,ddeltadq[ichannel][iq]*180.0/PI);
    Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
      -GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
      +GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
    Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
      /(4.0*PI*pow(epsilon,3));
  }
}

double Twavefunction_Xipi::calcpsisquared(int iq,double r,double ctheta){
  double psisquared,x,dpsi2,q;
  complex<double> psi,hstar1,psi0,delpsi;
  int ipartial,ichannel;

  if(iq>=nqmax){
    psisquared=1.0;
  }
  else{
    q=qarray[iq];
    x=q*r/HBARC;
    psi0=planewave[iq]->planewave(r,ctheta);

    if(r<epsilon){
      psisquared=real(psi0*conj(psi0));
      for(ichannel=0;ichannel<nchannels;ichannel++){
	dpsi2=(2.0*ell[ichannel]+1.0)*2.0*PI*Wepsilon[ichannel][iq]
	  *pow(HBARC,3)/(q*q);
	psisquared+=channelweight[ichannel]*dpsi2;
      }
    }
    else{
      hstar1=partwave[1][iq]->GetPhiIncoming(r)/x;
      delpsi=0.5*hstar1*(exp(-2.0*ci*delta[0][iq])-1.0)
	*ci*(3.0)*ctheta;
      // m_s equals 1/2, ell=1
      psi=psi0+(4.0/9.0)*delpsi;
      psisquared=real(psi*conj(psi));

      //This accounts for the other three channels, where the incoming
      // wave was either in the Xi0-pi0 channel or if m_ell=1, m_s=-1/2.
      psisquared+=(19.0/81.0)*real(delpsi*conj(delpsi));
      psisquared=real(psi0*conj(psi0));

    }
  }
  return psisquared;

}

void Twavefunction_Xipi::get_phaseshifts_Xistar(double q,double &delta,double &ddeltadq){
  double answer,tandelta,M,Epi,EXi,Lambda;
  const double MXi=1321.3,Lambda0=9.1,MXistar=1530.0;
  const double q0=sqrt(triangle(MXistar,MXi,MPI));
  EXi=sqrt(MXi*MXi+q*q);
  Epi=sqrt(MPI*MPI+q*q);
  Lambda=Lambda0*pow(q/q0,3);
  tandelta=0.5*Lambda/(MXistar-EXi-Epi);
  delta=atan(tandelta);
  if(delta<0.0) delta+=PI;
  ddeltadq=pow(cos(delta),2)*tandelta*((3.0/q)+(q/(MXistar-EXi-Epi))
				       *((1.0/EXi)+(1.0/Epi)));
}
