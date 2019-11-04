bool comparestrings(char *s1,char *s2){
  int length1,length2,ic;
  bool answer;
  answer=0;
  length1=strlen(s1);
  length2=strlen(s2);
  if(length1==length2){
    answer=1;
    for(ic=0;ic<length1;ic++){
      if(s1[ic]!=s2[ic]){
	answer=0;
	goto NOMATCH;
      }
    }
  }
 NOMATCH:
  return answer;
}

double legendre(int ell,double ctheta){
  double answer,a0,a1;
  int i;
  if(ell==0) answer=1.0;
  if(ell==1) answer=ctheta;
  if(ell>1){
    a0=1.0;
    a1=ctheta;
    for(i=1;i<ell;i++){
      answer=(1.0/double(i+1))*((2.0*i+1)*ctheta*a1-double(i)*a0);
      a0=a1;
      a1=answer;
    }
  }
  return answer;
}
 
double dgamma(int mm){
  /* This calc.s gamma functions which are in the form gamma(n)
     where n is an int > 0. */
  double dg;
  int j;
  dg=1.0;
  if(mm<1) {
    for(j=1;j<=-mm+1;j++){
      dg=dg/(1.0+double(-j));
    }
  }
  if(mm>1){
    for(j=1;j<=mm-1;j++){
      dg=dg*double(j);
    }
  }
  return dg;
}

double triangle(double m0,double m1,double m2){
  double answer,m0sq,m1sq,m2sq;
  if(m0<m1+m2) {
    printf("Disaster with triangle\n");
    exit(1);
  }
  m0sq=m0*m0;m1sq=m1*m1;m2sq=m2*m2;
  answer=m0sq*m0sq+m1sq*m1sq+m2sq*m2sq;
  answer=answer-2.0*(m0sq*m1sq+m0sq*m2sq+m1sq*m2sq);
  answer=answer/(4.0*m0sq);
  return answer;
}
