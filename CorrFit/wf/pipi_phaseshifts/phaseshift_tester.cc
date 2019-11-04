#include <cstdlib>
#include <cmath>
#include <cstdio>

const double PI=3.14159265358979323844;
const double MPI=139.58;
const double HBARC=197.3269602;

#include "pipi_phaseshifts.cc"

using namespace std;

int main(){
  int iq,I,ell;
  double q,e,delq=50.0,delta,ddeltadq,olddelta;
  printf("Enter I and ell: ");
  scanf("%d %d",&I,&ell);
  olddelta=0.0;
  for(q=0.5*delq;q<2000;q+=delq){
    e=2*sqrt(q*q+MPI*MPI);
    getphaseshift_pipi(I,ell,q,&delta,&ddeltadq);
    delta*=(180.0/PI);
    ddeltadq*=(180.0/PI);
    printf("%1d %1d %6.1f %6.2f %6.3f %6.3f\n",
	   I,ell,e,delta,ddeltadq,(delta-olddelta)/delq);
    olddelta=delta;
  }

  return 0;
}


