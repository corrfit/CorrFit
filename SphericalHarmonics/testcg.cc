#include "mmatrix.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
  double jot1 = atof(argv[1]);
  double em1  = atof(argv[2]);
  double jot2 = atof(argv[3]);
  double em2  = atof(argv[4]);
  double jot3 = atof(argv[5]);
  double em3  = atof(argv[6]);

  double cgc  = MMatrix::ClebschGordan(jot1, em1, jot2, em2, jot3, em3);
  double wgs  = MMatrix::WignerSymbol(jot1, em1, jot2, em2, jot3, em3);

  cout << "j1 m1 j2 m2 j3 m3 cgc wgs " 
       << jot1 << " " << em1 << " "
       << jot2 << " " << em2 << " "
       << jot3 << " " << em3 << " "
       << cgc  << " " << wgs << endl;

}
