// Parent class is Twavefunction

class Twavefunction_pk : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  void read_phaseshifts();
  Twavefunction_pk(char *parsfilename);
};
// Parent class is Twavefunction

class Twavefunction_ppi : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  void read_phaseshifts();
  Twavefunction_ppi(char *parsfilename);
};

class Twavefunction_generic : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  void reset(int q1q2,double m1,double m2,double symmweight);
  Twavefunction_generic(char *parsfilename,int q1q2,double m1,double m2,double symmweight);
  
};

class Twavefunction_Xipi : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  void get_phaseshifts_Xistar(double q,double &delta,double &ddeltadq);
  Twavefunction_Xipi(char *parsfilename);
};

class Twavefunction_pipluspiminus : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  Twavefunction_pipluspiminus(char *parsfilename);
};

class Twavefunction_pipluspiplus : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  Twavefunction_pipluspiplus(char *parsfilename);
};

class Twavefunction_lambdalambda : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  Twavefunction_lambdalambda(char *parsfilename);
  void GetPhaseshifts();
};

class Twavefunction_lambdalambda_antiparspin : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  Twavefunction_lambdalambda_antiparspin(char *parsfilename);
  void GetPhaseshifts();
};

class Twavefunction_lambdalambda_parspin : public Twavefunction {
 public:
  double calcpsisquared(int iq,double r,double ctheta);
  Twavefunction_lambdalambda_parspin(char *parsfilename);
  void GetPhaseshifts();
};

class Twavefunction_kpi : public Twavefunction {
public:
  double calcpsisquared(int iq,double r,double ctheta);
  Twavefunction_kpi(char *parsfilename);
};

