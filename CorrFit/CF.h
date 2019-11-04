#ifndef _CORRFIT_CF_H_
#define _CORRFIT_CF_H_

class CF 
{
 public:
  CF();
  virtual ~CF();
  virtual void Write() = 0;
  //  virtual void Read(TFile *aFile, TKey *aKey) = 0;
  virtual void Normalize() = 0;
  
  virtual double* GetContent();
  virtual double* GetError2();
  virtual const double* GetContent() const;
  virtual const double* GetError2() const;

 protected:
  double *mContent;
  double *mErr2;
};


#endif
