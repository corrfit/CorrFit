#include "ExpCF.h"

ExpCF::ExpCF()
{
  mNormN = 0;
  mNormP = 0;
  mPurCorr = 0.0;
  mMomRes = 0.0;
  mPurity = 0;
}

ExpCF::~ExpCF()
{
}

int ExpCF::GetNFitBin()
{
  return mNFitBin;
}

double ExpCF::GetMomResScale()
{
  return mMomRes;
}

double *ExpCF::GetContent()
{
  return CF::GetContent();
}

double *ExpCF::GetPurity()
{
  return mPurity;
}
