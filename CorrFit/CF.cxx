#include "CF.h"

CF::CF()
{
  mContent = 0;
  mErr2 = 0;
}

CF::~CF()
{
}

double* CF::GetContent()
{
  return mContent;
}

double* CF::GetError2()
{
  return mErr2;
}

const double* CF::GetContent() const
{
  return mContent;
}

const double* CF::GetError2() const
{
  return mErr2;
}

