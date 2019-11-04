#include "StHbtCorrFctnDirectYlmPurity.h"
#include <TMath.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <iostream>

using namespace std;
#define MAXJM 25

#ifdef __ROOT__ 
ClassImp(StHbtCorrFctnDirectYlmPurity)
#endif

  StHbtCorrFctnDirectYlmPurity::StHbtCorrFctnDirectYlmPurity(const char *name, int maxl, int ibin=30, double vmin=0.0, double vmax=0.3, int p1Type=1, int p2Type=2):
  StHbtCorrFctnDirectYlm(name, maxl, ibin, vmin, vmax)
{
  mP1Type = p1Type;
  mP2Type = p2Type;
}


StHbtCorrFctnDirectYlmPurity::StHbtCorrFctnDirectYlmPurity()
{
  StHbtCorrFctnDirectYlm::StHbtCorrFctnDirectYlm("StHbtCorrFctnDirectYlmPurity",2, 30, 0.0, 0.3);
  mP1Type = 1;
  mP2Type = 2;
}

StHbtCorrFctnDirectYlmPurity::~StHbtCorrFctnDirectYlmPurity()
{
  this->StHbtCorrFctnDirectYlm::~StHbtCorrFctnDirectYlm();
}

void StHbtCorrFctnDirectYlmPurity::AddRealPair(const StHbtPair *aPair)
{
  double pPurity = 0.0;

  switch (mP1Type)
    {
    case 1:
      pPurity =  aPair->track1()->GetPionPurity();
      break;
    case 2:
      pPurity =  aPair->track1()->GetKaonPurity();
      break;
    case 3:
      pPurity =  aPair->track1()->GetProtonPurity();
      break;
    case 4:
      pPurity =  aPair->track1()->Track()->PidProbElectron();
      break;
    }

  //  cout << "First particle purity: " << pPurity << endl;
  
  switch (mP2Type)
    {
    case 1:
      pPurity *=  aPair->track2()->GetPionPurity();
      break;
    case 2:
      pPurity *=  aPair->track2()->GetKaonPurity();
      break;
    case 3:
      pPurity *=  aPair->track2()->GetProtonPurity();
      break;
    case 4:
      pPurity *=  aPair->track2()->Track()->PidProbElectron();
      break;
    }

  StHbtCorrFctnDirectYlm::AddRealPair(aPair->dKOut(), aPair->dKSide(), aPair->dKLong(), pPurity);
  StHbtCorrFctnDirectYlm::AddMixedPair(aPair->dKOut(), aPair->dKSide(), aPair->dKLong(), 1.0);  
}

void StHbtCorrFctnDirectYlmPurity::AddMixedPair(const StHbtPair *aPair)
{
}

void StHbtCorrFctnDirectYlmPurity::Finish()
{
  StHbtCorrFctnDirectYlm::Finish();
}

StHbtString StHbtCorrFctnDirectYlmPurity::Report()
{
  return "StHbtCorrFctnDirectYlmPurity::Finish";
}

