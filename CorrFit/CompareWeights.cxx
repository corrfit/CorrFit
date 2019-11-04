#include "PairManager.h"
#include "StandAloneFsiKisiel.h"
#include "StandAloneFsiLednicky.h"
#include "SourceModelGausRInvJacobian.h"
#include "SourceModelExpGausRInv.h"
#include "TH1D.h"

int main(int argc, char **argv)
{
  int pipi = 0;
  int pmax = 1000000;

  PairManager *tPairManager;
  if (pipi)
    tPairManager = new PairManager("pairpipi.root");
  else
    tPairManager = new PairManager("pairpp.root");

  Pair *tPair;
  double tParameters[3] = {6.0,0.0,0.12};

  SourceModelExpGausRInv *tSourceModel = new SourceModelExpGausRInv();
  tSourceModel->SetModelParameters(tParameters);
  Pair::SetSourceModel(tSourceModel);

  if (pmax > tPairManager->GetPairCount())
    pmax = tPairManager->GetPairCount();

  for (int iter=0; iter<pmax; iter++) {
    tPair = tPairManager->ReadPair();
    tPairManager->StoreFitPair(tPair);
    tPair = tPairManager->GetFitPair(iter);
    tPair->InitRandVar();
  }

  TH1D *wsumled = new TH1D("wsumled","wsumled",50, 0.0, 0.2);
  TH1D *wsumkis = new TH1D("wsumkis","wsumkis",50, 0.0, 0.2);
  TH1D *wsumden = new TH1D("wsumden","wsumden",50, 0.0, 0.2);
  TH1D *dndeta  = new TH1D("dNdeta ","dEdeta ",50, 0.0, 0.2);
  
  wsumled->Sumw2();
  wsumkis->Sumw2();

  StandAloneFsiLednicky *tWeightLed = new StandAloneFsiLednicky();
  StandAloneFsiKisiel *tWeightKisiel = new StandAloneFsiKisiel();
  if (pipi) {
    tWeightLed->setPairType(PiPlusPiPlus);
    tWeightLed->setCoulOff();
    tWeightLed->setQuantumOn();
    tWeightLed->setStrongOff();
    tWeightLed->set3BodyOff();
    
    tWeightKisiel->setPairType(PiPlusPiPlus);
    tWeightKisiel->setCoulOff();
    tWeightKisiel->setQuantumOn();
    tWeightKisiel->setStrongOff();
  }
  else {
    tWeightLed->setPairType(ProtonProton);
    tWeightLed->setCoulOn();
    tWeightLed->setQuantumOn();
    tWeightLed->setStrongOn();
    tWeightLed->set3BodyOff();
    
    tWeightKisiel->setPairType(ProtonProton);
    tWeightKisiel->setCoulOn();
    tWeightKisiel->setQuantumOn();
    tWeightKisiel->setStrongOn();
  }
  double tWLed, tWKis;
  for (int iterp=0; iterp<pmax; iterp++) {
    tPair = tPairManager->GetFitPair(iterp);
    tPair->SetPosition();
    tWLed = tWeightLed->getWeight(*tPair);
    tWKis = tWeightKisiel->getWeight(*tPair);
    if ((fabs((tWLed - tWKis)/tWLed) > 1.0e-2) || (isnan(tWLed))) {
      cout << "w l k d p " << tWLed << "   \t" << tWKis << "   \t" << (tWLed - tWKis) << "   \t" << ((tWLed - tWKis)/tWLed) << endl; 
      cout << "KStar out side long " << tPair->GetKStar() << " " << tPair->GetKStarTimeOutSign() << " " << tPair->GetKStarTimeSideSign() << " " << tPair->GetKStarTimeLongSign() << endl;
      tWeightLed->getWeight(*tPair);
      tWeightKisiel->PrintNext();
      tWeightKisiel->getWeight(*tPair);
    }
    else {
      cout << tWLed << "   \t" << tWKis << "   \t" << (tWLed - tWKis) << "   \t" << ((tWLed - tWKis)/tWLed) << "   \t" << tPair->GetKStar() << endl; 
//       tWeightKisiel->PrintNext();
//       tWeightKisiel->getWeight(*tPair);
    }

    wsumled->Fill(tPair->GetKStar(), tWLed);
    wsumkis->Fill(tPair->GetKStar(), tWKis);
    wsumden->Fill(tPair->GetKStar());
  }

  TFile *ofile = new TFile("wdiff.root","RECREATE");
  ofile->cd();
  wsumled->Write();
  wsumkis->Write();
  wsumden->Write();
}
