#include "PairManager.h"
#include "StandAloneFsiKisiel.h"
#include "StandAloneFsiLednicky.h"
#include "SourceModelGausRInvJacobian.h"
#include "SourceModelExpGausRInv.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

int main(int argc, char **argv)
{
  PairManager *tPairManager = new PairManager("pairpipi.root");
  Pair *tPair;
//   double tParameters[1] = {6.0};

//   SourceModelGausRInvJacobian *tSourceModel = new SourceModelGausRInvJacobian();
  double tParameters[3] = {6.0,0.0,0.12};

  SourceModelExpGausRInv *tSourceModel = new SourceModelExpGausRInv();
  tSourceModel->SetModelParameters(tParameters);
  Pair::SetSourceModel(tSourceModel);

  for (int iter=0; iter<1000000; iter++) {
    tPair = tPairManager->ReadPair();
    tPairManager->StoreFitPair(tPair);
    tPair = tPairManager->GetFitPair(iter);
    tPair->InitRandVar();
  }

  TH2D *coulqsc = new TH2D("coulqsc","coulqsc",50, 0.0, 0.2, 50, 0.0, 20.0);
  TH2D *denoqsc = new TH2D("denoqsc","denoqsc",50, 0.0, 0.2, 50, 0.0, 20.0);

  TH3D *qsvscoul = new TH3D("qsvscoul", "qsvscoul", 50, 0.0, 2.0, 45, -15.0, 0.0, 50, -4.0*TMath::Pi(), 4.0*TMath::Pi());
  TH3D *coulvsal = new TH3D("coulvsal", "coulvsal", 
			    50, -4.0*TMath::Pi(), 4.0*TMath::Pi(),
			    50, 0.0, 8.0*TMath::Pi(),
			    50, 0.0, 0.2);

  TH3D *coulvsald = new TH3D("coulvsald", "coulvsald", 
			     50, -4.0*TMath::Pi(), 4.0*TMath::Pi(),
			     50, 0.0, 8.0*TMath::Pi(),
			     50, 0.0, 0.2);
  
  TH3D *coulvsalr = new TH3D("coulvsalr", "coulvsalr", 
			     50, -4.0*TMath::Pi(), 4.0*TMath::Pi(),
			     50, 0.0, 8.0*TMath::Pi(),
			     50, 0.0, 0.2);
  
  TH2D *qsvscoul2d = new TH2D("qsvscoul2d", "qsvscoul2d", 45, -15.0, 0.0, 50, -4.0*TMath::Pi(), 4.0*TMath::Pi()); 

  TH1D *wsumled = new TH1D("wsumled","wsumled",50, 0.0, 0.2);
  TH1D *wsumkis = new TH1D("wsumkis","wsumkis",50, 0.0, 0.2);
  
  wsumled->Sumw2();
  wsumkis->Sumw2();

  StandAloneFsiLednicky *tWeightLed = new StandAloneFsiLednicky();
  tWeightLed->setPairType(PiPlusPiPlus);
  tWeightLed->setCoulOff();
  tWeightLed->setQuantumOn();
  tWeightLed->setStrongOff();
  tWeightLed->set3BodyOff();

  StandAloneFsiKisiel *tWeightKisiel = new StandAloneFsiKisiel();
  tWeightKisiel->setPairType(PiPlusPiPlus);
  tWeightKisiel->setCoulOn();
  tWeightKisiel->setQuantumOn();
  tWeightKisiel->setStrongOff();

  double tWLed, tWKis;
  double tWC, tWQ;

  for (int iterp=0; iterp<1000000; iterp++) {
    tPair = tPairManager->GetFitPair(iterp);
    tPair->SetPosition();
    //    tWLed = tWeightLed->getWeight(*tPair);
    tWKis = tWeightKisiel->getWeight(*tPair);

    tWeightKisiel->setCoulOn();
    tWeightKisiel->setQuantumOff();
    tWC = tWeightKisiel->getWeight(*tPair);

    tWeightKisiel->setCoulOff();
    tWeightKisiel->setQuantumOn();
    tWQ = tWeightKisiel->getWeight(*tPair);

//     if ((fabs((tWLed - tWKis)/tWLed) > 1.0e-2) || (isnan(tWLed))) {
//       cout << tWLed << "   \t" << tWKis << "   \t" << (tWLed - tWKis) << "   \t" << ((tWLed - tWKis)/tWLed) << endl; 
//       cout << "KStar out side long " << tPair->GetKStar() << " " << tPair->GetKStarTimeOutSign() << " " << tPair->GetKStarTimeSideSign() << " " << tPair->GetKStarTimeLongSign() << endl;
//       tWeightKisiel->PrintNext();
//       tWeightKisiel->getWeight(*tPair);
//     }
//     else {
//       cout << tWLed << "   \t" << tWKis << "   \t" << (tWLed - tWKis) << "   \t" << ((tWLed - tWKis)/tWLed) << endl; 
//     }

    double tPx = tPair->p1().x+tPair->p2().x;
    double tPy = tPair->p1().y+tPair->p2().y;
    double tPz = tPair->p1().z+tPair->p2().z;
    double tE  = tPair->p1().t+tPair->p2().t;
    double tPt = tPx*tPx + tPy*tPy;
    double tMt = tE*tE - tPz*tPz;
    double tM  = sqrt(tMt - tPt);
    tMt = sqrt(tMt);
    tPt = sqrt(tPt);
    
    double tDX = tPair->x1().x-tPair->x2().x;
    double tDY = tPair->x1().y-tPair->x2().y;
    double mRLong     = tPair->x1().z-tPair->x2().z;
    double mDTime     = tPair->x1().t-tPair->x2().t;
  
    double mRTrans = tDX>0.? sqrt(tDX*tDX+tDY*tDY) : -1.*sqrt(tDX*tDX+tDY*tDY);
    double mROut = (tDX*tPx + tDY*tPy)/tPt;
    double mRSide = (-tDX*tPy + tDY*tPx)/tPt;
    double mRSidePairCMS = mRSide;
    
    // Lab -> LCMS -> PRF method
    double tBeta = tPz/tE;
    double tGamma = tE/tMt;
    double mRLongPairCMS = tGamma*(mRLong - tBeta* mDTime);
    double mDTimePairLCMS = tGamma*(mDTime - tBeta* mRLong);
    tBeta = tPt/tMt;
    tGamma = tMt/tM;
    double mROutPairCMS = tGamma*(mROut - tBeta* mDTimePairLCMS);
    double mDTimePairCMS = tGamma*(mDTimePairLCMS - tBeta* mROut);
    double mRStar = sqrt(mROutPairCMS*mROutPairCMS + mRSidePairCMS*mRSidePairCMS +
		  mRLongPairCMS*mRLongPairCMS);

    double mKStarRStar = tPair->GetKStarTimeOutSign() * mROutPairCMS +
      tPair->GetKStarTimeSideSign() * mRSidePairCMS +
      tPair->GetKStarTimeLongSign() * mRLongPairCMS;
    mKStarRStar /= 0.197327;

    coulqsc->Fill(tPair->GetKStar(), mRStar, tWKis);
    denoqsc->Fill(tPair->GetKStar(), mRStar, 1.0);

    wsumled->Fill(tPair->GetKStar(), tWLed);
    wsumkis->Fill(tPair->GetKStar(), tWKis);

    if (fabs(mKStarRStar) < 4.0*TMath::Pi()) {
      qsvscoul->Fill(tWQ, TMath::Log(1.0-tWC), mKStarRStar);
      qsvscoul2d->Fill( TMath::Log(1.0-tWC), mKStarRStar);
      cout << "KStar WKis " << " " << tPair->GetKStar() << " " << tWKis << endl;
      coulvsal ->Fill(mKStarRStar, tPair->GetKStar() * mRStar/0.197327, tPair->GetKStar(), tWKis);
      coulvsald->Fill(mKStarRStar, tPair->GetKStar() * mRStar/0.197327, tPair->GetKStar(), 1.0);
    }
  }

  coulvsalr = new TH3D(*coulvsal);
  coulvsalr->Divide(coulvsald);

  TFile *ofile = new TFile("wdiff.root","RECREATE");
  ofile->cd();
  wsumled->Write();
  wsumkis->Write();
  coulqsc->Write();
  denoqsc->Write();
  qsvscoul->Write();
  qsvscoul2d->Write();
  coulvsalr->Write();

  TH2D *ratqsc = new TH2D(*coulqsc);
  ratqsc->SetName("ratqsc");
  ratqsc->Divide(denoqsc);
  
  ratqsc->Write();

}
