#include "CFGlobal.h"
#include "Chi2Map.h"
#include <TFile.h>
#include <TGraphErrors.h>
#include <stdlib.h>
#include <iostream>
#ifdef RANDOMC2MORDER
#include <sys/time.h>
#include <stdlib.h>
#endif

extern TFile   *sOutFile;

using namespace std;

Chi2Map::Chi2Map()
{
  mDimension = 0;
  mMapStorage = NULL;
  mPurityBins = 1;
  mPurityMax = 1.0;
  mPurityMin = 1.0;
  mNormBins = 1;
  mNormMax = 1.0;
  mNormMin = 1.0;
  mChiMin = 1e30;
  mParMins = 0;
  mParMaxs = 0;
  mParBins = 0;
  mParSteps = 0;
  mBestPars = 0;
  mCurCoords = 0;
}

Chi2Map::Chi2Map(Int_t aDimension, Int_t *aNBins, Double_t *aBinMins, Double_t *aBinMaxs)
{
  mChiMin = 1e30;
  mMapStorage = NULL;
  mPurityBins   = 1;
  mPurityMax    = 1.0;
  mPurityMin    = 1.0;
  mNormBins = 1;
  mNormMax = 1.0;
  mNormMin = 1.0;
  mParMins = 0;
  mParMaxs = 0;
  mParBins = 0;
  mParSteps = 0;
  mBestPars = 0;
  mCurCoords = 0;
  SetMapParams(aDimension, aNBins, aBinMins, aBinMaxs);
}

Chi2Map::~Chi2Map()
{
  if (mMapStorage) free(mMapStorage);
  if (mParMins)    free(mParMins);
  if (mParMaxs)    free(mParMaxs);
  if (mParBins)    free(mParBins);
  if (mParSteps)   free(mParSteps);
  if (mBestPars)   free(mBestPars);
  if (mCurCoords)  free(mCurCoords);
}

void      
Chi2Map::SetPurityParams(Int_t aPurityBins, Double_t aPurityMin, Double_t aPurityMax)
{
  mPurityBins = aPurityBins;
  mPurityMin  = aPurityMin;
  mPurityMax  = aPurityMax;
  mPurityStep = (mPurityBins > 1) ? (mPurityMax - mPurityMin) / (mPurityBins - 1.0) : 0.0;
}

void      
Chi2Map::SetNormParams(Int_t aNormBins, Double_t aNormMin, Double_t aNormMax)
{
  mNormBins = aNormBins;
  mNormMin  = aNormMin;
  mNormMax  = aNormMax;
  mNormStep = (mNormBins > 1) ? (mNormMax - mNormMin) / (mNormBins - 1.0) : 0.0;
}

void      
Chi2Map::SetMapParams(Int_t aDimension, Int_t *aNBins, Double_t *aBinMins, Double_t *aBinMaxs)
    {
  if (mParBins && (aDimension!=mDimension)) {
    // Resetting parameters. Free the previous ones
    if (mMapStorage) free(mMapStorage);
    mMapStorage = NULL;
    if (mParMins)    free(mParMins);
    if (mParMaxs)    free(mParMaxs);
    if (mParBins)    free(mParBins);
    mParBins = NULL;
    if (mParSteps)   free(mParSteps);
    if (mBestPars)   free(mBestPars);
  }
  mDimension = aDimension;
  if (!mParBins) {
    mBestPars     = (Double_t *) malloc(sizeof(Double_t) * mDimension);
    mParMins      = (Double_t *) malloc(sizeof(Double_t) * mDimension);
    mParMaxs      = (Double_t *) malloc(sizeof(Double_t) * mDimension); 
    mParSteps     = (Double_t *) malloc(sizeof(Double_t) * mDimension); 
    mParBins      = (Int_t *)    malloc(sizeof(Int_t) * mDimension);
  }

  for (int ti=0; ti<mDimension; ti++)
    {
      mBestPars[ti] = 1e6;
      mParMins[ti]  = aBinMins[ti];
      mParMaxs[ti]  = aBinMaxs[ti];
      mParBins[ti]  = aNBins[ti];
      mParSteps[ti] = (aNBins[ti] > 1) ? (aBinMaxs[ti] - aBinMins[ti]) / (aNBins[ti] - 1.0) : 0.0;
    }
}
 
void 
Chi2Map::InitFill() throw(CorrFit_Chi2Map_Exception)
{
  if (!mParBins) throw (CHI2MAP_EXCEPTION_PARAMETERS_NOT_SET);
  PRINT_DEBUG("Initializing Chi2 map of " << mDimension << " dimension(s): ");
  
  // Create the map controlling variables
  mCurCoords    = (int *) malloc(sizeof(int) * mDimension);
  mCurCoordsInc = (int *) malloc(sizeof(int) * mDimension);

  // Generate a coarse Chi2 Map

  mMapSize = 1;
  for (int ti=0; ti<mDimension; ti++) {
    mCurCoordsInc[ti] = mMapSize;
    mMapSize *= mParBins[ti];
    mCurCoords[ti] = 0;
  }
  mCurCoords[0] = -1;
  mMapStorage = (double *) malloc(sizeof(double) * mMapSize * mPurityBins * mNormBins);

#ifdef RANDOMC2MORDER
  // Initialize a random filling order for the C2M

  struct timeval tv;
  int RAND_THIRD = RAND_MAX/3;

  int res = gettimeofday(&tv,NULL);
  srand(tv.tv_usec);
  
  list<int>::iterator tInsert = mOrder.begin();
  /*
  int rnum;
  int begnum = 1.0*rand()*mMapSize/RAND_MAX;
  mOrder.push_back(begnum);

  for (int iter=1; iter<mMapSize; iter++)
    {
      rnum = rand();
      if (rnum < (RAND_THIRD))
	mOrder.push_back((iter+begnum) % mMapSize);
      else if (rnum <(RAND_THIRD*2))
	mOrder.push_front((iter+begnum) % mMapSize);
      else
	mOrder.insert(tMiddle,(iter+begnum) % mMapSize);
    }
  */
  int ndiff;
  for (int iter=0; iter<mMapSize; iter++) {
    ndiff = rand() % (mMapSize/2 + 1) ;
    while (ndiff) {
      tInsert++;
      ndiff--;
      if (tInsert == mOrder.end())
	tInsert = mOrder.begin();
    }
    mOrder.insert(tInsert, iter);
  }

  mOrderCur = mOrder.begin();
#endif

  mSubMapSizes = (int *) malloc(sizeof(int) * mDimension);
  mSubMapSizes[mDimension-1] = 1;
  for (int ti=mDimension-2; ti>=0; ti--)
    mSubMapSizes[ti] = mSubMapSizes[ti+1]*mParBins[ti+1];

}
 
void      
Chi2Map::SetCellContent(Int_t *aCoords, Double_t *values) throw (CorrFit_Chi2Map_Exception)
{
  Int_t tMapPlace;
  
  if (!mMapStorage) {
    if (!mParBins) throw (CHI2MAP_EXCEPTION_PARAMETERS_NOT_SET);
    else throw (CHI2MAP_EXCEPTION_NO_MAP);
  }
  tMapPlace = 0;
  PRINT_DEBUG_3("Setting cell contents for coords");
  for (Int_t ti=0; ti<mDimension; ti++) {
    tMapPlace += aCoords[ti]*mCurCoordsInc[ti];
    PRINT_DEBUG_3(ti << " " << aCoords[ti] << " " << tMapPlace);
  }
  for (int tN=0; tN<mNormBins; tN++) {
    for (int tP=0; tP<mPurityBins; tP++) {
      mMapStorage[tMapPlace+mMapSize*(tP+tN*mPurityBins)] = values[tP];
      if (mChiMin > values[(tP+tN*mPurityBins)]) {
	mChiMin = values[(tP+tN*mPurityBins)];
	for (Int_t ti=0; ti<mDimension; ti++) {
	  mBestPars[ti] = mParMins[ti] + aCoords[ti]*mParSteps[ti];
	}    
	mBestPurity = mPurityMin+mPurityStep*tP;
	mBestNorm = mNormMin+mNormStep*tN;
	mBestChi = mChiMin;
	
	PRINT_MESSAGE("Setting minimum to:");
	for (int ti=0; ti<mDimension; ti++)
	  PRINT_MESSAGE(mBestPars[ti] << " ");
	PRINT_MESSAGE((mPurityMin+mPurityStep*tP) << " " << mChiMin);
      }
      
    }
  }
}

void      
Chi2Map::GetBestParams(Double_t *aParams, Double_t *aPurity, Double_t *aNorm, Double_t *aChiMin) throw(CorrFit_Chi2Map_Exception)
{
  if (!mMapStorage) {
    if (!mParBins) throw (CHI2MAP_EXCEPTION_PARAMETERS_NOT_SET);
    else throw (CHI2MAP_EXCEPTION_NO_MAP);
  }

  for (Int_t ti=0; ti<mDimension; ti++)
    aParams[ti] = mBestPars[ti];
  *aPurity = mBestPurity;
  *aNorm = mBestNorm;
  *aChiMin = mBestChi;
}


void      
Chi2Map::GetNextCellParams(Double_t *aPars) throw(CorrFit_Chi2Map_Exception)
{
  if (!mMapStorage) {
    if (!mParBins) throw (CHI2MAP_EXCEPTION_PARAMETERS_NOT_SET);
    else throw (CHI2MAP_EXCEPTION_NO_MAP);
  }

#ifdef RANDOMC2MORDER
  if (mOrderCur == mOrder.end())
    throw (CHI2MAP_EXCEPTION_END_OF_MAP);  
  int tCurCell = *mOrderCur;
  mOrderCur++;
  PRINT_DEBUG_3("Setting cell " << tCurCell);
  for (int ti=0; ti<mDimension; ti++)
    {
      mCurCoords[ti] = tCurCell / mSubMapSizes[ti];
      tCurCell = tCurCell % mSubMapSizes[ti];
      PRINT_DEBUG_3("Coordinate " << ti << "  " << mCurCoords[ti]);
    }
  
  
#else
  Int_t theEnd = 0;
  Int_t tinc   = 0;
  while ((++mCurCoords[tinc] >= mParBins[tinc]) && !theEnd) {
    mCurCoords[tinc] = 0;
    tinc++;
    if (tinc >= mDimension) {
      theEnd = 1;
      tinc = 0;
      mCurCoords[tinc] = 0;
    }
    }
  if (theEnd) throw (CHI2MAP_EXCEPTION_END_OF_MAP);
#endif
  
  for (Int_t ti=0; ti<mDimension; ti++)
    {
      aPars[ti] = mParMins[ti] + mCurCoords[ti]*mParSteps[ti];
    }
}

void      
Chi2Map::SetCurrentCellContent(Double_t *values) throw(CorrFit_Chi2Map_Exception)
{
  Int_t tMapPlace;
  TString tOutput;
  
  if (!mMapStorage) {
    if (!mParBins) throw (CHI2MAP_EXCEPTION_PARAMETERS_NOT_SET);
    else throw (CHI2MAP_EXCEPTION_NO_MAP);
  }
  tMapPlace = 0;
  for (Int_t ti=0; ti<mDimension; ti++) {
    tMapPlace += mCurCoords[ti]*mCurCoordsInc[ti];
  }
  for (int tP=0; tP<mPurityBins; tP++) {
    mMapStorage[tMapPlace+mMapSize*tP] = values[tP];
    if (mChiMin > values[tP]) {
      mChiMin = values[tP];
      for (Int_t ti=0; ti<mDimension; ti++) {
	mBestPars[ti] = mParMins[ti] + mCurCoords[ti]*mParSteps[ti];
      }    
      mBestPurity = mPurityMin+mPurityStep*tP;
      mBestChi = mChiMin;
      
      tOutput = "Setting minimum to: [";
      for (int ti=0; ti<mDimension; ti++) {
        if(ti) tOutput += ", ";
	tOutput += Form("%.3lf", mBestPars[ti]);
      }
      tOutput += "] for Purity: ";
      tOutput += Form("%.2lf", mPurityMin+mPurityStep*tP);
      tOutput += " and Chi2: ";
      tOutput += Form("%.2lf", mChiMin);
      PRINT_MESSAGE(tOutput);
    }
    
  }
}

void      
Chi2Map::MinuitMinimize() throw(CorrFit_Chi2Map_Exception)
{
}

int      
Chi2Map::GetMapSize() 
{
  return mMapSize;
}

// Input/output functions
void      
Chi2Map::WriteHisto()
{
  char tHName[30];
  int tBestPurBin=0;

  TGraphErrors *tG;
  double tBinXHalf, tBinYHalf;
  double tPurStep = (mPurityBins>1) ? (mPurityMax - mPurityMin) / (mPurityBins - 1.0) : 0.0;
  
  if (mDimension > 1) {
    TH2D *tH[mPurityBins];

    if (mParBins[0] > 1)
      tBinXHalf = (mParMaxs[0] - mParMins[0]) / (mParBins[0] - 1) / 2;
    else
      tBinXHalf = 0.01;
    if (mParBins[1] > 1)
      tBinYHalf = (mParMaxs[1] - mParMins[1]) / (mParBins[1] - 1) / 2;
    else
      tBinYHalf = 0.01;

    int tBeginX, tBeginY, tMultX, tMultY;

    sOutFile->cd();
    for (int tP=0; tP<mPurityBins; tP++) {
      sprintf(tHName, "Chi2MapP%.2lf", mPurityMin + tPurStep * tP);

      // The TH2D must have the smaller value as minimum
      // and larger value as maximum
      if (tBinXHalf > 0.0) {
	tBeginX = 1;
	tMultX = 1;
      }
      else {
	tBeginX = mParBins[0];
	tMultX = -1;
      }
      if (tBinYHalf > 0.0) {
	tBeginY = 1;
	tMultY = 1;
      }
      else {
	tBeginY = mParBins[1];
	tMultY = -1;
      }
      
      tH[tP] = new TH2D(tHName, tHName, 
			mParBins[0], 
			(tBinXHalf > 0) ? mParMins[0]-tBinXHalf : mParMaxs[0]+tBinXHalf, 
			(tBinXHalf > 0) ? mParMaxs[0]+tBinXHalf : mParMins[0]-tBinXHalf, 
			mParBins[1], 
			(tBinYHalf > 0) ? mParMins[1]-tBinYHalf : mParMaxs[1]+tBinYHalf, 
			(tBinYHalf > 0) ? mParMaxs[1]+tBinYHalf : mParMins[1]-tBinYHalf);
      for (int tx=0; tx<mParBins[0]; tx++)
	for (int ty=0; ty<mParBins[1]; ty++) {
	  PRINT_DEBUG_3("Setting cell " << tx << " " << ty << " to bin " << (tx+1) << " " <<(ty+1) << " " << (ty*mSubMapSizes[0]+tx*mSubMapSizes[1]));
	  tH[tP]->SetBinContent(tBeginX+tMultX*tx, tBeginY+tMultY*ty, mMapStorage[tx*mCurCoordsInc[0]+ty*mCurCoordsInc[1]+tP*mMapSize]);
	}
      
      tH[tP]->Write();
      if ((mPurityMin+tPurStep*tP) == mBestPurity) tBestPurBin = tP;
    }
  }
  else {
    TH1D *tH[mPurityBins];
    int tBegin, tMult;

    if (mParBins[0] > 1)
      tBinXHalf = (mParMaxs[0] - mParMins[0]) / (mParBins[0] - 1) / 2;
    else
      tBinXHalf = 0.01;

      if (tBinXHalf > 0.0) {
	tBegin = 1;
	tMult = 1;
      }
      else {
	tBegin = mParBins[1];
	tMult = -1;
      }

    sOutFile->cd();
    for (int tP=0; tP<mPurityBins; tP++) {
      sprintf(tHName, "Chi2MapP%.2lf", mPurityMin + tPurStep * tP);
      tH[tP] = new TH1D(tHName, tHName, mParBins[0], 
			(tBinXHalf > 0) ? mParMins[0]-tBinXHalf : mParMaxs[0]+tBinXHalf, 
			(tBinXHalf > 0) ? mParMaxs[0]+tBinXHalf : mParMins[0]-tBinXHalf);
      for (int tx=0; tx<mParBins[0]; tx++) {
	PRINT_DEBUG_3("Setting cell " << tx << " to bin " << (tx+1) << " " << (tx*mSubMapSizes[1]));
	tH[tP]->SetBinContent(tx+1, mMapStorage[tx*mCurCoordsInc[0]+tP*mMapSize]);
      }
      
      tH[tP]->Write();
      if ((mPurityMin+tPurStep*tP) == mBestPurity) tBestPurBin = tP;
      PRINT_DEBUG_3("Best purity bin " << tBestPurBin);
    }
  }
  if (mDimension == 3) {
    double tBinX, tBinY, tBinZ;
    int    tBestX, tBestY, tBestZ;

    if (mParBins[0] > 1) {
      tBinX = (mParMaxs[0] - mParMins[0]) / (mParBins[0] - 1);
      tBestX = (mBestPars[0] - mParMins[0])/tBinX;
    }
    else {
      tBinX = 0.01;
      tBestX = 0;
    }

    if (mParBins[1] > 1) {
      tBinY = (mParMaxs[1] - mParMins[1]) / (mParBins[1] - 1);
      tBestY = (mBestPars[1] - mParMins[1])/tBinY;
    }
    else {
      tBinY = 0.01;
      tBestY = 0;
    }

    if (mParBins[2] > 1) {
      tBinZ = (mParMaxs[2] - mParMins[2]) / (mParBins[2] - 1);
      tBestZ = (mBestPars[2] - mParMins[2])/tBinZ;
    }
    else {
      tBinZ = 0.01;
      tBestZ = 0;
    }

    double *xy0;
    double *yk0;
    double *xy1;
    double *yk1;
    double *xy2;
    double *yk2;
    
    xy0 = (double *) malloc(sizeof(double) * mParBins[0]);
    yk0 = (double *) malloc(sizeof(double) * mParBins[0]);

    for (int iterx=0; iterx<mParBins[0]; iterx++) {
      xy0[iterx] = mParMins[0] + iterx*tBinX;
      yk0[iterx] = mMapStorage[tBestPurBin*mMapSize + 
			       iterx *mCurCoordsInc[0] + 
			       tBestY*mCurCoordsInc[1] + 
			       tBestZ*mCurCoordsInc[2]];
    }

    TGraph *grx = new TGraph(mParBins[0], xy0, yk0);
    grx->SetName("CoordsMin0");
    grx->SetTitle("Chi2 plot along x axis");

    xy1 = (double *) malloc(sizeof(double) * mParBins[1]);
    yk1 = (double *) malloc(sizeof(double) * mParBins[1]);

    for (int iterx=0; iterx<mParBins[1]; iterx++) {
      xy1[iterx] = mParMins[1] + iterx*tBinY;
      yk1[iterx] = mMapStorage[tBestPurBin*mMapSize + 
			       tBestX*mCurCoordsInc[0] + 
			       iterx *mCurCoordsInc[1] + 
			       tBestZ*mCurCoordsInc[2]];
    }

    TGraph *gry = new TGraph(mParBins[1], xy1, yk1);
    gry->SetName("CoordsMin1");
    gry->SetTitle("Chi2 plot along y axis");

    xy2 = (double *) malloc(sizeof(double) * mParBins[2]);
    yk2 = (double *) malloc(sizeof(double) * mParBins[2]);

    for (int iterx=0; iterx<mParBins[2]; iterx++) {
      xy2[iterx] = mParMins[2] + iterx*tBinZ;
      yk2[iterx] = mMapStorage[tBestPurBin*mMapSize + 
			       tBestX*mCurCoordsInc[0] + 
			       tBestY*mCurCoordsInc[1] + 
			       iterx *mCurCoordsInc[2]];
    }

    TGraph *grz = new TGraph(mParBins[2], xy2, yk2);
    grz->SetName("CoordsMin2");
    grz->SetTitle("Chi2 plot along z axis");

    grx->Write();
    gry->Write();
    grz->Write();
  }
}



