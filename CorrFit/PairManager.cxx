#include <stdlib.h>
#include "PairManager.h"

PairManager::PairManager()
{
  mFitPairCount = 0;
  mNormPairCount = 0;
  mFPBufSize = 0;
  mNPBufSize = 0;
  mFitPairs = NULL;
  mNormPairs = NULL;
  mPairReader = NULL;
  mPairsRead = 0;
}

PairManager::PairManager(const char *filename)
{
  mPairReader = new PairReader(filename, PAIR_READER_READ);
  mFitPairs = NULL;
  mNormPairs = NULL;
  mFitPairCount = 0;
  mNormPairCount = 0;
  mFPBufSize = 0;
  mNPBufSize = 0;
  mPairsRead = 0;
}

PairManager::~PairManager()
{
  if (mPairReader) delete mPairReader;
  if (mFitPairs) {
    for (int iter=0; iter<mFitPairCount; iter++)
      delete mFitPairs[iter];
    delete [] mFitPairs;
  }
  if (mNormPairs) {
    for (int iter=0; iter<mNormPairCount; iter++)
      delete mNormPairs[iter];
    delete [] mNormPairs;
  }
}

Pair* PairManager::GetFitPair(Int_t aNum)
{
  return mFitPairs[aNum];
}

Pair* PairManager::GetNormPair(Int_t aNum)
{
  return mNormPairs[aNum];
}

Int_t PairManager::GetFitPairCount() 
{
  return mFitPairCount;
}

Int_t PairManager::GetNormPairCount()
{
  return mNormPairCount;
}

Pair* PairManager::ReadPair()
{
  Float_t *tPair;
  Pair *tRealPair;
  if (!mPairReader) 
    {
      D_("PairManager::ReadPair - Reading pairs without a Reader!");
      return NULL;
    }

  try {
//     tPair = mPairReader->ReadPair();
    mPairReader->ReadPair(&tPair);
    
    tRealPair = new Pair();
    tRealPair->SetMomentum(tPair);
    tRealPair->SetPairNum(mPairsRead++);
  }
  catch (Pair_Reader_Exception pe) {
    D_("PairManager::ReadPairs - Caught " << pe);
  }

  return tRealPair;
}

Int_t PairManager::GetPairCount()
{
  return mPairReader->GetPairCount();
}

void 
PairManager::StoreFitPair(Pair *aPair)
{
  //  PRINT_DEBUG("Size of pair class is " << sizeof(Pair));
  
  if (!mFPBufSize) {
    mFPBufSize = mPairReader->GetPairCount() / 5;
    mFitPairs = (Pair **) malloc(sizeof(Pair *) * mFPBufSize);
  }
  
  if (mFitPairCount >= mFPBufSize) {
    mFPBufSize += mFPBufSize/2;
    mFitPairs = (Pair **) realloc(mFitPairs,sizeof(Pair *) * mFPBufSize);
  }
      
  mFitPairs[mFitPairCount++] = aPair;
}

void
PairManager::StoreNormPair(Pair *aPair)
{
  if (!mNPBufSize) {
    mNPBufSize = mPairReader->GetPairCount() / 5;
    mNormPairs = (Pair **) malloc(sizeof(Pair *) * mNPBufSize);
  }
  
  if (mNormPairCount >= mNPBufSize) {
    mNPBufSize += mNPBufSize/2;
    mNormPairs = (Pair **) realloc(mNormPairs,sizeof(Pair *) * mNPBufSize);
  }

//   if (mNormPairCount < )
  mNormPairs[mNormPairCount++] = aPair;
}
			   
char      *
PairManager::GetFirstPairHash()
{
  Pair  *tPair;
  double tHash;
  char  *ret;
  
  ret = (char *) malloc(sizeof(char)*20);
  tPair = GetFitPair(0);
  tHash = (tPair->p1().x+2*tPair->p1().y+tPair->p1().z*tPair->p1().t -
	  tPair->p2().x*tPair->p2().y+3*tPair->p2().z+tPair->p2().t);
  tPair = GetFitPair(mFitPairCount-2);
  tHash += (tPair->p1().x+2*tPair->p1().y+tPair->p1().z*tPair->p1().t -
	  tPair->p2().x*tPair->p2().y+3*tPair->p2().z+tPair->p2().t);
  tPair = GetNormPair(mNormPairCount-3);
  tHash += (tPair->p1().x+2*tPair->p1().y+tPair->p1().z*tPair->p1().t -
	  tPair->p2().x*tPair->p2().y+3*tPair->p2().z+tPair->p2().t);

  sprintf(ret, "%lf\0", tHash);
  return ret;
}

void 
PairManager::SetPairFileName(const char *filename)
{
  mPairReader = new PairReader(filename, PAIR_READER_READ);
}

void 
PairManager::CloseFile()
{
  mPairReader->CloseFile();
}

