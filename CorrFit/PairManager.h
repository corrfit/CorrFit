#ifndef _CORRFIT_PAIRMANAGER_H_
#define _CORRFIT_PAIRMANAGER_H_

#include "CFGlobal.h"
#include "PairReader.h"
#include "Pair.h"

//#define MAXNORPAIR 80000

class PairManager
{
 public: 
  PairManager();
  PairManager(const char *filename);
  ~PairManager();
  
  void SetPairFileName(const char *filename);
  void CloseFile();
  
  Pair* ReadPair();
  Int_t GetPairCount();

  Pair* GetFitPair(Int_t aNum);
  Pair* GetNormPair(Int_t aNum);
  Int_t GetFitPairCount();
  Int_t GetNormPairCount();

  void StoreFitPair(Pair *aPair);
  void StoreNormPair(Pair *aPair);
  char *GetFirstPairHash();

 private:
  PairReader *mPairReader;
  Pair **mFitPairs;
  Pair **mNormPairs;

  Int_t mFitPairCount;
  Int_t mNormPairCount;

  Int_t mFPBufSize;
  Int_t mNPBufSize;

  Int_t mPairsRead;
};

#endif 
