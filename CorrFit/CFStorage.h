#ifndef _CORRFIT_CFSTORAGE_H_
#define _CORRFIT_CFSTORAGE_H_

#include "TFile.h"
#include "CFGlobal.h"
#include "ReadPar.h"
#include "TVectorD.h"

enum _CF_Storage_Exception{
  CF_STORAGE_EXCEPTION_NO_ENTRY,
  CF_STORAGE_EXCEPTION_NO_FILE
};

typedef enum _CF_Storage_Exception CF_Storage_Exception;

class CFStorage 
{
 public:
  CFStorage();
  CFStorage(const char *aFileName);
  ~CFStorage();

  int       SetPairName(const char *aPairName);
  void      SetCalculationParameters(int aPairNumber, int aStr, int aCoul, int aQS);
  void      GetStoredFunction(int aPairNumber, const char *aCFIdentifier, TVectorD **aContent, TVectorD **aError) throw(CF_Storage_Exception);
  void      StoreFunction(int aPairNumber, const char *aCFIdentifier, TVectorD *aContent, TVectorD *aError) throw(CF_Storage_Exception);
  
 protected:
  TFile      *mStorage;
  TDirectory *mDirPair[4];
  TDirectory *mDirCalc[4];
  STR mPairName[4];
  STR mCalcPar[4];
  STR mFileName;
  int mPairCount;

  void OpenFile();
  void ChangeToCurrentDir(int aPairNumber);
  
};

#endif
