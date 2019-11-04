#include "stdlib.h"
#include "CFStorage.h"

CFStorage::CFStorage()
{
  mFileName = "CalcStore.root";
  OpenFile();
  for (Int_t i =0; i<4; i++)
   { 
     mDirPair[i] = 0x0;
     mDirCalc[i] = 0x0;  
   }
//   mPairName = "";
//   mCalcPar = "";
  mPairCount = 0;
}

CFStorage::CFStorage(const char *aFileName)
{
  mFileName = aFileName;
  OpenFile();
  
  for (Int_t i =0; i<4; i++)
   { 
     mDirPair[i] = 0x0;
     mDirCalc[i] = 0x0;  
   }

//   mPairName = "";
//   mCalcPar = "";
  mPairCount = 0;
}

CFStorage::~CFStorage()
{
  mStorage->Close();
}

int      
CFStorage::SetPairName(const char *aPairName)
{
  int tCurCount=0;

  while (tCurCount < mPairCount) {
    if (!(strcmp(mPairName[tCurCount], aPairName))) {
      return tCurCount;
    }
    tCurCount++;
  }

  if (mPairCount > 3)
    {
      PRINT_MESSAGE("Why do You want to store more than four pair types at the same time?");
      PRINT_MESSAGE("This is unsupported - Aborting.");
      exit(0);
    }
  
  mPairName[mPairCount] = aPairName;
  mCalcPar[mPairCount] = "";

  ChangeToCurrentDir(mPairCount);
  
  mPairCount++;
  
  return (mPairCount-1);
}

void      
CFStorage::SetCalculationParameters(int aPairNumber, int aStr, int aCoul, int aQS)
{
  TString tDir;
  tDir += "S";
  tDir += aStr;
  tDir += "C";
  tDir += aCoul;
  tDir += "Q";
  tDir += aQS;
  
  mCalcPar[aPairNumber] = tDir;
  
  ChangeToCurrentDir(aPairNumber);
}

void
CFStorage::GetStoredFunction(int aPairNumber, const char *aCFIdentifier, TVectorD **aContent, TVectorD **aError) 
  throw(CF_Storage_Exception)
{
  TDirectory *tDFun;
  
  if (!mStorage || !mStorage->IsOpen()) 
    throw CF_STORAGE_EXCEPTION_NO_FILE;
    
  if ( ((mDirPair[aPairNumber]) && (mDirCalc[aPairNumber]) ) == kFALSE)
   {
     ::Info("CFStorage::GetStoredFunction","%d %#x %#x",aPairNumber,mDirPair[aPairNumber],mDirCalc[aPairNumber]);
     throw CF_STORAGE_EXCEPTION_NO_ENTRY;
   }   
  TDirectory* dir = mDirCalc[aPairNumber];
  if (dir == 0x0)
   {
    ::Info("CFStorage::GetStoredFunction","dir is null");
    throw CF_STORAGE_EXCEPTION_NO_ENTRY;
   }
   
  tDFun = (TDirectory*)dir->Get(aCFIdentifier); 
  
  if ( tDFun == 0x0)
   {
    ::Info("CFStorage::GetStoredFunction","tDFun is null");
    throw CF_STORAGE_EXCEPTION_NO_ENTRY;
   } 

  *aContent = (TVectorD *) tDFun->Get("TVectorD;1");
  *aError   = (TVectorD *) tDFun->Get("TVectorD;2");
  
  if (!aContent)
    throw CF_STORAGE_EXCEPTION_NO_ENTRY;
}


void      
CFStorage::StoreFunction(int aPairNumber, const char *aCFIdentifier, TVectorD *aContent, TVectorD *aError) 
  throw(CF_Storage_Exception)
{
  TDirectory *tDFun;

  ChangeToCurrentDir(aPairNumber);
  if (!(tDFun = (TDirectory *) mDirCalc[aPairNumber]->Get(aCFIdentifier)))
    {
      tDFun = mDirCalc[aPairNumber]->mkdir(aCFIdentifier);
    }
  tDFun->cd();
  aContent->Write();
  aError->Write();
  mDirCalc[aPairNumber]->Write();
  mDirPair[aPairNumber]->Write();
  mStorage->Flush();
  //  mStorage->Close();
  //   delete mDirCalc;
  //   delete mDirPair;
  //   delete mStorage;
  //   OpenFile();
  //   ChangeToCurrentDir();
  
}

void CFStorage::OpenFile()
{
  mStorage = new TFile(mFileName.Data(), "UPDATE");
  if (!(mStorage && mStorage->IsOpen()))
    mStorage = new TFile(mFileName.Data(), "RECREATE");
}


void CFStorage::ChangeToCurrentDir(int aPairNumber)
{
  if (mPairCount)
   {
    if (mPairName[aPairNumber].IsNull() == kFALSE) 
     {
       mStorage->cd();

       if (mDirPair[aPairNumber] == 0x0)
        {
          mDirPair[aPairNumber] = (TDirectory *) mStorage->Get(mPairName[aPairNumber]);
          if ( mDirPair[aPairNumber] == 0x0 ) 
            {
              ::Info("\n\n\n","Creating %s", mPairName[aPairNumber].Data());
              mDirPair[aPairNumber] = mStorage->mkdir(mPairName[aPairNumber]);
            }
        }

       if (mCalcPar[aPairNumber].IsNull() == kFALSE) 
        {
          PRINT_DEBUG_3("Getting " << mCalcPar[aPairNumber].Data());
          if (mDirCalc[aPairNumber] == 0x0)
           {
             mDirCalc[aPairNumber] = (TDirectory *) mDirPair[aPairNumber]->Get(mCalcPar[aPairNumber]);
             if (mDirCalc[aPairNumber] ==0x0) 
              {
                ::Info("\n\n\n","Creating %s", mCalcPar[aPairNumber].Data());
                mDirCalc[aPairNumber] = mDirPair[aPairNumber]->mkdir(mCalcPar[aPairNumber]);
              }
           } 
        }
    }
  }  
}


