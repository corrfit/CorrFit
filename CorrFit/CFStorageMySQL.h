#ifndef _CORRFIT_CFSTORAGEMYSQL_H_
#define _CORRFIT_CFSTORAGEMYSQL_H_

#ifdef MYSQLSTORAGE

#include "TFile.h"
#include "CFGlobal.h"
#include "ReadPar.h"
#include "TVectorD.h"
#include <mysql.h>

enum _CF_Storage_Exception{
  CF_STORAGE_EXCEPTION_NO_ENTRY,
  CF_STORAGE_EXCEPTION_NO_FILE
};

typedef enum _CF_Storage_Exception CF_Storage_Exception;

class CFStorageMySQL 
{
 public:
  CFStorageMySQL();
  ~CFStorageMySQL();

  int       SetPairName(const char *aPairName);
  void      SetCalculationParameters(int aPairNumber, int aStr, int aCoul, int aQS);
  void      SetStoreParameters(int aPairNumber,
			       const char *aModelName,
			       int aNModelParameters,
			       int aNFitBin,
			       double aFitMin,
			       double aFitMax,
			       double aNormMin,
			       double aNormMax,
			       double aMomResScale,
			       const char *aPairHash,
			       int aPairCount);
  void      GetStoredFunction(int aPairNumber, 
			      double *aParameterValues,
			      TVectorD **aContent, 
			      TVectorD **aError) throw(CF_Storage_Exception);
  void      StoreFunction(int aPairNumber, 
			  double *aParameterValues,
			  TVectorD *aContent, 
			  TVectorD *aError) throw(CF_Storage_Exception);
  
 protected:
  MYSQL *connection;
  int    mPairSystemID[4];
  int    mInteractionID;
  int    mSourceModelID;
  int    mPairHashID[4];
  int    mFitParameterID;
  int    mPairCount;

  int    mNParameters;
  int    mNFitBins;
  STR    mPairName[4];
  STR    mDBHost;
  STR    mDBUser;
  STR    mDBPass;

  void ReadParameters();
  int  GetSingleID(const char *aTable, const char *aIDField, const char *aField, const char *aValue);
  int  GetSMParametersID(int aPairNumber, double *aParameterValues);
  
};

#endif

#endif
