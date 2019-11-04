#ifndef _CORRFIT_PAIRREADER_H_
#define _CORRFIT_PAIRREADER_H_

#include <iostream>
#include <TFile.h>
#include <TNtuple.h>

using namespace std;
#define D_(_mes) cout << _mes << endl;

enum _Pair_Reader_Mode{
  PAIR_READER_READ,
  PAIR_READER_WRITE
};

typedef enum _Pair_Reader_Mode Pair_Reader_Mode;

enum _Pair_Reader_Exception{
  PAIR_READER_EXCEPTION_END_OF_ENTRIES,
  PAIR_READER_EXCEPTION_NO_FILE,
  PAIR_READER_EXCEPTION_NO_NTUPLE
};

typedef enum _Pair_Reader_Exception Pair_Reader_Exception;

class PairReader 
{
 public:
  PairReader();
  PairReader(const char *filename, Pair_Reader_Mode mode);
  ~PairReader();
  void             SetFileName(const char *filename);
  void             SetMode(Pair_Reader_Mode mode);
  Float_t         *ReadPair() throw(Pair_Reader_Exception);
  void             ReadPair(Float_t **aBuffer) throw(Pair_Reader_Exception);
  Int_t            GetPairCount() throw(Pair_Reader_Exception);
  void             WritePair(Float_t *) throw(Pair_Reader_Exception);
  void             WritePair(Float_t px1, Float_t py1, Float_t pz1, Float_t e1,
			     Float_t px2, Float_t py2, Float_t pz2, Float_t e2) throw(Pair_Reader_Exception);
  const char      *GetFileName();
  Pair_Reader_Mode GetMode();
  void CloseFile();
 protected:
  TFile           *mFile;
  TNtuple         *mNtuple;
  Pair_Reader_Mode mMode;
  
  Float_t          mRelBuffer[8];
  int              mRelocateTable[8];
  int              mRelocated;
  int              mCurrentIndex;
  void OpenFile(const char *filename);
  void SetRelocate(int number, const char *name);
};


#endif
