#include "PairReader.h"

PairReader::PairReader()
{
  mMode = PAIR_READER_READ;
  mFile = 0;
  mRelocated = 0;
}

PairReader::PairReader(const char *filename, Pair_Reader_Mode mode)
{
  mMode = mode;
  OpenFile(filename);
  mRelocated = 0;
}

PairReader::~PairReader()
{
//   mFile->cd();
//   if (mMode == PAIR_READER_WRITE)
//     mNtuple->Write();
//   mFile->Close();

  CloseFile();
//   if (mNtuple)
//     delete mNtuple;
  if (mFile)
    delete mFile;
}
  
void 
PairReader::SetFileName(const char *filename)
{
  if (mFile->IsOpen()) {
    mFile->Close();
  }
  OpenFile(filename);
}
    
void 
PairReader::SetMode(Pair_Reader_Mode mode)
{
  mMode = mode;
}

Float_t *
PairReader::ReadPair() throw (Pair_Reader_Exception)
{
  Float_t *tBuffer;
  
  if (!mCurrentIndex) {
    if (!mFile || !mFile->IsOpen()) 
      throw PAIR_READER_EXCEPTION_NO_FILE;
    if (!mNtuple)
      throw PAIR_READER_EXCEPTION_NO_NTUPLE;
  }
  else if (mCurrentIndex >= mNtuple->GetEntries())
    throw PAIR_READER_EXCEPTION_END_OF_ENTRIES;
  mNtuple->GetEntry(mCurrentIndex++);
  tBuffer = mNtuple->GetArgs();
  for (int ti=0; ti<8; ti++)
    mRelBuffer[ti] = tBuffer[mRelocateTable[ti]];
  for (int ti=0; ti<8; ti++)
    tBuffer[ti] = mRelBuffer[ti];
  return tBuffer;
}

void
PairReader::ReadPair(Float_t **aBuffer) throw (Pair_Reader_Exception)
{
  if (!mCurrentIndex) {
    if (!mFile || !mFile->IsOpen()) 
      throw PAIR_READER_EXCEPTION_NO_FILE;
    if (!mNtuple)
      throw PAIR_READER_EXCEPTION_NO_NTUPLE;
  }
  else if (mCurrentIndex >= mNtuple->GetEntries())
    throw PAIR_READER_EXCEPTION_END_OF_ENTRIES;
  mNtuple->GetEntry(mCurrentIndex++);
  (*aBuffer) = mNtuple->GetArgs();
  for (int ti=0; ti<8; ti++)
    mRelBuffer[ti] = (*aBuffer)[mRelocateTable[ti]];
  for (int ti=0; ti<8; ti++)
    (*aBuffer)[ti] = mRelBuffer[ti];
}

Int_t
PairReader::GetPairCount()  throw (Pair_Reader_Exception)
{
  if (!mCurrentIndex) {
    if (!mFile || !mFile->IsOpen()) 
      throw PAIR_READER_EXCEPTION_NO_FILE;
    if (!mNtuple)
      throw PAIR_READER_EXCEPTION_NO_NTUPLE;
  }
  return mNtuple->GetEntries();
}

void 
PairReader::WritePair(Float_t *aArg) throw (Pair_Reader_Exception)
{
  if (!mCurrentIndex) {
    if (!mFile || !mFile->IsOpen()) 
      throw PAIR_READER_EXCEPTION_NO_FILE;
    if (!mNtuple)
      throw PAIR_READER_EXCEPTION_NO_NTUPLE;
  }
  mNtuple->Fill(aArg[0], aArg[1], aArg[2], aArg[3],
		aArg[4], aArg[5], aArg[6], aArg[7]);
  mCurrentIndex++;
}

void 
PairReader::WritePair(Float_t px1, Float_t py1, Float_t pz1, Float_t e1,
		      Float_t px2, Float_t py2, Float_t pz2, Float_t e2) 
  throw (Pair_Reader_Exception)
{
  if (!mCurrentIndex) {
    if (!mFile || !mFile->IsOpen()) 
      throw PAIR_READER_EXCEPTION_NO_FILE;
    if (!mNtuple)
      throw PAIR_READER_EXCEPTION_NO_NTUPLE;
  }
  mNtuple->Fill(px1, py1, pz1, e1,
		px2, py2, pz2, e2);
  mCurrentIndex++;
}

const char *
PairReader::GetFileName()
{
  return mFile->GetName();
}

Pair_Reader_Mode 
PairReader::GetMode()
{
  return mMode;
}

void 
PairReader::OpenFile(const char *aFileName)
{
  mCurrentIndex = 0;
  if (mMode == PAIR_READER_READ) {
    mFile = new TFile(aFileName);
    if (mFile->IsOpen()) {
      mNtuple = (TNtuple *) mFile->Get("Pair");
      if (!mNtuple) {
	D_("No pairs tree in file " << aFileName);
      }

      TObjArray *branches = mNtuple->GetListOfBranches();
      for (int ti=0; ti<branches->GetEntries(); ti++)
	{
	  SetRelocate(ti, branches->At(ti)->GetName());
	}
    }
    else {
      D_("No file " << aFileName);
    }
  }
  else if (mMode == PAIR_READER_WRITE) {
    mFile = new TFile(aFileName, "recreate");
    mFile->cd();
    mNtuple = new TNtuple("pair", "pair", "px1:py1:pz1:e1:px2:py2:pz2:e2");
  }
}

void 
PairReader::CloseFile()
{
  mFile->cd();
  if (mMode == PAIR_READER_WRITE)
    mNtuple->Write();
  //  mFile->Close();
  
//   if (mFile)
//     delete mFile;
  mFile=0;
  D_("Pair file closed");
}


void PairReader::SetRelocate(int number, const char *name)
{
  if (!strcmp(name,"px1"))
    {
      D_("Relocating " << name << " " << number << " number as px1 0");
      mRelocateTable[0] = number;
      mRelocated &= (1 >> 0);
    }
  else if (!strcmp(name,"py1"))
    {
      D_("Relocating " << name << " " << number << " number as py1 1");
      mRelocateTable[1] = number;
      mRelocated &= (1 >> 1);
    }
  else if (!strcmp(name,"pz1"))
    {
      D_("Relocating " << name << " " << number << " number as pz1 2");
      mRelocateTable[2] = number;
      mRelocated &= (1 >> 2);
    }
  else if (!strcmp(name,"e1"))
    {
      D_("Relocating " << name << " " << number << " number as e1 3");
      mRelocateTable[3] = number;
      mRelocated &= (1 >> 3);
    }
  else if (!strcmp(name,"px2"))
    {
      D_("Relocating " << name << " " << number << " number as px2 4");
      mRelocateTable[4] = number;
      mRelocated &= (1 >> 4);
    }
  else if (!strcmp(name,"py2"))
    {
      D_("Relocating " << name << " " << number << " number as py2 5");
      mRelocateTable[5] = number;
      mRelocated &= (1 >> 5);
    }
  else if (!strcmp(name,"pz2"))
    {
      D_("Relocating " << name << " " << number << " number as pz2 6");
      mRelocateTable[6] = number;
      mRelocated &= (1 >> 6);
    }
  else if (!strcmp(name,"e2"))
    {
      D_("Relocating " << name << " " << number << " number as e2 7");
      mRelocateTable[7] = number;
      mRelocated &= (1 >> 7);
    }
  else {
    D_("Unknown branch");
  }  
}

