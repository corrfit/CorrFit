#include <stdlib.h>
#include "SourceModelGaus.h"

TRandom *SourceModelGaus::mRandom = new TRandom();

SourceModelGaus::SourceModelGaus()
{
  for (int iter=0; iter<8; iter++)
    mParameters[iter] = 0.0;
  mNPar = 2;
  mType = SOURCEMODELGAUS_ONEGAUS_ONERADIUS;
  mNRandVar = 8;
}

SourceModelGaus::~SourceModelGaus()
{
}

SourceModelGaus::SourceModelGaus(SourceModelGausType aType)
{
  for (int iter=0; iter<8; iter++)
    mParameters[iter] = 0.0;
  mType = aType;
  mNRandVar = 8;
  switch(mType){
  case SOURCEMODELGAUS_ONEGAUS_THREERADII:
    mNPar = 4;
    break;
  case SOURCEMODELGAUS_TWOGAUS_ONERADIUS:
    mNPar = 4;
    break;
  case SOURCEMODELGAUS_TWOGAUS_THREERADII:
    mNPar = 8;
  case SOURCEMODELGAUS_ONEGAUS_ONERADIUS:
  default:
    mNPar = 2;
    break;
  }
}

void       
SourceModelGaus::SetModelParameters(double *aParTable)
{
  for (int ti=0; ti<mNPar; ti++)
    mParameters[ti] = aParTable[ti];
}

double*
SourceModelGaus::GetModelParameters()
{
  return mParameters;
}

void        
SourceModelGaus::InitRandVar(double** aRandVar)
{
  (*aRandVar) = (double *) malloc(sizeof(double) * mNRandVar);
  for (int ti = 0; ti<mNRandVar-1; ti++) {
    (*aRandVar)[ti] = mRandom->Gaus(0.0,1.0);
  }
  (*aRandVar)[mNRandVar-1] = mRandom->Exp(1.0);
 switch (mType) 
   {
   case SOURCEMODELGAUS_TWOGAUS_ONERADIUS:
     (*aRandVar)[1] = mRandom->Exp(1.0);
     break;
   case SOURCEMODELGAUS_TWOGAUS_THREERADII:
     (*aRandVar)[3] = mRandom->Exp(1.0);
     break;
   }
}

void
SourceModelGaus::GeneratePos(const VectorPair *aMom, VectorPair *aPos, double *aRandVar)
{
  switch (mType) {
  case SOURCEMODELGAUS_ONEGAUS_ONERADIUS:
    aPos->part1.x = aRandVar[0]*mParameters[0];
    aPos->part1.y = aRandVar[1]*mParameters[0];
    aPos->part1.z = aRandVar[2]*mParameters[0];
    aPos->part1.t = aRandVar[3]*mParameters[1];
    aPos->part2.x = aRandVar[4]*mParameters[0];
    aPos->part2.y = aRandVar[5]*mParameters[0];
    aPos->part2.z = aRandVar[6]*mParameters[0];
    aPos->part2.t = aRandVar[7]*mParameters[1];
    break;
  case SOURCEMODELGAUS_ONEGAUS_THREERADII:
    aPos->part1.x = aRandVar[0]*mParameters[0];
    aPos->part1.y = aRandVar[1]*mParameters[1];
    aPos->part1.z = aRandVar[2]*mParameters[2];
    aPos->part1.t = aRandVar[3]*mParameters[3];
    aPos->part2.x = aRandVar[4]*mParameters[0];
    aPos->part2.y = aRandVar[5]*mParameters[1];
    aPos->part2.z = aRandVar[6]*mParameters[2];
    aPos->part2.t = aRandVar[7]*mParameters[3];
    break;
  case SOURCEMODELGAUS_TWOGAUS_ONERADIUS:
    aPos->part1.x = aRandVar[0]*mParameters[0];
    aPos->part1.y = aRandVar[1]*mParameters[0];
    aPos->part1.z = aRandVar[2]*mParameters[0];
    aPos->part1.t = aRandVar[3]*mParameters[1];
    aPos->part2.x = aRandVar[4]*mParameters[2];
    aPos->part2.y = aRandVar[5]*mParameters[2];
    aPos->part2.z = aRandVar[6]*mParameters[2];
    aPos->part2.t = aRandVar[7]*mParameters[3];
    break;
  case SOURCEMODELGAUS_TWOGAUS_THREERADII:
    aPos->part1.x = aRandVar[0]*mParameters[0];
    aPos->part1.y = aRandVar[1]*mParameters[1];
    aPos->part1.z = aRandVar[2]*mParameters[2];
    aPos->part1.t = aRandVar[3]*mParameters[3];
    aPos->part2.x = aRandVar[4]*mParameters[4];
    aPos->part2.y = aRandVar[5]*mParameters[5];
    aPos->part2.z = aRandVar[6]*mParameters[6];
    aPos->part2.t = aRandVar[7]*mParameters[7];
    break;
  }  
}

const char*       
SourceModelGaus::GetParameterName(int aPar)
{
  switch (mType) {
  case SOURCEMODELGAUS_ONEGAUS_ONERADIUS:
    switch (aPar) {
    case 0:
      return "Radius";
      break;
    case 1:
      return "Time spread";
      break;
    }
    break;
  case SOURCEMODELGAUS_ONEGAUS_THREERADII:
    switch (aPar) {
    case 0:
      return "Out Radius";
      break;
    case 1:
      return "Side Radius";
      break;
    case 2:
      return "Long Radius";
      break;
    case 3:
      return "Time spread";
      break;
    }
    break;
  case SOURCEMODELGAUS_TWOGAUS_ONERADIUS:
    switch (aPar) {
    case 0:
      return "Radius 1";
      break;
    case 1:
      return "Time spread 1";
      break;
    case 2:
      return "Radius 2";
      break;
    case 3:
      return "Time spread 2";
      break;
    }
    break;
  case SOURCEMODELGAUS_TWOGAUS_THREERADII:
    switch (aPar) {
    case 0:
      return "Out Radius 1";
      break;
    case 1:
      return "Side Radius 1";
      break;
    case 2:
      return "Long Radius 1";
      break;
    case 3:
      return "Time spread 1";
      break;
    case 4:
      return "Out Radius 2";
      break;
    case 5:
      return "Side Radius 2";
      break;
    case 6:
      return "Long Radius 2";
      break;
    case 7:
      return "Time spread 2";
      break;
    }
    break;
  }
}

const char* SourceModelGaus::GetModelIdentifier()
{
  switch(mType){
  case SOURCEMODELGAUS_ONEGAUS_ONERADIUS:
    return "SMGaus11";
    break;
  case SOURCEMODELGAUS_ONEGAUS_THREERADII:
    return "SMGaus13";
    break;
  case SOURCEMODELGAUS_TWOGAUS_ONERADIUS:
    return "SMGaus21";
    break;
  case SOURCEMODELGAUS_TWOGAUS_THREERADII:
    return "SMGaus23";
    break;
  }
}


