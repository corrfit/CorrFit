/***************************************************************************
 *
 * $Id: 
 *
 * Author: Laurent Conin, Fabrice Retiere, Subatech, France
 ***************************************************************************
 *
 * Description : Create pair from StHbtFsiHiddenInfo
 *
 ***************************************************************************
 *
 * $Log: 
 *
 ***************************************************************************/

#ifndef ST_HBT_DATA_PAIR_HH
#define ST_HBT_DATA_PAIR_HH

class TRandom2;
#include <math.h>
#include <iostream>
#include "CFGlobal.h"
#include "SourceModel.h"
#include "TRandom2.h"

struct MomResParameters{
  double PhiA;
  double PhiB;
  double PhiAlpha;
  double PhiMu;
  double ThetaA;
  double ThetaB;
  double ThetaAlpha;
  double PMeanA;
  double PMeanB;
  double PMeanAlpha;
  double PResA;
  double PResB;
  double PResAlpha;
  double PResC;
};

class Pair{ 
public:
  Pair();
  virtual ~Pair();
  Pair& operator=(Pair& aPair);

  virtual void SetMomentum(const float* aVect);
  virtual void Set(const float* aVect);

  static void SetSourceModel(SourceModel *aSourceModel);
  static void SetMomRes(double aMomResScale, int aPairSystem);
  static SourceModel* GetSourceModel();
  static double GetMomRes();
  
  void SetPairNum(int aNum);
  int  GetPairNum();
  static void GetSeed(UInt_t &s1, UInt_t &s2);
  static void SetSeed(UInt_t s1, UInt_t s2);

  static int nPar();

  void SetPosition();

  void SetRelocate(int number, const char *name);

  double GetKStarTimeOutSign();
  double GetKStarTimeSideSign();
  double GetKStarTimeLongSign();
  double GetKStar();
  double GetBetat();
  double GetKt();

  double GetKStarTimeOutSignSmeared();
  double GetKStarTimeSideSignSmeared();
  double GetKStarTimeLongSignSmeared();
  double GetKStarSmeared();
  
  double GetKStarOut();
  double GetKStarSide();
  double GetKStarLong();

  FourVector p1();
  FourVector p2();
  FourVector x1();
  FourVector x2();

  VectorPair PairMomentum();
  VectorPair PairPosition();
  
  int bad(){return mBad;}
  int GetBin();
  void SetBin(int aBin);  

  void SetPosMom();
  void SetPosMomEFirstNotation();
  void InitRandVar();
  
  static double *mom1;
  static double *mom2;
  static double *pos1;
  static double *pos2;

  VectorPair mMom;
  VectorPair mPos;

 private:
  int mBad;
  static TRandom2* mRand;

  static int mRelocateTable[8];
  static int mRelocated;
  static double mMomResScale;

  double* mRandVar;

  double mKStarSigned;

  double mKStarOut;
  double mKStarSide;
  double mKStarLong;
  
  int mBin;

  int mPairNum;
  
  static MomResParameters* GetPionMomResParameters();
  static MomResParameters* GetKaonMomResParameters();
  static MomResParameters* GetProtonMomResParameters();
  static MomResParameters* GetAntiProtonMomResParameters();
  static MomResParameters* mMomRes1;
  static MomResParameters* mMomRes2;  

  static SourceModel *mSourceModel;
};

inline int Pair::GetBin(){return mBin;}

inline void Pair::SetBin(int aBin){ mBin = aBin; }

inline double Pair::GetKStarTimeOutSign(){return mKStarOut;}
inline double Pair::GetKStarTimeSideSign(){return mKStarSide;}
inline double Pair::GetKStarTimeLongSign(){return mKStarLong; }
inline double Pair::GetKStar(){return mKStarSigned>0.0 ? mKStarSigned : -mKStarSigned;}

inline double Pair::GetKStarOut() { return mKStarOut; };
inline double Pair::GetKStarSide() { return mKStarSide; };
inline double Pair::GetKStarLong() { return mKStarLong; };

inline FourVector Pair::p1(){ return mMom.part1; }
inline FourVector Pair::p2(){ return mMom.part2; }
inline FourVector Pair::x1(){ return mPos.part1; }
inline FourVector Pair::x2(){ return mPos.part2; }


inline VectorPair Pair::PairMomentum(){ return mMom; }
inline VectorPair Pair::PairPosition(){ return mPos; }

#endif
