#ifndef _CORRFIT_PAIRSYSTEMS_H_
#define _CORRFIT_PAIRSYSTEMS_H_

enum _CorrFitPairType {
  NeutronNeutron       = 1,
  ProtonProton         = 2,
  NeutronProton        = 3,
  AlfaAlfa             = 4,
  PiPlusPiMinus        = 5,
  PiZeroPiZero         = 6,
  PiPlusPiPlus         = 7,
  NeutronDeuteron      = 8,
  ProtonDeuteron       = 9,
  PiPlusKaonMinus      = 10,
  PiPlusKaonPlus       = 11,
  PiPlusProton         = 12,
  PiMinusProton        = 13,
  KaonPlusKaonMinus    = 14,
  KaonPlusKaonPlus     = 15,
  KaonPlusProton       = 16,
  KaonMinusProton      = 17,
  DeuteronDeuteron     = 18,
  DeuteronAlfa         = 19,
  TritonTriton         = 20,
  TritonAlfa           = 21,
  KaonZeroKaonZero     = 22,
  KaonZeroAntiKaonZero = 23,
  DeuteronTriton       = 24,
  ProtonTriton         = 25,
  ProtonAlfa           = 26,
  ProtonLambda         = 27,
  NeutronLambda        = 28,
  LambdaLambda         = 29,
  ProtonAntiProton     = 30
};

typedef enum _CorrFitPairType CorrFitPairType;

enum _CorrFitPairSystem {
  PiplusKplus    = 1,
  PiminusKminus  = 2,
  PiplusKminus   = 3,
  PiminusKplus   = 4,
  PiplusP        = 11,
  PiminusAntiP   = 12,
  PiplusAntiP    = 13,
  PiminusP       = 14,
  KplusP         = 21,
  KminusAntiP    = 22,
  KplusAntiP     = 23,
  KminusP        = 24,
  PiplusPiplus   = 31,
  PiminusPiminus = 32,
  PiplusPiminus  = 33,
  KplusKplus     = 41,
  KminusKminus   = 42,
  KplusKminus    = 43,
  KzeroKzero     = 45,
  PP             = 51,
  PAntiP         = 53,
  NN             = 55,
  NP             = 56,
  PD             = 61,
  ND             = 66,
  PT             = 71,
  PA             = 81,
  DD             = 91,
  DT             = 101,
  DA             = 111,
  TT             = 121,
  TA             = 131,
  AA             = 141,
  KzeroSKzeroS   = 155,
  LastSystem     = 156
};

typedef enum _CorrFitPairSystem CorrFitPairSystem;

#endif
