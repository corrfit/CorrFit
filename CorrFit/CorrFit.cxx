#include "fstream"
#include "CorrFit.h"
#include "CFFitter1DHBT.h"
#include "CFFitterNonId.h"
#include "CFFitterNonIdMult.h"
#include "CFFitterQuad.h"
#include "CFFitterSH.h"
#include "SourceModelCMSHyper.h"
#include "SourceModelCMSHyperBlast.h"
#include "SourceModelCMSHyperBt.h"
#include "SourceModelCMSHyperExpRQMD.h"
#include "SourceModelCMSHyperRQMD.h"
#include "SourceModelEHRInvShift.h"
#include "SourceModelExpGausRInv.h"
#include "SourceModelGaus.h"
#include "SourceModelGausCMS.h"
#include "SourceModelGausLCMS.h"
#include "SourceModelGausRInv.h"
#include "SourceModelGausRInvJacobian.h"
#include "SourceModelGausRInvShift.h"
#include "SourceModelGausROut.h"
#include "SourceModelGausROutDouble.h"
#include "SourceModelGausROutFile.h"
#include "SourceModelGausROutSourceCMS.h"
#include "SourceModelGausROutTwo.h"
#include "SourceModelGausRStar.h"
#include "ReadPar.h"
#include "StandAloneFsiLednicky.h"
#include "StandAloneFSIPratt.h"
#include "StandAloneFsiKisiel.h"

using namespace std;

ReadPar *sRPInstance;
STR      sRPFileName;
TFile   *sOutFile;
STR      sOutFileName;
STR      sCFFitterName;
STR      sSourceModelName;
STR      sWeightCalcName;

void ReadParameters()
{
  try {
    sOutFileName = sRPInstance->getPar("OutResultFileNamePrefix");
    PRINT_MESSAGE("of " << (sOutFileName.Data()));
    sOutFileName    += ".root";
    PRINT_MESSAGE("of " << (sOutFileName.Data()));
    sCFFitterName    = sRPInstance->getPar("CorrelationFitter");
    sSourceModelName = sRPInstance->getPar("SourceModel");
    sWeightCalcName  = sRPInstance->getPar("FsiCalculator");
  }
  catch (int e) {
    PRINT_DEBUG("CorrFit::ReadParameters - Caught exception " << e);
    PRINT_MESSAGE("Did not find one of the neccessary parameters in the parameters file.");
    PRINT_MESSAGE("Please see edit the file \"" << sRPFileName.Data() << "\" and add the missing data");
    exit(0);
  }
}

void GenerateParameterSub()
{
  ofstream *os;
  
  PRINT_MESSAGE("Parameter file not found.");
  PRINT_MESSAGE("Generating a stub parameter file \"" << sRPFileName.Data() << "\"");
  PRINT_MESSAGE("Please edit it and fill all the fileds.");
  os = new ofstream(sRPFileName.Data());
  
  (*os) 
    << "# This is a stub parameter file" << endl 
    << "# Fill the values and run the program again" << endl << endl;
  (*os) 
    << "# This is the list of available correlation fitters" << endl
 << "# CFFitter1DHBT" << endl
 << "# CFFitterNonId" << endl
 << "# CFFitterNonIdMult" << endl
 << "# CFFitterQuad" << endl
 << "# CFFitterSH" << endl
    << "CorrelationFitter = " << endl;
  (*os) 
    << "# This is the list of available source models" << endl
 << "# SourceModelCMSHyper" << endl
 << "# SourceModelCMSHyperBlast" << endl
 << "# SourceModelCMSHyperBt" << endl
 << "# SourceModelCMSHyperExpRQMD" << endl
 << "# SourceModelCMSHyperRQMD" << endl
 << "# SourceModelEHRInvShift" << endl
 << "# SourceModelExpGausRInv" << endl
 << "# SourceModelGaus" << endl
 << "# SourceModelGausCMS" << endl
 << "# SourceModelGausLCMS" << endl
 << "# SourceModelGausRInv" << endl
 << "# SourceModelGausRInvJacobian" << endl
 << "# SourceModelGausRInvShift" << endl
 << "# SourceModelGausROut" << endl
 << "# SourceModelGausROutDouble" << endl
 << "# SourceModelGausROutFile" << endl
 << "# SourceModelGausROutSourceCMS" << endl
 << "# SourceModelGausROutTwo" << endl
 << "# SourceModelGausRStar" << endl
    << "SourceModel = " << endl;
  (*os) 
    << "# This is the list of available weight generators" << endl
    << "# FsiLednicky - an FSI calculator from Richard Lednicky" << endl
    << "# FsiPratt    - an FSI calculator from Scott Pratt" << endl
    << "# FsiKisiel   - a simple FSI calculator for Coulomb and QS" << endl
    << "FsiCalculator = " << endl;
  (*os) 
    << "# The name of the output file, which will contain all the results" << endl
    << "OutResultFileNamePrefix = " << endl;
}

int main(int argc, char **argv)
{
  CFFitter    *tFitter;
  SourceModel *tModel;

  try {
    if (argc > 1) 
      sRPFileName = argv[1];
    else
      sRPFileName = "par.file";
    sRPInstance = new ReadPar(sRPFileName.Data());
    if (argc > 2)
      sRPInstance->parseOptString(argv[2]);
  }
  catch (int e) {
    PRINT_DEBUG("Caught: " << e);
    PRINT_MESSAGE("No Par File. Generating the stub...");
    GenerateParameterSub();
    exit(3);
  }

  ReadParameters();
  
  sOutFile = new TFile(sOutFileName.Data(),"recreate");
  
if (sSourceModelName == "SourceModelCMSHyper") { tModel = (SourceModelCMSHyper *) new SourceModelCMSHyper(); }
if (sSourceModelName == "SourceModelCMSHyperBlast") { tModel = (SourceModelCMSHyperBlast *) new SourceModelCMSHyperBlast(); }
if (sSourceModelName == "SourceModelCMSHyperBt") { tModel = (SourceModelCMSHyperBt *) new SourceModelCMSHyperBt(); }
if (sSourceModelName == "SourceModelCMSHyperExpRQMD") { tModel = (SourceModelCMSHyperExpRQMD *) new SourceModelCMSHyperExpRQMD(); }
if (sSourceModelName == "SourceModelCMSHyperRQMD") { tModel = (SourceModelCMSHyperRQMD *) new SourceModelCMSHyperRQMD(); }
if (sSourceModelName == "SourceModelEHRInvShift") { tModel = (SourceModelEHRInvShift *) new SourceModelEHRInvShift(); }
if (sSourceModelName == "SourceModelExpGausRInv") { tModel = (SourceModelExpGausRInv *) new SourceModelExpGausRInv(); }
if (sSourceModelName == "SourceModelGaus") { tModel = (SourceModelGaus *) new SourceModelGaus(); }
if (sSourceModelName == "SourceModelGausCMS") { tModel = (SourceModelGausCMS *) new SourceModelGausCMS(); }
if (sSourceModelName == "SourceModelGausLCMS") { tModel = (SourceModelGausLCMS *) new SourceModelGausLCMS(); }
if (sSourceModelName == "SourceModelGausRInv") { tModel = (SourceModelGausRInv *) new SourceModelGausRInv(); }
if (sSourceModelName == "SourceModelGausRInvJacobian") { tModel = (SourceModelGausRInvJacobian *) new SourceModelGausRInvJacobian(); }
if (sSourceModelName == "SourceModelGausRInvShift") { tModel = (SourceModelGausRInvShift *) new SourceModelGausRInvShift(); }
if (sSourceModelName == "SourceModelGausROut") { tModel = (SourceModelGausROut *) new SourceModelGausROut(); }
if (sSourceModelName == "SourceModelGausROutDouble") { tModel = (SourceModelGausROutDouble *) new SourceModelGausROutDouble(); }
if (sSourceModelName == "SourceModelGausROutFile") { tModel = (SourceModelGausROutFile *) new SourceModelGausROutFile(); }
if (sSourceModelName == "SourceModelGausROutSourceCMS") { tModel = (SourceModelGausROutSourceCMS *) new SourceModelGausROutSourceCMS(); }
if (sSourceModelName == "SourceModelGausROutTwo") { tModel = (SourceModelGausROutTwo *) new SourceModelGausROutTwo(); }
if (sSourceModelName == "SourceModelGausRStar") { tModel = (SourceModelGausRStar *) new SourceModelGausRStar(); }
  CFFitter::SetSourceModel(tModel);
if (sCFFitterName == "CFFitter1DHBT") { tFitter = (CFFitter1DHBT *) new CFFitter1DHBT(); }
if (sCFFitterName == "CFFitterNonId") { tFitter = (CFFitterNonId *) new CFFitterNonId(); }
if (sCFFitterName == "CFFitterNonIdMult") { tFitter = (CFFitterNonIdMult *) new CFFitterNonIdMult(); }
if (sCFFitterName == "CFFitterQuad") { tFitter = (CFFitterQuad *) new CFFitterQuad(); }
if (sCFFitterName == "CFFitterSH") { tFitter = (CFFitterSH *) new CFFitterSH(); }
  tFitter->Initialize();
  tFitter->Fit();
  sOutFile->cd();

  tFitter->Write();

  sOutFile->Close();
  
  return 0;
};


