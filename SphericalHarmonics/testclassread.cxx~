#include <stdlib.h>
#include <iostream>
#include "CorrFctnDirectYlm.h"
#include <TH1D.h>
#include <TRandom2.h>
#include <math.h>
#include <TFile.h>
#include <complex>
#include <TObject.h>
#include <TKey.h>

using namespace std;

int main(int argc, char **argv)
{
  TFile *ofile = new TFile("shmout.recalc.root","RECREATE");
  TFile *infile = new TFile(argv[1]);
  CorrFctnDirectYlm *cylm;

  int maxl = -1;
  int maxlset;

  if (argc>3)
    maxlset = atoi(argv[3]);
  else
    maxlset = -1;
  cout << "Maximum L is " << maxlset << endl;
  char  suff[500];
  if (1) {
    TH1D *bufh;

    TIter gdiriter(infile->GetListOfKeys());
    TObject *cur;
    char  buf[1000];
    char  lbuf[500];
    char  nbuf[500];
    suff[0] = '\0';
    
    while (cur = gdiriter()) {
      //      cout << "buf is " << cur->GetName() << endl;
      strcpy(buf, cur->GetName());
      if (strstr(buf, "NumReYlm00")) {
	cout << "Found suffix " << buf+10 << endl;
	strcpy(suff, buf+10);
	bufh = (TH1D *) ((TKey *) cur)->ReadObj();
	if (argc > 2)
	  if (strcmp(suff, argv[2]))
	    cout << "Suffix does not match parameter" << endl;
	  else
	    break;
	else
	  break;
      }
      if (strstr(buf, "CfnReYlm00")) {
	cout << "Found suffix " << buf+10 << endl;
	strcpy(suff, buf+10);
	bufh = (TH1D *) ((TKey *) cur)->ReadObj();
	if (argc > 2)
	  if (strcmp(suff, argv[2]))
	    cout << "Suffix does not match parameter" << endl;
	  else
	    break;
	else
	  break;
      }
    }
    
    if (!suff[0]) {
      cout << "Did not find suffix" << endl;
      return -1;
    }
    
    gdiriter.Reset();
    while (cur = gdiriter()) {
      strcpy(buf, cur->GetName());
      if ((strstr(buf, "NumReYlm")) && (strstr(buf, suff))){
	cout << "Found l " << buf+8 << endl;
	strncpy(lbuf, buf+8,1);
	lbuf[1] = '\0';
	int lcur = atoi(lbuf);
	cout << "Found l " << lcur << endl;
	if (lcur > maxl) maxl = lcur;
      }
      if ((strstr(buf, "CfnReYlm")) && (strstr(buf, suff))){
	cout << "Found l " << buf+8 << endl;
	strncpy(lbuf, buf+8,1);
	lbuf[1] = '\0';
	int lcur = atoi(lbuf);
	cout << "Found l " << lcur << endl;
	if (lcur > maxl) maxl = lcur;
      }
    }
    
    cout << "Found max l L " << maxl << " " << maxlset << endl;

    if (maxlset>=0)
      if ((maxlset<maxl) || (maxl < 0)) maxl = maxlset;
    cout << "Maximum l " << maxl << endl;
    cout << "Hist has bins min max " << bufh->GetNbinsX() << " " << bufh->GetXaxis()->GetXmin() << " " << bufh->GetXaxis()->GetXmax() << endl;
    

    cylm = new CorrFctnDirectYlm(argv[2],maxl,bufh->GetNbinsX(),bufh->GetXaxis()->GetXmin(),bufh->GetXaxis()->GetXmax());
  }
  else {
    cylm = new CorrFctnDirectYlm(argv[2],maxl,60,0.0,0.1);
    strcpy (suff, argv[2]);
  }

#ifdef _SH_EXP_
  if (argc>2) {
    Double_t radius = atof(argv[3])/0.197327;
    Double_t bohrac = atof(argv[4])/0.197327;
    
    cylm->SetAdvancedNormalization(radius/0.197327, bohrac/0.197327, 0.72, 21,40);
  }
#endif
  cylm->ReadFromFile(infile, suff, maxl);

  ofile->cd();
  cylm->Write();

  delete cylm;
}

  
