void plotshrel()
{
  gStyle->SetOptStat(1000000001);
  TIter gdiriter(gDirectory->GetListOfKeys());
  TObject *cur;
  char *buf;
  char  suff[500];
  char  lbuf[500];
  char  nbuf[500];
  suff[0] = '\0';

  while (cur = gdiriter()) {
    buf = cur->GetName();
    if (strstr(buf, "CfnReYlm00")) {
      cout << "Found suffix " << buf+10 << endl;
      strcpy(suff, buf+10);
      break;
    }
  }

  if (!suff[0]) {
    cout << "Did not find suffix" << endl;
    return;
  }

  int maxl = 0;
  gdiriter->Reset();
  while (cur = gdiriter()) {
    buf = cur->GetName();
    if ((strstr(buf, "CfnReYlm")) && (strstr(buf, suff))){
      cout << "Found l " << buf+8 << endl;
      strncpy(lbuf, buf+8,1);
      lbuf[1] = '\0';
      int lcur = atoi(lbuf);
      cout << "Found l " << lcur << endl;
      if (lcur > maxl) maxl = lcur;
    }
  }
  cout << "Maximum l " << maxl << endl;
  

  int ilmcount= 0;
  TH1D **cfnsr;
  cfnsr = malloc(sizeof(TH1D *) * (maxl+1)*(maxl+1));
  TH1D **cfnsi;
  cfnsi = malloc(sizeof(TH1D *) * (maxl+1)*(maxl+1));
  for (int il=0; il<=maxl; il++)
    for (int im=0; im<=il; im++) {
      sprintf(nbuf, "CfnReYlm%i%i%s", il, im, suff);
      cfnsr[ilmcount] = new TH1D(*((TH1D *) gDirectory->Get(nbuf)));

      sprintf(nbuf, "CfnImYlm%i%i%s", il, im, suff);
      cfnsi[ilmcount] = new TH1D(*((TH1D *) gDirectory->Get(nbuf)));

      ilmcount++;
    }

  TLine *l1 = new TLine(0.0, 1.0, cfnsr[0]->GetXaxis()->GetXmax(), 1.0);
  l1->SetLineColor(14);
  l1->SetLineStyle(2);
  
  TLine *l0 = new TLine(0.0, 0.0, cfnsr[0]->GetXaxis()->GetXmax(), 0.0);
  l0->SetLineColor(14);
  l0->SetLineStyle(2);

  
  
  TCanvas *canshrel = new TCanvas("canshrel","canshrel",600,1000);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  canshrel->Divide(1,3,0.0001,0.0001);

  canshrel->cd(1);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  
  cfnsr[0]->SetTitle("");
  cfnsr[0]->SetMarkerColor(2);
  cfnsr[0]->SetMarkerSize(0.8);
  cfnsr[0]->SetMarkerStyle(8);
  cfnsr[0]->SetMaximum(1.8);
  cfnsr[0]->SetMinimum(0.9);
  cfnsr[0]->Draw();
  l1->Draw();

  canshrel->cd(2);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  
  cfnsr[3]->SetTitle("");
  cfnsr[3]->SetMarkerColor(2);
  cfnsr[3]->SetMarkerSize(0.8);
  cfnsr[3]->SetMarkerStyle(8);
  cfnsr[3]->SetMaximum(0.11);
  cfnsr[3]->SetMinimum(-0.11);
  cfnsr[3]->Draw();
  l0->Draw();

  canshrel->cd(3);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  
  cfnsr[5]->SetTitle("");
  cfnsr[5]->SetMarkerColor(2);
  cfnsr[5]->SetMarkerSize(0.8);
  cfnsr[5]->SetMarkerStyle(8);
  cfnsr[5]->SetMaximum(0.11);
  cfnsr[5]->SetMinimum(-0.11);
  cfnsr[5]->Draw();
  l0->Draw();

}
