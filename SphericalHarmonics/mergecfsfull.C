void mergecfsfull(const char *in1, const char *in2, const char *out, const char *pair="PipPp")
{
  TFile *inf1 = new TFile(in1);
  TFile *inf2 = new TFile(in2);
  TFile *outf = new TFile(out, "RECREATE");

  char fname[200];

  sprintf(fname, "NumReYlm%i%i%sNonIdYlms", 0, 0, pair);
  cout << "Reading " << fname << endl;
      
  TH1D *nm1 = inf1->Get(fname);
  TH1D *nm2 = inf2->Get(fname);

  sprintf(fname, "CovCfc%sNonIdYlms", pair);
  TH3D *cv1 = inf1->Get(fname);
  TH3D *cv2 = inf2->Get(fname);

  TH3D *cvo = new TH3D(*cv1);

  TH1D *cfs1r[32];
  TH1D *cfs2r[32];
  
  int lcount = 0;

  for (int el=0; el<=3; el++) {
    for (int em=0; em<el*2+1; em++) {

      sprintf(fname, "CfnReYlm%i%i%sNonIdYlms", el, em, pair);
      cout << "Reading " << fname << endl;
      
      TH1D *cf1 = inf1->Get(fname);
      TH1D *cf2 = inf2->Get(fname);

      cfs1r[lcount] = cf1;
      cfs2r[lcount] = cf2;
      
      TH1D *cfo = new TH1D(*cf1);
      for (int iter=1; iter<=cfo->GetNbinsX(); iter++) {
	double w1 = nm1->GetBinContent(iter);
	if (w1>0.0) {
	  w1 *= w1;
	  w1 = 1/w1;
	}
	else 
	  w1 = 0;
	double w2 = nm2->GetBinContent(iter);
	if (w2 > 0.0) {
	  w2 *= w2;
	  w2 = 1/w2;
	}
	else
	  w2 = 0;

	if ((w1>0.0) && (w2>0.0)) {
	  double sum = w1+w2;
	  w1/=sum;
	  w2/=sum;

	  double mean1 = cf1->GetBinContent(iter);
	  double mean2 = cf2->GetBinContent(iter);

	  double sigm1 = cf1->GetBinError(iter);
	  double sigm2 = cf2->GetBinError(iter);

	  cfo->SetBinContent(iter, mean1*w1 + mean2*w2);
	  cfo->SetBinError(iter, TMath::Sqrt(w1*sigm1*sigm1 + w2*sigm2*sigm2 + w1*w2*(mean1-mean2)*(mean1-mean2)));
	}
	else {
	  cfo->SetBinContent(iter,0.0);
	  cfo->SetBinError(iter, 0.0);
	}
      }
      outf->cd();
      cfo->Write();
      
      lcount++;
      
      sprintf(fname, "CfnImYlm%i%i%sNonIdYlms", el, em, pair);
      cout << "Reading " << fname << endl;
      
      cf1 = (TH1D *) inf1->Get(fname);
      cf2 = (TH1D *) inf2->Get(fname);

      cfs1r[lcount] = cf1;
      cfs2r[lcount] = cf2;
      
      cfo = new TH1D(*cf1);

      for (int iter=1; iter<=cfo->GetNbinsX(); iter++) {
	double w1 = nm1->GetBinContent(iter);
	if (w1>0.0) {
	  w1 *= w1;
	  w1 = 1/w1;
	}
	else 
	  w1 = 0;
	double w2 = nm2->GetBinContent(iter);
	if (w2 > 0.0) {
	  w2 *= w2;
	  w2 = 1/w2;
	}
	else
	  w2 = 0;

	if ((w1>0.0) && (w2>0.0)) {
	  double sum = w1+w2;
	  w1/=sum;
	  w2/=sum;

	  double mean1 = cf1->GetBinContent(iter);
	  double mean2 = cf2->GetBinContent(iter);

	  double sigm1 = cf1->GetBinError(iter);
	  double sigm2 = cf2->GetBinError(iter);

	  cfo->SetBinContent(iter, mean1*w1 + mean2*w2);
	  cfo->SetBinError(iter, TMath::Sqrt(w1*sigm1*sigm1 + w2*sigm2*sigm2 + w1*w2*(mean1-mean2)*(mean1-mean2)));
	}
	else {
	  cfo->SetBinContent(iter,0.0);
	  cfo->SetBinError(iter, 0.0);
	}
      }
      outf->cd();
      cfo->Write();
      
      lcount++;
    }
  }
  
  for (int ibin=1; ibin<cfo->GetNbinsX(); ibin++) {
    double w1 = nm1->GetBinContent(ibin);
    double w2 = nm2->GetBinContent(ibin);
    for (int ilm1=0; ilm1<32; ilm1++)
      for (int ilm2=0; ilm2<32; ilm2++) {
	double c1 = cv1->GetBinContent(ibin, ilm1+1, ilm2+1);
	double c2 = cv2->GetBinContent(ibin, ilm1+1, ilm2+1);

	double xm1 = cfs1r[ilm1]->GetBinContent(ibin);
	double ym1 = cfs1r[ilm2]->GetBinContent(ibin);

	double xm2 = cfs2r[ilm1]->GetBinContent(ibin);
	double ym2 = cfs2r[ilm2]->GetBinContent(ibin);

	cvo->SetBinContent(ibin, ilm1+1, ilm2+1, 
			   w1*c1 + w2*c2 + w1*w2*(xm1-xm2)*(ym1-ym2));
      }
  }

  cvo->SetBinContent(0,0,0,cv1->GetBinContent(0,0,0));

  outf->cd();
  cvo->Write();
}
