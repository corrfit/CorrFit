// This macro reads two files:
// a file with the correctly calculated correlation function C
// a file with the correctly calculated purity correction histogram P
// both histograms need to be in spherical harmonics
// Then it proceeds to subtract 1.0 from C
// and save it in a new file as numerator
// and save the P histogram as denominator

void PreparePurityCorrectionFiles(const char *inCfile, const char *inPfile, const char *inCname, const char *inPname, int lmax= 3)
{
  TFile *infC = new TFile(inCfile);
  TFile *infP = new TFile(inPfile);
  TFile *outC = new TFile("shmout.forPC.root","RECREATE");

  TH1D *numsre[16];
  TH1D *densre[16];
  TH1D *numsim[16];
  TH1D *densim[16];
  
  TH3D *covnum;
  TH3D *covden;

  int ih=0;
  char bufname[200];
  char outCname[200];
  sprintf(outCname, "%sPC", inCname);

  for (int il=0; il<=lmax; il++) {
    for (int im=0; im<=il*2; im++) {
      sprintf(bufname, "CfnReYlm%i%i%s", il, im, inCname);
      numsre[ih] = new TH1D (*((TH1D *) infC->Get(bufname)));
      sprintf(bufname, "NumReYlm%i%i%s", il, im, outCname);
      numsre[ih]->SetName(bufname);
      numsre[ih]->SetTitle(bufname);
      
      if ((il == 0) && (im ==0)) {
	for (int ib=1; ib <=numsre[ih]->GetNbinsX(); ib++)
	  numsre[ih]->SetBinContent(ib, numsre[ih]->GetBinContent(ib)-1.0);
	numsre[ih]->SetEntries(numsre[ih]->GetNbinsX());
      }

      sprintf(bufname, "CfnImYlm%i%i%s", il, im, inCname);
      numsim[ih] = new TH1D (*((TH1D *) infC->Get(bufname)));
      sprintf(bufname, "NumImYlm%i%i%s", il, im, outCname);
      numsim[ih]->SetName(bufname);
      numsim[ih]->SetTitle(bufname);
      
      sprintf(bufname, "CfnReYlm%i%i%s", il, im, inPname);
      densre[ih] = new TH1D (*((TH1D *) infP->Get(bufname)));
      sprintf(bufname, "DenReYlm%i%i%s", il, im, outCname);
      densre[ih]->SetName(bufname);
      densre[ih]->SetTitle(bufname);
      
      sprintf(bufname, "CfnImYlm%i%i%s", il, im, inPname);
      densim[ih] = new TH1D (*((TH1D *) infP->Get(bufname)));
      sprintf(bufname, "DenImYlm%i%i%s", il, im, outCname);
      densim[ih]->SetName(bufname);
      densim[ih]->SetTitle(bufname);
      
      outC->cd();
      numsre[ih]->Write();
      numsim[ih]->Write();
      densre[ih]->Write();
      densim[ih]->Write();

      ih++;
    }
  }

  sprintf(bufname, "CovCfc%s", inCname);
  covnum = new TH3D(*((TH3D *) infC->Get(bufname)));
  sprintf(bufname, "CovNum%s", outCname);
  covnum->SetName(bufname);
  covnum->SetTitle(bufname);
  
  sprintf(bufname, "CovCfc%s", inPname);
  covden = new TH3D(*((TH3D *) infP->Get(bufname)));
  sprintf(bufname, "CovDen%s", outCname);
  covden->SetName(bufname);
  covden->SetTitle(bufname);
  
  outC->cd();
  covnum->Write();
  covden->Write();

}
