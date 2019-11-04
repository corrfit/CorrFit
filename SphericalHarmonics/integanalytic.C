TF1 *funexcc00phi;
TF1 *funaprc00phi;

TF1 *funexcc22phi;
TF1 *funaprc22phi;
TF1 *funaprc22phib;

TGraph *grpexcc00phi;

TGraph *grpexcc22phi;

TGraph *getgrpfromfun(TF1 *fun, Double_t rmin, Double_t rmax) {
  Double_t xy[100];
  Double_t yk[100];
  
  for (int iter=0; iter<100; iter++) {
    xy[iter] = rmin + ((rmax-rmin)*iter)/100.0 + (rmax-rmin)/200.0;
    yk[iter] = fun->Eval(xy[iter]);

    //    cout << iter << " " << xy[iter] << " " << yk[iter] << endl;
  }

  TGraph * gra = new TGraph(100,xy,yk);

  return gra;
}

void integanalytic()
{
  double Rout  = 1.2/0.197327;
  double Rside = 1.0/0.197327;
  double Rlong = 1.5/0.197327;
  double qinv  = 0.25;

  funexcc00phi = new TF1("funexcc00phi","exp(-[0]*([1]*[1]*cos(x)*cos(x) + [2]*[2]*sin(x)*sin(x)))");
  funexcc00phi->SetRange(0.0,TMath::Pi()*2);
  
  funexcc00phi->SetParameter(0,qinv*qinv);
  funexcc00phi->SetParameter(1,Rout);
  funexcc00phi->SetParameter(2,Rside);

  grpexcc00phi = getgrpfromfun(funexcc00phi, 0.0, TMath::Pi()*2);

  funaprc00phi = new TF1("funaprc00phi","[0]+[1]*cos([2]*x)");
  funaprc00phi->SetParameters(1.0,-0.1,1.0);
  funaprc00phi->SetLineStyle(2);
  funaprc00phi->SetLineColor(2);
  funaprc00phi->SetRange(0.0,TMath::Pi()*2);

  //  grpexcc00phi->Fit(funaprc00phi);

  Double_t kk = qinv*qinv;
  Double_t av = 0.5*(exp(-kk*Rout*Rout) + exp(-kk*Rside*Rside));
  Double_t bv = 0.5*(exp(-kk*Rout*Rout) - exp(-kk*Rside*Rside));

  funaprc00phi->SetParameters(av, bv, 2.0);

  funexcc22phi = new TF1("funexcc22phi","exp(-[0]*([1]*[1]*cos(x)*cos(x) + [2]*[2]*sin(x)*sin(x)))*cos(2*x)");
  funexcc22phi->SetRange(0.0,TMath::Pi()*2);
  
  funexcc22phi->SetParameter(0,qinv*qinv);
  funexcc22phi->SetParameter(1,Rout);
  funexcc22phi->SetParameter(2,Rside);

  grpexcc22phi = getgrpfromfun(funexcc22phi, 0.0, TMath::Pi()*2);

  funaprc22phi = new TF1("funaprc22phi","[0]+[1]*cos([2]*x)");
  funaprc22phi->SetParameters(1.0,-0.1,1.0);
  funaprc22phi->SetLineStyle(2);
  funaprc22phi->SetLineColor(2);
  funaprc22phi->SetRange(0.0,TMath::Pi()*2);

  funaprc22phib = new TF1("funaprc22phib","[0]+[1]*exp(-[2]*[2]*x*x)+[1]*exp(-[2]*[2]*(x-[3])*(x-[3]))+[1]*exp(-[2]*[2]*(x-[4])*(x-[4]))+[1]*exp(-[2]*[2]*(x-[5])*(x-[5]))+[1]*exp(-[2]*[2]*(x-[6])*(x-[6]))");
  funaprc22phib->SetParameters(-0.2,0.3,1.2,TMath::Pi(),2*TMath::Pi());
  funaprc22phib->SetLineStyle(2);
  funaprc22phib->SetLineColor(2);
  funaprc22phib->FixParameter(3,TMath::Pi());
  funaprc22phib->FixParameter(4,2*TMath::Pi());
  funaprc22phib->FixParameter(5,3*TMath::Pi());
  funaprc22phib->FixParameter(6,-TMath::Pi());
  funaprc22phib->SetRange(0.0,TMath::Pi()*2);

  grpexcc22phi->Fit(funaprc22phib,"","",0.0,2*TMath::Pi());

//   Double_t kk = qinv*qinv;
//   Double_t av = 0.5*(exp(-kk*Rout*Rout) + exp(-kk*Rside*Rside));
//   Double_t bv = 0.5*(exp(-kk*Rout*Rout) - exp(-kk*Rside*Rside));

//   funaprc22phi->SetParameters(av, bv, 2.0);
}
