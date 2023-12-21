#ifndef NAGAO_MACRO_H_
#define NAGAO_MACRO_H_

void SetTH1(TH1 *h, TString name, TString xname, TString yname, Int_t LColor=1, Int_t FStyle=0, Int_t FColor=0);
void SetTH2(TH2 *h, TString name, TString xname, TString yname, Double_t min=0.8);
void SetGr(TGraph *gr, TString name, TString xname, TString yname, Int_t LColor=1, Int_t LWidth=1, Int_t LStyle=1, Int_t MColor=1, Int_t MStyle=20, Double_t MSize=1.);
void SetTF1(TF1 *f, Int_t LColor=2, Int_t LWidth=1, Int_t LStyle=1, Int_t Npx=1000);

#endif  // NAGAO_MACRO_H_