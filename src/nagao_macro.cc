#include <iostream>

#include "TMath.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TColor.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TGaxis.h"

const Double_t PI = 4.*atan(1.);
const Double_t deg_to_rad = PI / 180.;
const Double_t rad_to_deg = 180. / PI;
const Double_t sigma_to_fwhm = 2.*sqrt(2.*log(2.));
const Double_t fwhm_to_sigma = 1./sigma_to_fwhm;
const Double_t cm_to_barn = 1e+24;
const Double_t alpha = 1./137.035999074; // fine structure constant
const Double_t hbar = 6.58211928*1e-22;  // Planc constant (reduced) (MeV x s)
const Double_t hbarc = 197.3269718;      // conversion constant (MeV x fm)
const Double_t kb = 8.6173324*1e-5;      // Boltzmann constant
const Double_t e = 1.602176565*1e-19;    // electron charge magnitude (C)
const Double_t c = 0.299792458;          // speed of light in vacuum (m/ns)
const Double_t re = 2.817e-13;           // classical electron radius (cm)
const Double_t Na = 6.02214129*1e+23;    // Avogadro constant
const Double_t Me = 0.510998928;         // electron     mass (MeV/c2)
const Double_t Mmu = 105.6583715;        // muon         mass (MeV/c2)
const Double_t Mpi = 139.57018;          // charged pion mass (MeV/c2)
const Double_t Mpi0 = 134.9766;          // charged pion mass (MeV/c2)
const Double_t MK = 493.677;             // charged Kaon mass (MeV/c2)
const Double_t Mp = 938.272046;          // proton       mass (MeV/c2)
const Double_t Mn = 939.565379;          // proton       mass (MeV/c2)
const Double_t Mu = 931.494061;          // proton       mass (MeV/c2)
const Double_t ML = 1115.683;            // Lambda       mass (MeV/c2)
const Double_t MS0 = 1192.642;           // Sigma Zero   mass (MeV/c2)
const Double_t MSm = 1197.449;           // Sigma Minus  mass (MeV/c2)
const Double_t MSp = 1189.37;            // Sigma Plus   mass (MeV/c2)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH1(TH1 *h, TString name, TString xname, TString yname, Int_t LColor=1, Int_t FStyle=0, Int_t FColor=0){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetTitleSize(0.04,"");
  h->SetTitleFont(42,"");
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(5);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH2(TH2 *h, TString name, TString xname, TString yname, Double_t min=0.8){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(1);
  h->SetTitleSize(0.05,"");
  // h->SetMarkerStyle(20);
  // h->SetMarkerSize(1.5);
  // h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(5);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetGr(TGraph *gr, TString name, TString xname, TString yname, Int_t LColor=1, Int_t LWidth=1, Int_t LStyle=1, Int_t MColor=1, Int_t MStyle=20, Double_t MSize=1.){
  gr->SetTitle(name);
  gr->SetName(name);

  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetTitleOffset(0.90);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelOffset(0.01);

  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleOffset(1.00);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)gr->GetYaxis())->SetMaxDigits(4);

  gr->SetLineColor(LColor);
  gr->SetLineStyle(LStyle);
  gr->SetLineWidth(LWidth);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTF1(TF1 *f, Int_t LColor=2, Int_t LWidth=1, Int_t LStyle=1, Int_t Npx=1000){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}
