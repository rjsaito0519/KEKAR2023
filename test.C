#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>

#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TColor.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TTreeReader.h"
#include "TParticle.h"
#include "TLatex.h"
#include "TH2Poly.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include "include/macro.hh"
#include "include/nagao_macro.hh"
#include "include/variable.hh"

void analyze(TString path)
{
    // いろんな設定
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.05, "XY");
    gStyle->SetTitleSize(1, "XY");
    gStyle->SetTitleFontSize(0.08);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // -- load file ------------------------------------------------------------------
    auto *f = new TFile(path.Data());

    // -- prepare reader and histogram -----------------------------------------------------------------------------------
    TTreeReader reader("tree", f);
    TTreeReaderArray<Double_t> kvca(reader, "kvca");
    TTreeReaderArray<Double_t> kvcsuma(reader, "kvcsuma");

    // -- prepare histogram ------------------------------------------------------------------
    TH1D *hKVCa[maxKVCch];
    for (Int_t ch = 0; ch < maxKVCch; ch++) hKVCa[ch] = new TH1D(Form("kvca_%d", ch), "", adc_bin_num, adc_min, adc_max);
    TH1D *hKVCSUMa[maxKVCSUMch];
    for (Int_t ch = 0; ch < maxKVCSUMch; ch++) hKVCSUMa[ch] = new TH1D(Form("kvcsuma_%d", ch), "", adc_bin_num, adc_min, adc_max);

    // -- fill ------------------------------------------------------------------
    reader.Restart();
    while (reader.Next()){
        for (Int_t ch = 0; ch < maxKVCch; ch++) {
            hKVCa[ch]->Fill( kvca[ch] );
        }
        for (Int_t ch = 0; ch < maxKVCSUMch; ch++) {
            hKVCSUMa[ch]->Fill( kvcsuma[ch] );
        }
    }

    TCanvas *c1 = new TCanvas("", "", 1500, 900);
    c1->Divide(4, 2);
    for (Int_t ch = 0; ch < maxKVCch; ch++) {
        c1->cd(ch+1);
        hKVCa[ch]->GetXaxis()->SetRangeUser(0, 350);
        hKVCa[ch]->Draw();
    }
    for (Int_t ch = 0; ch < maxKVCSUMch; ch++) {
        c1->cd(ch+5);
        hKVCSUMa[ch]->GetXaxis()->SetRangeUser(0, 350);
        hKVCSUMa[ch]->Draw();
    }
    
    std::cout << "finish" << std::endl;

}

Int_t main(int argc, char** argv) {
    TString path = argv[1];
    TApplication *theApp = new TApplication("App", &argc, argv);    
    analyze(path);
    theApp->Run();
    return 0;
}
