#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <random>

// ROOT libraries
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TGTab.h>
#include <TColor.h>

// Custom headers
#include "config.h"
#include "param.h"
#include "ana_helper.h"
#include "paths.h"


void analyze(Int_t run_num, Int_t ch, TVirtualPad *c, Int_t n_c)
{   
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0/256, 160.0/256, 44.0/256);
    
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x"); // x軸のタイトルサイズ
    gStyle->SetTitleSize(0.06, "y"); // y軸のタイトルサイズ
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // -- parameter -----
    Config& conf = Config::getInstance();

    // +----------------+
    // | load root file |
    // +----------------+
    TString root_file_path = Form("%s/kekar_run%05d.root", DATA_DIR.Data(), run_num);
    auto *f = new TFile( root_file_path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return;
    }
    TTreeReader reader("tree", f);
    TTreeReaderArray<Double_t> baca(reader, "kvca");

    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    auto *h = new TH1D(Form("BACa_%d_%d", run_num, ch+1), Form("run%05d BAC(ADC) ch%d;ADC;", run_num, ch+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);

    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    reader.Restart();
    while (reader.Next()){
        h->Fill(baca[ch]);
    }

    Double_t peak_pos = h->GetMean();
    Double_t stdev = h->GetStdDev();
    c->cd(n_c);
    h->GetXaxis()->SetRangeUser(peak_pos - 3.0*stdev, peak_pos + 5.0*stdev);
    h->Draw();

}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

    // +-------------+
    // | dev version |
    // +-------------+
    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
        return 1;
    }
    Int_t run_num = std::atoi(argv[1]);

    TApplication *theApp = new TApplication("App", &argc, argv);

    // -- create window -----
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    TGTab *tab = new TGTab(main, 1000, 800);

    // -- test -----
    TCanvas *c = ana_helper::add_tab(tab, "kvc");
    c->Divide(2, 2);
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) analyze(run_num, ch, c, ch+1);

    // -- add tab and draw window -----
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();
    
    theApp->Run();

    return 0;
}