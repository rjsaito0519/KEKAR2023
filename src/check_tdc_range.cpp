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

static TH1D *h_sacsumt;
static TH1D *h_bacsumt;
static TH1D *h_kvcsumt[4];


void analyze(Int_t run_num)
{   
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x"); // x軸のタイトルサイズ
    gStyle->SetTitleSize(0.06, "y"); // y軸のタイトルサイズ
    gStyle->SetTitleFontSize(0.06);
    gROOT->GetColor(0)->SetAlpha(0.01);

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
    TTreeReaderArray<Double_t> sacsumt(reader, "sacsumt");
    TTreeReaderArray<Double_t> bacsumt(reader, "bacsumt");
    TTreeReaderArray<Double_t> kvcsumt(reader, "kvcsumt");
    

    // +------------+
    // | Fill event |
    // +------------+
    reader.Restart();
    while (reader.Next()){
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            h_sacsumt->Fill(sacsumt[n_hit]);
            h_bacsumt->Fill(bacsumt[n_hit]);
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) h_kvcsumt[ch]->Fill(kvcsumt[conf.max_nhit_tdc*ch+n_hit]);
        }
    }

    std::cout << run_num << std::endl;

    delete f;
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    TApplication *theApp = new TApplication("App", &argc, argv);

    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    h_sacsumt = new TH1D("SACSUMt", "SACSUM(TDC);TDC;", conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    h_bacsumt = new TH1D("BACSUMt", "BACSUM(TDC);TDC;", conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        h_kvcsumt[ch] = new TH1D(Form("KVCSUMt_%d", ch+1), Form("KVCSUM(TDC) ch%d;TDC;", ch+1), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    std::vector<Int_t> ana_run_num{
        // 269, 270, 271, 272, 273,
        // 274, 275, 276, 277, 278,
        // 280, 281, 282, 283, 284,

        412, 414, 415, 416, 417, 445,
        418, 419, 420, 421, 422,
        423, 424, 425, 426, 427

    };

    for (const auto& run_num : ana_run_num) analyze(run_num);

    // -- create window -----
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    TGTab *tab = new TGTab(main, 1000, 800);

    TCanvas *c_sac = ana_helper::add_tab(tab, "sac");
    c_sac->cd(1);
    h_sacsumt->Draw();

    TCanvas *c_bac = ana_helper::add_tab(tab, "bac");
    c_bac->cd(1);
    h_bacsumt->Draw();

    TCanvas *c_kvc = ana_helper::add_tab(tab, "kvc");
    c_kvc->Divide(2, 2);
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        c_kvc->cd(ch+1);
        h_kvcsumt[ch]->Draw();
    }

    // -- add tab and draw window -----
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();

    theApp->Run();

    return 0;
}