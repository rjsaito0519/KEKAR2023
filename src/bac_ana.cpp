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
#include "one_photon_gain_helper.h"
#include "ana_helper.h"

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

    // -- parameter -----
    Config& conf = Config::getInstance();

    // +----------------+
    // | load root file |
    // +----------------+
    TString root_file_path = Form( "../root/kekar_run%05d.root", run_num );
    auto *f = new TFile( root_file_path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return;
    }
    TTreeReader reader("tree", f);
    TTreeReaderArray<Double_t> t1a(reader, "t1a");
    TTreeReaderArray<Double_t> t1t(reader, "t1t");
    TTreeReaderArray<Double_t> t2a(reader, "t2a");
    TTreeReaderArray<Double_t> t2t(reader, "t2t");
    TTreeReaderArray<Double_t> t3a(reader, "t3a");
    TTreeReaderArray<Double_t> t3t(reader, "t3t");
    TTreeReaderArray<Double_t> t4a(reader, "t4a");
    TTreeReaderArray<Double_t> t4t(reader, "t4t");
    TTreeReaderArray<Double_t> baca(reader, "baca");
    TTreeReaderArray<Double_t> bacsuma(reader, "bacsuma");
    TTreeReaderArray<Double_t> bacsumt(reader, "bacsumt");
    TTreeReaderArray<Double_t> kvcsuma(reader, "kvcsuma");
    TTreeReaderArray<Double_t> kvcsumt(reader, "kvcsumt");
    

    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    auto *h_t1a = new TH1D(Form("T1a_%d", run_num), Form("run%05d T1(ADC);ADC;", run_num), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t1t = new TH1D(Form("T1t_%d", run_num), Form("run%05d T1(TDC);TDC;", run_num), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    auto *h_t2a = new TH1D(Form("T2a_%d", run_num), Form("run%05d T2(ADC);ADC;", run_num), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t2t = new TH1D(Form("T2t_%d", run_num), Form("run%05d T2(TDC);TDC;", run_num), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    auto *h_t3a = new TH1D(Form("T3a_%d", run_num), Form("run%05d T3(ADC);ADC;", run_num), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t3t = new TH1D(Form("T3t_%d", run_num), Form("run%05d T3(TDC);TDC;", run_num), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    auto *h_t4a = new TH1D(Form("T4a_%d", run_num), Form("run%05d T4(ADC);ADC;", run_num), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t4t = new TH1D(Form("T4t_%d", run_num), Form("run%05d T4(TDC);TDC;", run_num), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);

    TH1D *h_baca[conf.max_bac_ch];
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) h_baca[ch] = new TH1D(Form("BACa_%d_%d", run_num, ch+1), Form("run%05d BAC(ADC) ch%d;ADC;", run_num, ch+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_bacsumt = new TH1D(Form("BACSUMt_%d", run_num), Form("run%05d BACSUM(TDC);TDC;", run_num), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);

    TH1D *h_kvcsuma[conf.max_kvc_ch];
    TH1D *h_kvcsumt[conf.max_kvc_ch];
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        h_kvcsuma[ch] = new TH1D(Form("KVCSUMa_%d_%d", run_num, ch+1), Form("run%05d KVCSUM(ADC) ch%d;ADC;", run_num, ch+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);
        h_kvcsumt[ch] = new TH1D(Form("KVCSUMt_%d_%d", run_num, ch+1), Form("run%05d KVCSUM(TDC) ch%d;TDC;", run_num, ch+1), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    reader.Restart();
    while (reader.Next()){
        h_t1a->Fill(t1a[0]);
        h_t2a->Fill(t2a[0]);
        h_t3a->Fill(t3a[0]);
        h_t4a->Fill(t4a[0]);
        for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) h_baca[ch]->Fill(baca[ch]);
        for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) h_kvcsuma[ch]->Fill(kvcsuma[ch]);
               
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            h_t1t->Fill(t1t[n_hit]);
            h_t2t->Fill(t2t[n_hit]);
            h_t3t->Fill(t3t[n_hit]);
            h_t4t->Fill(t4t[n_hit]);
            h_bacsumt->Fill(bacsumt[n_hit]);
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) h_kvcsumt[ch]->Fill(kvcsumt[n_hit + ch*conf.max_kvc_ch]);
        }
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- create window -----
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    TGTab *tab = new TGTab(main, 1000, 800);

    // -- test -----
    TCanvas *c = ana_helper::add_tab(tab, "bac");
    c->Divide(3, 2);
    ana_helper::trig_counter_tdc_fit(h_t1t, c, 1);
    ana_helper::trig_counter_tdc_fit(h_t2t, c, 2);
    ana_helper::trig_counter_tdc_fit(h_t3t, c, 3);
    ana_helper::trig_counter_tdc_fit(h_t4t, c, 4);
    ana_helper::trig_counter_tdc_fit(h_bacsumt, c, 5);
    

    // -- add tab and draw window -----
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();


    // // -- prepare -----
    // c->cd(n_c);
    
    // std::vector<Double_t> result, result_err;
    // std::pair<Double_t, Double_t> mean{ h->GetMean(),  h->GetMeanError() };
    // Double_t n_total = h->GetEntries();

    // // -- load and set parameter -----
    // TString key = Form("%05d-%d", run_num, ch);
    // std::vector<Double_t> par = param::bac_opg.count(key.Data()) ? param::bac_opg.at(key.Data()) : param::bac_opg.at("default");
    // Double_t n_gauss         = par[0];
    // Double_t first_peak_pos  = par[1];
    // Double_t fit_range_left  = par[2];
    // Double_t fit_range_right = par[3];
    // h->GetXaxis()->SetRangeUser(fit_range_left-5.0, fit_range_right+5.0);
    // Double_t half_width = 5.0;

    // // -- first fit -----
    // TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", first_peak_pos-half_width, first_peak_pos+half_width);
    // f_prefit->SetParameter(1, first_peak_pos);
    // f_prefit->SetParameter(2, half_width);
    // h->Fit(f_prefit, "0Q", "", first_peak_pos-half_width, first_peak_pos+half_width );
    // for (Int_t i = 0; i < 3; i++) result.push_back(f_prefit->GetParameter(i));
    // delete f_prefit;

    // // -- second fit -----
    // TString func_str = "[0]*TMath::Gaus(x, [1], [2], true)";
    // for (Int_t i = 1; i < n_gauss; i++) func_str += Form(" + [%d]*TMath::Gaus(x, [1]+%d*[3], TMath::Sqrt(%d)*[4], true)", i+4, i, i+1);
    // TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), func_str.Data(), fit_range_left, fit_range_right);
    // f_fit->SetParameter(0, result[0]);
    // f_fit->SetParameter(1, result[1]);
    // f_fit->SetParameter(2, result[2]*0.9);
    // Double_t peak_to_peak = (fit_range_right-fit_range_left-2*result[2])/n_gauss;
    // f_fit->SetParameter(3, peak_to_peak);
    // f_fit->SetParameter(4, result[2]);
    // for (Int_t i = 1; i < n_gauss; i++) {
    //     h->GetBinCenter( h->GetMaximumBin() );
    //     Double_t height = h->GetBinContent(h->FindBin( result[1]+i*peak_to_peak ));
    //     f_fit->SetParameter(i+4, height*TMath::Sqrt(i+1)*result[2]);
    //     f_fit->SetParLimits(i+4, 0, 0x10000000);
    // }
    // f_fit->SetLineColor(kOrange);
    // f_fit->SetLineWidth(2);
    // h->Fit(f_fit, "0", "", fit_range_left, fit_range_right);
    // result.clear();
    // Int_t n_par = f_fit->GetNpar();
    // for (Int_t i = 0; i < n_par; i++) {
    //     result.push_back(f_fit->GetParameter(i));
    //     result_err.push_back(f_fit->GetParError(i));
    // }

    // // -- for adjust fit range -----
    // std::cout << result[3]*n_gauss + 2*result[2] + fit_range_left << std::endl;

    // // -- cal one photon gain -----
    // std::pair<Double_t, Double_t> pedestal{ f_fit->GetParameter(1), f_fit->GetParError(1) };
    // std::pair<Double_t, Double_t> n_pedestal{ f_fit->GetParameter(0), f_fit->GetParError(0) };
    // std::pair<Double_t, Double_t> one_photon_gain = cal_one_photon_gain( mean, pedestal, n_pedestal, n_total);
    // // std::cout << one_photon_gain.first << ", " << one_photon_gain.second << std::endl;
    // result.push_back(one_photon_gain.first);
    // result_err.push_back(one_photon_gain.second);

    // // -- draw -----
    // h->Draw();
    // f_fit->Draw("same");

    // TF1 *f_first_gaus = new TF1(Form("gauss_%s_0", h->GetName()), "gausn", fit_range_left, fit_range_right);
    // f_first_gaus->SetParameters(result[0], result[1], result[2]);
    // f_first_gaus->SetLineColor(kBlue);
    // f_first_gaus->SetLineStyle(2);
    // // f_first_gaus->SetFillColor(kBlue);
    // // f_first_gaus->SetFillStyle(3003);
    // f_first_gaus->Draw("same");

    // for (Int_t i = 1; i < n_gauss; i++) {
    //     TF1 *f_single_gaus = new TF1(Form("gauss_%s_%d", h->GetName(), i), "gausn", fit_range_left, fit_range_right);
    //     f_single_gaus->SetParameters(result[i+4], result[1]+i*result[3], TMath::Sqrt(i+1)*result[4]);
    //     f_single_gaus->SetLineColor(kBlue);
    //     f_single_gaus->SetLineStyle(2);
    //     f_single_gaus->Draw("same");    
    // }    
    // c->Update();

    // // // -- delete -----
    // // delete f;

    // std::vector<std::vector<Double_t>> result_container;
    // result_container.push_back(result);
    // result_container.push_back(result_err);

    // return result_container;
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

    analyze(run_num);
    
    theApp->Run();

    return 0;
}