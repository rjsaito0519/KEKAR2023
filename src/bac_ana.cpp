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
    
    HistPair h_bacsuma(Form("BACSUMa_%d", run_num), Form("run%05d BACSUM(ADC);ADC;", run_num), conf.adc_bin_num, conf.adc_min, conf.adc_max);

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
               
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            h_t1t->Fill(t1t[n_hit]);
            h_t2t->Fill(t2t[n_hit]);
            h_t3t->Fill(t3t[n_hit]);
            h_t4t->Fill(t4t[n_hit]);
            h_bacsumt->Fill(bacsumt[n_hit]);
        }

        h_bacsuma.raw->Fill(bacsuma[0]);
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
    FitResult tmp_fit_result;

    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t1t, c, 1);
    Double_t t1_tdc_min = tmp_fit_result.additional[0];
    Double_t t1_tdc_max = tmp_fit_result.additional[1];

    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t2t, c, 2);
    Double_t t2_tdc_min = tmp_fit_result.additional[0];
    Double_t t2_tdc_max = tmp_fit_result.additional[1];

    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t3t, c, 3);
    Double_t t3_tdc_min = tmp_fit_result.additional[0];
    Double_t t3_tdc_max = tmp_fit_result.additional[1];

    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t4t, c, 4);
    Double_t t4_tdc_min = tmp_fit_result.additional[0];
    Double_t t4_tdc_max = tmp_fit_result.additional[1];

    tmp_fit_result = ana_helper::cherenkov_tdc_fit(h_bacsumt, c, 5);
    Double_t bac_tdc_min = tmp_fit_result.additional[0];
    Double_t bac_tdc_max = tmp_fit_result.additional[1];

    tmp_fit_result = ana_helper::trig_counter_adc_fit(h_bacsuma.raw, c, 6);
    
    auto *c3 = new TCanvas("", "", 800, 800);
    c3->cd(1);
    h_bacsuma.raw->Draw();


    TCanvas *c2 = ana_helper::add_tab(tab, "adc");
    c2->Divide(2, 2);
    tmp_fit_result = ana_helper::trig_counter_adc_fit(h_t1a, c2, 1);
    Double_t t1_adc_min = tmp_fit_result.additional[0];

    tmp_fit_result = ana_helper::trig_counter_adc_fit(h_t2a, c2, 2);
    Double_t t2_adc_min = tmp_fit_result.additional[0];

    tmp_fit_result = ana_helper::trig_counter_adc_fit(h_t3a, c2, 3);
    Double_t t3_adc_min = tmp_fit_result.additional[0];

    tmp_fit_result = ana_helper::trig_counter_adc_fit(h_t4a, c2, 4);
    Double_t t4_adc_min = tmp_fit_result.additional[0];

        
    
    // -- prepare kvc veto flag -----
    std::vector<Double_t> shower_adc_min;
    std::vector<std::pair<Double_t, Double_t>> shower_tdc_gate;
    if ( run_num <= conf.separate_run_num ) {
        shower_adc_min  = conf.adc_min_condition1;
        shower_tdc_gate = conf.tdc_gate_condition1;
    } else {
        shower_adc_min  = conf.adc_min_condition2;
        shower_tdc_gate = conf.tdc_gate_condition2;    
    }
    // Bool_t veto






    // // +------------------+
    // // | Fill event (2nd) |
    // // +------------------+
    // Int_t n_trig = 0, n_sachit = 0;
    // std::vector<Double_t> online_sum_container{}, trig_online_sum_container{}; // for online sum npe
    // Double_t sum_pedestal_pos = sacsum_pedestal_pos[ pedestal_index(run_num) ]; // sumのpedestalは一番近い別のpedestal runからとってくる
    
    
    // reader.Restart();
    // while (reader.Next()){
    //     // -- set up flag -----
    //     Bool_t do_hit_t1a = false, do_hit_t2a = false, do_hit_t3a = false, do_hit_t4a = false;
    //     Bool_t do_hit_t1t = false, do_hit_t2t = false, do_hit_t3t = false, do_hit_t4t = false, do_hit_bac = false;
    //     if ( t1_adc_min < t1a[0] ) do_hit_t1a = true;
    //     if ( t2_adc_min < t2a[0] ) do_hit_t2a = true;
    //     if ( t3_adc_min < t3a[0] ) do_hit_t3a = true;
    //     if ( t4_adc_min < t4a[0] ) do_hit_t4a = true;
    //     for (Int_t n_hit = 0; n_hit < max_nhit_tdc; n_hit++) {
    //         if ( t1_tdc_min < t1t[n_hit] && t1t[n_hit] < t1_tdc_max ) do_hit_t1t = true;
    //         if ( t2_tdc_min < t2t[n_hit] && t2t[n_hit] < t2_tdc_max ) do_hit_t2t = true;
    //         if ( t3_tdc_min < t3t[n_hit] && t3t[n_hit] < t3_tdc_max ) do_hit_t3t = true;
    //         if ( t4_tdc_min < t4t[n_hit] && t4t[n_hit] < t4_tdc_max ) do_hit_t4t = true
    //         if ( bac_tdc_min < bacsumt[n_hit] && bacsumt[n_hit] < bac_tdc_max ) do_hit_bac = true;
    //     }
    //     Bool_t trig_flag_adc = do_hit_t4a && do_hit_t2a && do_hit_t3a && do_hit_t4a;
    //     Bool_t trig_flag_tdc = do_hit_t4t && do_hit_t2t && do_hit_t3t && do_hit_t4t;

    //     // -- event selection and fill data -----
    //     if ( trig_flag_tdc && trig_flag_adc ) {
    //         n_trig++;
    //         online_sum_container.push_back(sacsuma[0]);
    //         if (do_hit_bac) {
    //             n_sachit++;
    //             Double_t offline_sum_a   = 0.;
    //             Double_t offline_sum_npe = 0.;
    //             for (Int_t ch = 0; ch < max_sac_ch; ch++) {
    //                 h_trig_saca[ch]->Fill( saca[ch] );
    //                 h_npe[ch]->Fill( (saca[ch] - result_pedestal[ch][1]) / sac_one_photon_gain[ch][2] );
    //                 offline_sum_a   +=  saca[ch] - result_pedestal[ch][1];
    //                 offline_sum_npe += (saca[ch] - result_pedestal[ch][1]) / sac_one_photon_gain[ch][2];
    //             }
    //             h_offline_sum_a->Fill(offline_sum_a);
    //             h_online_sum_a ->Fill(sacsuma[0] - sum_pedestal_pos);
    //             h_offline_sum_npe->Fill(offline_sum_npe);
    //             trig_online_sum_container.push_back(sacsuma[0]);

    //             // 最初にOnline sumのone photon gainをoffline sumとの相関より推測する
    //             h_correlation->Fill( offline_sum_npe, sacsuma[0] - sum_pedestal_pos );
    //         }
    //     }
    // }









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

    // -- add tab and draw window -----
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();

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