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
#include <TH2.h>
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

static const TString pdf_name =  OUTPUT_DIR + "/img/trig_explain.pdf";

std::unordered_map<std::string, std::vector<FitResult>> analyze(Int_t run_num, Int_t start_or_end = 0) // 0: mid_page, 1:start_page, 2: end_page, 3: both
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

    // -- prepare container -----
    std::unordered_map<std::string, std::vector<FitResult>> result_container;

    // +----------------+
    // | load root file |
    // +----------------+
    TString root_file_path = Form("%s/kekar_run%05d.root", DATA_DIR.Data(), run_num);
    auto *f = new TFile( root_file_path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return {};
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

    std::vector<HistPair> h_baca;
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        TString name  = Form("BACa_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d BAC(ADC) ch%d;ADC;", run_num, ch + 1);
        h_baca.emplace_back(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);
    }
    HistPair h_onsum_adc(Form("BAConSUMa_%d", run_num), Form("run%05d BAC online sum(ADC);ADC;", run_num), conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);
    HistPair h_offsum_adc(Form("BACoffSUMa_%d", run_num), Form("run%05d BAC offline sum(ADC);ADC;", run_num), conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);
    auto *h_bacsumt = new TH1D(Form("BACSUMt_%d", run_num), Form("run%05d BACSUM(TDC);TDC;", run_num), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    
    // -- NPE -----
    std::vector<HistPair> h_npe;
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        TString name  = Form("BACnpe_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d BAC(NPE) ch%d;NPE;", run_num, ch + 1);
        h_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);
    }
    HistPair h_onsum_npe(Form("BAConsumnpe_%d", run_num), Form("run%05d BAC online SUM(NPE);NPE;", run_num), conf.npe_bin_num, conf.npe_min, conf.npe_max);
    HistPair h_offsum_npe(Form("BACoffsumnpe_%d", run_num), Form("run%05d BAC offline sum(NPE);NPE;", run_num), conf.npe_bin_num, conf.npe_min, conf.npe_max);
   
    // shower event
    TH1D *h_npe_shower[conf.max_bac_ch];
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        TString name  = Form("BACnpe_shower_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d BAC(NPE) ch%d;NPE;", run_num, ch + 1);
        h_npe_shower[ch] = new TH1D(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);
    }
    auto *h_onsum_npe_shower = new TH1D(Form("BAConsumnpe_shower_%d", run_num), Form("run%05d BAC online SUM(NPE);NPE;", run_num), conf.npe_bin_num, conf.npe_min, conf.npe_max);
    auto *h_offsum_npe_shower = new TH1D(Form("BACoffsumnpe_shower_%d", run_num), Form("run%05d BAC offline SUM(NPE);NPE;", run_num), conf.npe_bin_num, conf.npe_min, conf.npe_max);

    // -- 2d histogram -----
    auto *h_correlation = new TH2D("on_off_correlation", ";online sum[adc];offline sum[npe]", conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max, conf.npe_bin_num, conf.npe_min, conf.npe_max);

    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    reader.Restart();
    while (reader.Next()){
        h_t1a->Fill(t1a[0]);
        h_t2a->Fill(t2a[0]);
        h_t3a->Fill(t3a[0]);
        h_t4a->Fill(t4a[0]);
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            h_t1t->Fill(t1t[n_hit]);
            h_t2t->Fill(t2t[n_hit]);
            h_t3t->Fill(t3t[n_hit]);
            h_t4t->Fill(t4t[n_hit]);
            h_bacsumt->Fill(bacsumt[n_hit]);
        }
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- make canvas and draw -----
    // -- create window -----
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    TGTab *tab = new TGTab(main, 1000, 800);


    TCanvas *c_trig1 = ana_helper::add_tab(tab, "trigger");
    c_trig1->Divide(2, 2);

    // -- T1 -----
    FitResult tmp_fit_result;
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t1t, c_trig1, 1);
    Double_t t1_tdc_min = tmp_fit_result.additional[0];
    Double_t t1_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t1a, c_trig1, 2);
    Double_t t1_adc_min = tmp_fit_result.additional[0];

    // -- T2 -----
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t2t, c_trig1, 3);
    Double_t t2_tdc_min = tmp_fit_result.additional[0];
    Double_t t2_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t2a, c_trig1, 4);
    Double_t t2_adc_min = tmp_fit_result.additional[0];

    
    TCanvas *c_trig2 = ana_helper::add_tab(tab, "trigger");
    c_trig2->Divide(2, 2);

    // -- T3 -----
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t3t, c_trig2, 1);
    Double_t t3_tdc_min = tmp_fit_result.additional[0];
    Double_t t3_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t3a, c_trig2, 2);
    Double_t t3_adc_min = tmp_fit_result.additional[0];
    Double_t t3_adc_max = tmp_fit_result.additional[1];
    
    // -- T4 -----
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t4t, c_trig2, 3);
    Double_t t4_tdc_min = tmp_fit_result.additional[0];
    Double_t t4_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t4a, c_trig2, 4);
    Double_t t4_adc_min = tmp_fit_result.additional[0];    
    Double_t t4_adc_max = tmp_fit_result.additional[1];    


    TCanvas *c_bact = ana_helper::add_tab(tab, "bac");
    
    // -- bac sum -----
    tmp_fit_result = ana_helper::cherenkov_tdc_fit(h_bacsumt, c_bact, 1);
    Double_t bac_tdc_min = tmp_fit_result.additional[0];
    Double_t bac_tdc_max = tmp_fit_result.additional[1];


    
    // +--------------------------+
    // | Preparation for 2nd fill |
    // +--------------------------+
    // -- load pedestal data -----
    TString pedestal_file_path = WORK_DIR + "/data/pedestal.root";
    auto *f_ped = new TFile( pedestal_file_path.Data() );
    if (!f_ped || f_ped->IsZombie()) {
        std::cerr << "Error: Could not open file : " << pedestal_file_path << std::endl;
        return {};
    }
    TTreeReader reader_ped("tree", f_ped);
    TTreeReaderValue<Int_t> run_num_ped(reader_ped, "run_num");
    TTreeReaderArray<Double_t> bac_ped_pos_val(reader_ped, "bac_ped_pos_val");
    TTreeReaderArray<Double_t> bacsum_ped_pos_val(reader_ped, "bacsum_ped_pos_val");
    TTreeReaderArray<Double_t> bac_ped_sig_val(reader_ped, "bac_ped_sig_val");
    TTreeReaderArray<Double_t> bacsum_ped_sig_val(reader_ped, "bacsum_ped_sig_val");
    reader_ped.Restart();
    while (reader_ped.Next()){
        if (*run_num_ped == ana_helper::get_pedestal_run_num(run_num)) break;
    }

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
    std::vector<Int_t> should_hit_ch = ana_helper::get_should_hit_ch(run_num);

    // +------------------+
    // | Fill event (2nd) |
    // +------------------+
    Int_t n_trig = 0, n_hitbac = 0;
    std::unordered_map<std::string, std::vector<std::pair<Bool_t, Double_t>>> online_sum_container{ {"raw", {}}, {"trig", {}} };
    
    reader.Restart();
    while (reader.Next()){
        // -- set up flag -----
        Bool_t do_hit_t1a = false, do_hit_t2a = false, do_hit_t3a = false, do_hit_t4a = false;
        Bool_t do_hit_t1t = false, do_hit_t2t = false, do_hit_t3t = false, do_hit_t4t = false, do_hit_bac = false;
        if ( t1_adc_min < t1a[0] ) do_hit_t1a = true;
        if ( t2_adc_min < t2a[0] ) do_hit_t2a = true;
        if ( t3_adc_min < t3a[0] ) do_hit_t3a = true;
        if ( t4_adc_min < t4a[0] ) do_hit_t4a = true;
        // if ( t3_adc_min < t3a[0] && t3a[0] < t3_adc_max ) do_hit_t3a = true;
        // if ( t4_adc_min < t4a[0] && t4a[0] < t4_adc_max ) do_hit_t4a = true;
        
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            if ( t1_tdc_min < t1t[n_hit] && t1t[n_hit] < t1_tdc_max ) do_hit_t1t = true;
            if ( t2_tdc_min < t2t[n_hit] && t2t[n_hit] < t2_tdc_max ) do_hit_t2t = true;
            if ( t3_tdc_min < t3t[n_hit] && t3t[n_hit] < t3_tdc_max ) do_hit_t3t = true;
            if ( t4_tdc_min < t4t[n_hit] && t4t[n_hit] < t4_tdc_max ) do_hit_t4t = true;
            if ( bac_tdc_min < bacsumt[n_hit] && bacsumt[n_hit] < bac_tdc_max ) do_hit_bac = true;
        }
        Bool_t trig_flag_adc = do_hit_t4a && do_hit_t2a && do_hit_t3a && do_hit_t4a;
        Bool_t trig_flag_tdc = do_hit_t4t && do_hit_t2t && do_hit_t3t && do_hit_t4t;

        // -- check shower -----
        Bool_t shower_flag = false;
        // for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        //     if (!std::binary_search(should_hit_ch.begin(), should_hit_ch.end(), ch) 
        //         && shower_adc_min[ch] < kvcsuma[ch] 
        //         && shower_tdc_gate[ch].first < kvcsumt[ch]
        //         && kvcsumt[ch] < shower_tdc_gate[ch].second 
        //     ) shower_flag = true;
        // }
        if ( t3_adc_max < t3a[0] || t4_adc_max < t4a[0] ) shower_flag = true;


        // -- event selection and fill data -----
        if ( trig_flag_tdc && trig_flag_adc ) {
            n_trig++;
            Double_t offline_sum_a   = 0.;
            Double_t offline_sum_npe = 0.;
            for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
                h_baca[ch].raw->Fill(baca[ch]);
                h_npe[ch].raw->Fill( (baca[ch] - bac_ped_pos_val[ch]) / conf.bac_opg[58][ch].first );
                offline_sum_a   +=  baca[ch] - bac_ped_pos_val[ch];
                offline_sum_npe += (baca[ch] - bac_ped_pos_val[ch]) / conf.bac_opg[58][ch].first;
            }
            h_offsum_adc.raw->Fill(offline_sum_a);
            h_onsum_adc.raw->Fill(bacsuma[0] - bacsum_ped_pos_val[0]);
            h_offsum_npe.raw->Fill(offline_sum_npe);
            online_sum_container["raw"].emplace_back(shower_flag, bacsuma[0]);
            
            if (do_hit_bac) {
                n_hitbac++;
                offline_sum_a   = 0.;
                offline_sum_npe = 0.;
                for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
                    h_baca[ch].trig->Fill( baca[ch] );
                    h_npe[ch].trig->Fill( (baca[ch] - bac_ped_pos_val[ch]) / conf.bac_opg[58][ch].first );
                    if (shower_flag) h_npe_shower[ch]->Fill( (baca[ch] - bac_ped_pos_val[ch]) / conf.bac_opg[58][ch].first );
                    offline_sum_a   +=  baca[ch] - bac_ped_pos_val[ch];
                    offline_sum_npe += (baca[ch] - bac_ped_pos_val[ch]) / conf.bac_opg[58][ch].first;
                }
                h_offsum_adc.trig->Fill(offline_sum_a);
                h_onsum_adc.trig->Fill(bacsuma[0] - bacsum_ped_pos_val[0]);
                h_offsum_npe.trig->Fill(offline_sum_npe);
                online_sum_container["trig"].emplace_back(shower_flag, bacsuma[0]);

                // 最初にOnline sumのone photon gainをoffline sumとの相関より推測する
                h_correlation->Fill( bacsuma[0] - bacsum_ped_pos_val[0], offline_sum_npe );

                if (shower_flag) h_offsum_npe_shower->Fill(offline_sum_npe);
            }
        }
    }

    FitResult eff_result;
    eff_result.additional.push_back( static_cast<Double_t>(n_trig) );
    eff_result.additional.push_back( static_cast<Double_t>(n_hitbac) );
    result_container["eff"].push_back(eff_result);

    // +----------------------------------------+
    // | Estimate one photon gain of online sum |
    // +----------------------------------------+
    // -- fitting and estimate -----
    TCanvas *c_corr = ana_helper::add_tab(tab, "corr");

    FitResult linear_fit_result = ana_helper::correlation_fit(h_correlation, c_corr, 1);
    result_container["linear"].push_back(linear_fit_result);

    // -- calc and fill -----
    for (const auto& pair : online_sum_container["raw"]) h_onsum_npe.raw->Fill( linear_fit_result.par[0]*(pair.second - bacsum_ped_pos_val[0]) + linear_fit_result.par[1] );
    for (const auto& pair : online_sum_container["trig"]) {
        h_onsum_npe.trig->Fill( linear_fit_result.par[0]*(pair.second - bacsum_ped_pos_val[0]) + linear_fit_result.par[1] );
        if (pair.first) h_onsum_npe_shower->Fill( linear_fit_result.par[0]*(pair.second - bacsum_ped_pos_val[0]) + linear_fit_result.par[1] );
    }


    //  +-----------------+
    //  | Draw histograms |
    //  +-----------------+
    // -- indiv adc -----
    TCanvas *c_indiva = ana_helper::add_tab(tab, "adc");
    c_indiva->Divide(2, 2);

    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        c_indiva->cd(ch+1);
        // gPad->SetLogy(1);
        Double_t peak_pos = h_baca[ch].raw->GetMean();
        Double_t stdev = h_baca[ch].raw->GetStdDev();
        // h_baca[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_baca[ch].raw->SetLineColor(kBlue);
        // h_baca[ch].raw->Draw();
        // h_baca[ch].trig->SetLineColor(kRed);
        // h_baca[ch].trig->SetFillColor(kRed);
        // h_baca[ch].trig->SetFillStyle(3003);
        h_baca[ch].trig->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_baca[ch].trig->Draw("same");
    }

    // -- indiv NPE -----
    TCanvas *c_indivnpe = ana_helper::add_tab(tab, "adc");
    c_indivnpe->Divide(2, 2);

    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        c_indivnpe->cd(ch+1);
        // gPad->SetLogy(1);
        // h_npe[ch].raw->GetXaxis()->SetRangeUser(-5, 100.0);
        // h_npe[ch].raw->SetLineColor(kBlue);
        // h_npe[ch].raw->Draw();
        // h_npe[ch].trig->SetLineColor(kRed);
        // h_npe[ch].trig->SetFillColor(kRed);
        // h_npe[ch].trig->SetFillStyle(3003);
        FitResult result = ana_helper::npe_gauss_fit(h_npe[ch].trig, c_indivnpe, ch+1);
        result_container["indiv_npe"].push_back(result);
        h_npe_shower[ch]->SetLineColor(kGreen);
        h_npe_shower[ch]->SetFillColor(kGreen);
        h_npe_shower[ch]->SetFillStyle(3003);
        h_npe_shower[ch]->Draw("same");
    }

    // -- sum adc -----
    TCanvas *c_suma = ana_helper::add_tab(tab, "adc");
    c_suma->Divide(2, 1);
    {
        c_suma->cd(1);
        Double_t peak_pos = h_offsum_adc.raw->GetMean();
        Double_t stdev = h_offsum_adc.raw->GetStdDev();
        // h_offsum_adc.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_offsum_adc.raw->SetLineColor(kBlue);
        // h_offsum_adc.raw->Draw();
        // h_offsum_adc.trig->SetLineColor(kRed);
        // h_offsum_adc.trig->SetFillColor(kRed);
        // h_offsum_adc.trig->SetFillStyle(3003);
        h_offsum_adc.trig->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_offsum_adc.trig->Draw("same");

        c_suma->cd(2);
        // h_onsum_adc.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_onsum_adc.raw->SetLineColor(kBlue);
        // h_onsum_adc.raw->Draw();
        // h_onsum_adc.trig->SetLineColor(kRed);
        // h_onsum_adc.trig->SetFillColor(kRed);
        // h_onsum_adc.trig->SetFillStyle(3003);
        h_onsum_adc.trig->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_onsum_adc.trig->Draw("same");
    }

    // -- sum npe -----
    TCanvas *c_sumnpe = ana_helper::add_tab(tab, "adc");
    c_sumnpe->Divide(2, 1);
    {        
        c_sumnpe->cd(1);
        // Double_t peak_pos = h_offsum_npe.raw->GetMean();
        // Double_t stdev = h_offsum_npe.raw->GetStdDev();
        // h_offsum_npe.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_offsum_npe.raw->SetLineColor(kBlue);
        // h_offsum_npe.raw->Draw();
        // h_offsum_npe.trig->SetLineColor(kRed);
        // h_offsum_npe.trig->SetFillColor(kRed);
        // h_offsum_npe.trig->SetFillStyle(3003);
        FitResult offsum_result = ana_helper::npe_gauss_fit(h_offsum_npe.trig, c_sumnpe, 1, 1.5);
        result_container["offsum_npe"].push_back(offsum_result);
        h_offsum_npe_shower->SetLineColor(kGreen);
        h_offsum_npe_shower->SetFillColor(kGreen);
        h_offsum_npe_shower->SetFillStyle(3003);
        h_offsum_npe_shower->Draw("same");

        c_sumnpe->cd(2);
        // h_onsum_npe.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_onsum_npe.raw->SetLineColor(kBlue);
        // h_onsum_npe.raw->Draw();
        // h_onsum_npe.trig->SetLineColor(kRed);
        // h_onsum_npe.trig->SetFillColor(kRed);
        // h_onsum_npe.trig->SetFillStyle(3003);
        FitResult onsum_result = ana_helper::npe_gauss_fit(h_onsum_npe.trig, c_sumnpe, 2, 1.5);
        result_container["onsum_npe"].push_back(onsum_result);
        h_onsum_npe_shower->SetLineColor(kGreen);
        h_onsum_npe_shower->SetFillColor(kGreen);
        h_onsum_npe_shower->SetFillStyle(3003);
        h_onsum_npe_shower->Draw("same");
    }

    // -- add tab and draw window -----
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();
    

    return result_container;
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    conf.bac_initialize();

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
    analyze(run_num, 3);
    theApp->Run();

    return 0;
}