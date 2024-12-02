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

static const TString pdf_name =  OUTPUT_DIR + "/img/kvc_2cm_pos_scan_analysis.pdf";

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
    TTreeReaderArray<Double_t> kvca(reader, "kvca");
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

    std::vector<HistPair> h_kvca;
    std::vector<HistPair> h_onsum_adc;
    std::vector<HistPair> h_offsum_adc;
    TH1D* h_kvcsumt[conf.max_kvc_ch];
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        TString name  = Form("KVCa_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d KVC(ADC) seg2 %s %d;ADC;", run_num, ch/2==0 ? "UP" : "DOWN", ch%2+1);
        h_kvca.emplace_back(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);

        name  = Form("KVConSUMa_%d_%d", run_num, ch+1);
        title = Form("run%05d KVC online sum(ADC) seg%d;ADC;", run_num, ch+1);
        h_onsum_adc.emplace_back(name, title, conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);

        name  = Form("KVCoffSUMa_%d_%d", run_num, ch+1);
        title = Form("run%05d KVC offline sum(ADC) seg%d;ADC;", run_num, ch+1);
        h_offsum_adc.emplace_back(name, title, conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);
    
        h_kvcsumt[ch] = new TH1D(Form("KVCSUMt_%d_%d", run_num, ch+1), Form("run%05d KVCSUM(TDC) seg%d;TDC;", run_num, ch+1), conf.adjust_tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }
    
    // -- NPE -----
    std::vector<HistPair> h_npe;
    std::vector<HistPair> h_onsum_npe;
    std::vector<HistPair> h_offsum_npe;
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        TString name  = Form("KVCnpe_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d KVC(NPE) seg2 %s %d;NPE;", run_num, ch/2==0 ? "UP" : "DOWN", ch%2+1);
        h_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);

        name  = Form("KVConsumnpe_%d_%d", run_num, ch+1);
        title = Form("run%05d KVC online SUM(NPE) seg%d;NPE;", run_num, ch+1);
        h_onsum_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);

        name  = Form("KVCoffsumnpe_%d_%d", run_num, ch+1);
        title = Form("run%05d KVC offline SUM(NPE) seg%d;NPE;", run_num, ch+1);
        h_offsum_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);
    }
   
    // shower event
    TH1D *h_npe_shower[conf.max_kvc_ch];
    TH1D *h_onsum_npe_shower[conf.max_kvc_ch];
    TH1D *h_offsum_npe_shower[conf.max_kvc_ch];
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        TString name  = Form("KVCnpe_shower_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d KVC(NPE) seg2 %s %d;NPE;", run_num, ch/2==0 ? "UP" : "DOWN", ch%2+1);
        h_npe_shower[ch] = new TH1D(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);

        name  = Form("KVConsumnpe_shower_%d_%d", run_num, ch+1);
        title =  Form("run%05d KVC online SUM(NPE) seg%d;NPE;", run_num, ch+1);
        h_onsum_npe_shower[ch] = new TH1D(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);

        name  = Form("KVCoffsumnpe_shower_%d_%d", run_num, ch+1);
        title =  Form("run%05d KVC offline SUM(NPE) seg%d;NPE;", run_num, ch+1);
        h_offsum_npe_shower[ch] = new TH1D(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);
    }

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
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
                h_kvcsumt[ch]->Fill(kvcsumt[conf.max_nhit_tdc*ch+n_hit]);
            }
        }
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- make canvas and draw -----
    Int_t nth_pad = 1;
    const Int_t rows = 2;
    const Int_t cols = 2;
    const Int_t max_pads = rows * cols;
    auto *c = new TCanvas("", "", 1500, 1200);
    c->Divide(cols, rows);
    if (start_or_end == 1 || start_or_end == 3) c->Print(pdf_name + "["); // start


    // -- T1 -----
    FitResult tmp_fit_result;
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t1t, c, nth_pad);
    Double_t t1_tdc_min = tmp_fit_result.additional[0];
    Double_t t1_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t1a, c, ++nth_pad);
    Double_t t1_adc_min = tmp_fit_result.additional[0];
    Double_t t1_adc_max = tmp_fit_result.additional[1];

    // -- T2 -----
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t2t, c, ++nth_pad);
    Double_t t2_tdc_min = tmp_fit_result.additional[0];
    Double_t t2_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t2a, c, ++nth_pad);
    Double_t t2_adc_min = tmp_fit_result.additional[0];
    Double_t t2_adc_max = tmp_fit_result.additional[1];

    c->Print(pdf_name);
    c->Clear();
    c->Divide(cols, rows);
    nth_pad = 1;

    // -- T3 -----
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t3t, c, nth_pad);
    Double_t t3_tdc_min = tmp_fit_result.additional[0];
    Double_t t3_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t3a, c, ++nth_pad);
    Double_t t3_adc_min = tmp_fit_result.additional[0];
    Double_t t3_adc_max = tmp_fit_result.additional[1];

    // -- T4 -----
    tmp_fit_result = ana_helper::trig_counter_tdc_fit(h_t4t, c, ++nth_pad);
    Double_t t4_tdc_min = tmp_fit_result.additional[0];
    Double_t t4_tdc_max = tmp_fit_result.additional[1];
    tmp_fit_result = ana_helper::trig_counter_adc_gauss_fit(h_t4a, c, ++nth_pad);
    Double_t t4_adc_min = tmp_fit_result.additional[0];
    Double_t t4_adc_max = tmp_fit_result.additional[1];
    
    c->Print(pdf_name);
    c->Clear();
    c->Divide(cols, rows);
    nth_pad = 1;


    // -- kvc sum -----
    Double_t kvc_tdc_min[conf.max_kvc_ch];
    Double_t kvc_tdc_max[conf.max_kvc_ch];
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        tmp_fit_result = ana_helper::cherenkov_tdc_fit(h_kvcsumt[ch], c, nth_pad);
        kvc_tdc_min[ch] = tmp_fit_result.additional[0];
        kvc_tdc_max[ch] = tmp_fit_result.additional[1];
        nth_pad++;
    }

    
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
    TTreeReaderArray<Double_t> kvc_ped_pos_val(reader_ped, "kvc_ped_pos_val");
    TTreeReaderArray<Double_t> kvcsum_ped_pos_val(reader_ped, "kvcsum_ped_pos_val");
    TTreeReaderArray<Double_t> kvc_ped_sig_val(reader_ped, "kvc_ped_sig_val");
    TTreeReaderArray<Double_t> kvcsum_ped_sig_val(reader_ped, "kvcsum_ped_sig_val");
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
    Int_t n_trig = 0, n_hitkvc = 0;
    std::unordered_map<std::string, std::vector<std::vector<std::pair<Bool_t, Double_t>>>> online_sum_container{ {"raw", {{}, {}, {}, {}}}, {"trig", {{}, {}, {}, {}}} };
    
    reader.Restart();
    while (reader.Next()){
        // -- set up flag -----
        Bool_t do_hit_t1a = false, do_hit_t2a = false, do_hit_t3a = false, do_hit_t4a = false;
        Bool_t do_hit_t1t = false, do_hit_t2t = false, do_hit_t3t = false, do_hit_t4t = false, do_hit_kvc = false;
        // if ( t1_adc_min < t1a[0] ) do_hit_t1a = true;
        // if ( t2_adc_min < t2a[0] ) do_hit_t2a = true;
        // if ( t3_adc_min < t3a[0] ) do_hit_t3a = true;
        // if ( t4_adc_min < t4a[0] ) do_hit_t4a = true;
        if ( t1_adc_min < t1a[0] && t1a[0] < t1_adc_max ) do_hit_t1a = true;
        if ( t2_adc_min < t2a[0] && t2a[0] < t2_adc_max ) do_hit_t2a = true;
        if ( t3_adc_min < t3a[0] && t3a[0] < t3_adc_max ) do_hit_t3a = true;
        if ( t4_adc_min < t4a[0] && t4a[0] < t4_adc_max ) do_hit_t4a = true;

        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            if ( t1_tdc_min < t1t[n_hit] && t1t[n_hit] < t1_tdc_max ) do_hit_t1t = true;
            if ( t2_tdc_min < t2t[n_hit] && t2t[n_hit] < t2_tdc_max ) do_hit_t2t = true;
            if ( t3_tdc_min < t3t[n_hit] && t3t[n_hit] < t3_tdc_max ) do_hit_t3t = true;
            if ( t4_tdc_min < t4t[n_hit] && t4t[n_hit] < t4_tdc_max ) do_hit_t4t = true;
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) if ( kvc_tdc_min[ch] < kvcsumt[conf.max_nhit_tdc*ch+n_hit] && kvcsumt[conf.max_nhit_tdc*ch+n_hit] < kvc_tdc_max[ch] ) do_hit_kvc = true;
        }
        Bool_t trig_flag_adc = do_hit_t4a && do_hit_t2a && do_hit_t3a && do_hit_t4a;
        Bool_t trig_flag_tdc = do_hit_t4t && do_hit_t2t && do_hit_t3t && do_hit_t4t;

        // -- check shower -----
        Bool_t shower_flag = false;
        // for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        //     if (!std::binary_search(should_hit_ch.begin(), should_hit_ch.end(), ch) && shower_adc_min[ch] < kvcsuma[ch] ) {
        //         for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
        //             if (shower_tdc_gate[ch].first < kvcsumt[conf.max_nhit_tdc*ch+n_hit] && kvcsumt[conf.max_nhit_tdc*ch+n_hit] < shower_tdc_gate[ch].second ) shower_flag = true;
        //         }
        //     }
        // }

        // -- event selection and fill data -----
        if ( trig_flag_tdc && trig_flag_adc ) {
            n_trig++;
            Double_t offline_sum_a   = 0.0;
            Double_t offline_sum_npe = 0.0;
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
                h_kvca[ch].raw->Fill(kvca[ch]);
                Int_t index_opg = conf.max_kvc_ch * ch + 1; // seg2
                h_npe[ch].raw->Fill( (kvca[ch] - kvc_ped_pos_val[ch]) / conf.kvc_thick_opg[56][index_opg].first );
                offline_sum_a   +=  kvca[ch] - kvc_ped_pos_val[ch];
                offline_sum_npe += (kvca[ch] - kvc_ped_pos_val[ch]) / conf.kvc_thick_opg[56][index_opg].first;

                h_onsum_adc[ch].raw->Fill(kvcsuma[ch] - kvcsum_ped_pos_val[ch]);
                online_sum_container["raw"][ch].emplace_back(shower_flag, kvcsuma[ch]);
            }
            h_offsum_adc[1].raw->Fill(offline_sum_a); // seg 2
            h_offsum_npe[1].raw->Fill(offline_sum_npe);
            
            if (do_hit_kvc) {
                n_hitkvc++;
                offline_sum_a   = 0.0;
                offline_sum_npe = 0.0;
                for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
                    h_kvca[ch].trig->Fill(kvca[ch]);
                    Int_t index_opg = conf.max_kvc_ch * ch + 1; // seg2
                    h_npe[ch].trig->Fill( (kvca[ch] - kvc_ped_pos_val[ch]) / conf.kvc_thick_opg[56][index_opg].first );
                    if (shower_flag) h_npe_shower[ch]->Fill( (kvca[ch] - kvc_ped_pos_val[ch]) / conf.kvc_thick_opg[56][index_opg].first );
                    offline_sum_a   +=  kvca[ch] - kvc_ped_pos_val[ch];
                    offline_sum_npe += (kvca[ch] - kvc_ped_pos_val[ch]) / conf.kvc_thick_opg[56][index_opg].first;

                    h_onsum_adc[ch].trig->Fill(kvcsuma[ch] - kvcsum_ped_pos_val[ch]);
                    online_sum_container["trig"][ch].emplace_back(shower_flag, kvcsuma[ch]);
                }
                h_offsum_adc[1].trig->Fill(offline_sum_a); // seg 2
                h_offsum_npe[1].trig->Fill(offline_sum_npe);

                // 最初にOnline sumのone photon gainをoffline sumとの相関より推測する
                h_correlation->Fill( kvcsuma[1] - kvcsum_ped_pos_val[1], offline_sum_npe );

                if (shower_flag) {
                    h_offsum_npe_shower[1]->Fill(offline_sum_npe);
                }
            }
        }
    }

    FitResult eff_result;
    eff_result.additional.push_back( static_cast<Double_t>(n_trig) );
    eff_result.additional.push_back( static_cast<Double_t>(n_hitkvc) );
    result_container["eff"].push_back(eff_result);

    // +----------------------------------------+
    // | Estimate one photon gain of online sum |
    // +----------------------------------------+
    // -- fitting and estimate -----
    c->Print(pdf_name);
    c->Clear();
    c->Divide(1, 1);
    nth_pad = 1;
    FitResult linear_fit_result = ana_helper::correlation_fit(h_correlation, c, nth_pad);
    result_container["linear"].push_back(linear_fit_result);

    // -- calc and fill -----
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        Double_t coeff_a, coeff_b;
        if (ch == 1) {
            coeff_a = linear_fit_result.par[0];
            coeff_b = linear_fit_result.par[1];
        } else {
            // coeff_a = linear_fit_result.par[0] 
            //           * (conf.kvc_thick_opg[56][1].first + conf.kvc_thick_opg[56][5].first + conf.kvc_thick_opg[56][9].first + conf.kvc_thick_opg[56][13].first)
            //           / (conf.kvc_thick_opg[56][ch].first + conf.kvc_thick_opg[56][ch+conf.max_kvc_ch].first + conf.kvc_thick_opg[56][ch+2*conf.max_kvc_ch].first + conf.kvc_thick_opg[56][ch+3*conf.max_kvc_ch].first);

            // coeff_b = linear_fit_result.par[1] 
            //           * (conf.kvc_thick_opg[56][1].first + conf.kvc_thick_opg[56][5].first + conf.kvc_thick_opg[56][9].first + conf.kvc_thick_opg[56][13].first)
            //           / (conf.kvc_thick_opg[56][ch].first + conf.kvc_thick_opg[56][ch+conf.max_kvc_ch].first + conf.kvc_thick_opg[56][ch+2*conf.max_kvc_ch].first + conf.kvc_thick_opg[56][ch+3*conf.max_kvc_ch].first);

            coeff_a = 4.0/(conf.kvc_thick_opg[56][ch].first + conf.kvc_thick_opg[56][ch+conf.max_kvc_ch].first + conf.kvc_thick_opg[56][ch+2*conf.max_kvc_ch].first + conf.kvc_thick_opg[56][ch+3*conf.max_kvc_ch].first);
            coeff_b = 0.0;

            // coeff_a = (1.0/conf.kvc_thick_opg[56][ch].first + 1.0/conf.kvc_thick_opg[56][ch+conf.max_kvc_ch].first + 1.0/conf.kvc_thick_opg[56][ch+2*conf.max_kvc_ch].first + 1.0/conf.kvc_thick_opg[56][ch+3*conf.max_kvc_ch].first)/4.0;
            // coeff_b = 0.0;

            // std::cout << linear_fit_result.par[0] << ", " << linear_fit_result.par[1] << std::endl;
            // std::cout << 4.0/(conf.kvc_thick_opg[56][1].first + conf.kvc_thick_opg[56][5].first + conf.kvc_thick_opg[56][9].first + conf.kvc_thick_opg[56][13].first) << std::endl;
            // std::cout << coeff_a << ", " << coeff_b << std::endl;
        }

        for (const auto& pair : online_sum_container["raw"][ch]) h_onsum_npe[ch].raw->Fill( coeff_a*(pair.second - kvcsum_ped_pos_val[ch]) + coeff_b );
        
        for (const auto& pair : online_sum_container["trig"][ch]) {
            h_onsum_npe[ch].trig->Fill( coeff_a*(pair.second - kvcsum_ped_pos_val[ch]) + coeff_b );
            if (pair.first) h_onsum_npe_shower[ch]->Fill( coeff_a*(pair.second - kvcsum_ped_pos_val[ch]) + coeff_b );
        }
    }


    //  +-----------------+
    //  | Draw histograms |
    //  +-----------------+
    // -- indiv adc -----
    c->Print(pdf_name);
    c->Clear();
    c->Divide(cols, rows);
    nth_pad = 1;

    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        c->cd(nth_pad);
        // gPad->SetLogy(1);
        Double_t peak_pos = h_kvca[ch].raw->GetMean();
        Double_t stdev = h_kvca[ch].raw->GetStdDev();
        // h_kvca[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_kvca[ch].raw->SetLineColor(kBlue);
        // h_kvca[ch].raw->Draw();
        // h_kvca[ch].trig->SetLineColor(kRed);
        // h_kvca[ch].trig->SetFillColor(kRed);
        // h_kvca[ch].trig->SetFillStyle(3003);
        h_kvca[ch].trig->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_kvca[ch].trig->Draw("same");
        nth_pad++;
    }

    // -- indiv NPE -----
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        c->cd(nth_pad);
        // gPad->SetLogy(1);
        // h_npe[ch].raw->GetXaxis()->SetRangeUser(-5, 100.0);
        // h_npe[ch].raw->SetLineColor(kBlue);
        // h_npe[ch].raw->Draw();
        // h_npe[ch].trig->SetLineColor(kRed);
        // h_npe[ch].trig->SetFillColor(kRed);
        // h_npe[ch].trig->SetFillStyle(3003);
        FitResult result = ana_helper::npe_gauss_fit(h_npe[ch].trig, c, nth_pad);
        result_container["indiv_npe"].push_back(result);
        h_npe_shower[ch]->SetLineColor(kGreen);
        h_npe_shower[ch]->SetFillColor(kGreen);
        h_npe_shower[ch]->SetFillStyle(3003);
        h_npe_shower[ch]->Draw("same");
        nth_pad++;
    }

    // -- sum adc -----
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        c->cd(nth_pad);
        Double_t peak_pos = h_offsum_adc[ch].raw->GetMean();
        Double_t stdev = h_offsum_adc[ch].raw->GetStdDev();
        // h_offsum_adc[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_offsum_adc[ch].raw->SetLineColor(kBlue);
        // h_offsum_adc[ch].raw->Draw();
        // h_offsum_adc[ch].trig->SetLineColor(kRed);
        // h_offsum_adc[ch].trig->SetFillColor(kRed);
        // h_offsum_adc[ch].trig->SetFillStyle(3003);
        h_offsum_adc[ch].trig->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_offsum_adc[ch].trig->Draw("same");
        nth_pad++;

        c->cd(nth_pad);
        // h_onsum_adc[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_onsum_adc[ch].raw->SetLineColor(kBlue);
        // h_onsum_adc[ch].raw->Draw();
        // h_onsum_adc[ch].trig->SetLineColor(kRed);
        // h_onsum_adc[ch].trig->SetFillColor(kRed);
        // h_onsum_adc[ch].trig->SetFillStyle(3003);
        h_onsum_adc[ch].trig->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_onsum_adc[ch].trig->Draw("same");
        nth_pad++;
    }

    // -- sum npe -----
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }

        Double_t cutoff_threshold = 0.0;
        if (std::binary_search(should_hit_ch.begin(), should_hit_ch.end(), ch)) cutoff_threshold = 50.0;

        c->cd(nth_pad);
        // Double_t peak_pos = h_offsum_npe[ch].raw->GetMean();
        // Double_t stdev = h_offsum_npe[ch].raw->GetStdDev();
        // h_offsum_npe[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_offsum_npe[ch].raw->SetLineColor(kBlue);
        // h_offsum_npe[ch].raw->Draw();
        // h_offsum_npe[ch].trig->SetLineColor(kRed);
        // h_offsum_npe[ch].trig->SetFillColor(kRed);
        // h_offsum_npe[ch].trig->SetFillStyle(3003);
        FitResult offsum_result = ana_helper::npe_gauss_fit(h_offsum_npe[ch].trig, c, nth_pad, 1.5, cutoff_threshold);
        result_container["offsum_npe"].push_back(offsum_result);
        h_offsum_npe_shower[ch]->SetLineColor(kGreen);
        h_offsum_npe_shower[ch]->SetFillColor(kGreen);
        h_offsum_npe_shower[ch]->SetFillStyle(3003);
        h_offsum_npe_shower[ch]->Draw("same");
        nth_pad++;

        c->cd(nth_pad);
        // h_onsum_npe[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        // h_onsum_npe[ch].raw->SetLineColor(kBlue);
        // h_onsum_npe[ch].raw->Draw();
        // h_onsum_npe[ch].trig->SetLineColor(kRed);
        // h_onsum_npe[ch].trig->SetFillColor(kRed);
        // h_onsum_npe[ch].trig->SetFillStyle(3003);
        FitResult onsum_result = ana_helper::npe_gauss_fit(h_onsum_npe[ch].trig, c, nth_pad, 1.5, cutoff_threshold);
        result_container["onsum_npe"].push_back(onsum_result);
        h_onsum_npe_shower[ch]->SetLineColor(kGreen);
        h_onsum_npe_shower[ch]->SetFillColor(kGreen);
        h_onsum_npe_shower[ch]->SetFillStyle(3003);
        h_onsum_npe_shower[ch]->Draw("same");
        nth_pad++;
    }

    c->Print(pdf_name);
    if (start_or_end == 2 || start_or_end == 3) c->Print(pdf_name + "]"); // end
    gROOT->Reset();

    return result_container;
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    conf.kvc_thick_initialize();

    // // +-------------+
    // // | dev version |
    // // +-------------+
    // // -- check argments -----
    // if (argc < 2) {
    //     std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
    //     return 1;
    // }
    // Int_t run_num = std::atoi(argv[1]);

    // analyze(run_num, 3);
    

    // +-------------+
    // | pro version |
    // +-------------+
    std::vector<Int_t> ana_run_num{ 
        // Condition 2
        498, 499, 500, 501, 502, 503, 504,
        488, 489, 490, 491, 492, 493, 494,
        448, 449, 450, 451, 452, 453, 454,
        458, 459, 460, 461, 462, 463, 464,
        468, 469, 470, 471, 472, 473, 474,
        478, 479, 480, 481, 482, 483, 484,
    };

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = OUTPUT_DIR + "/root/kvc_thick_pos_scan_analysis.root";
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", ""); 

    // -- prepare root file branch -----
    Int_t tmp_run_num, pos_x, pos_y;
    Double_t n_trig, n_hit;
    std::vector<Double_t> linear_a, linear_b;
    std::vector<Double_t> indiv_npe_val, indiv_npe_err, onsum_npe_val, onsum_npe_err, offsum_npe_val, offsum_npe_err;

    output_tree.Branch("run_num", &tmp_run_num, "run_num/I");
    output_tree.Branch("pos_x", &pos_x, "pos_x/I");
    output_tree.Branch("pos_y", &pos_y, "pos_y/I");
    output_tree.Branch("n_trig", &n_trig, "n_trig/D");
    output_tree.Branch("n_hit", &n_hit, "n_hit/D");
    output_tree.Branch("linear_a", &linear_a);
    output_tree.Branch("linear_b", &linear_b);
    output_tree.Branch("indiv_npe_val", &indiv_npe_val);
    output_tree.Branch("indiv_npe_err", &indiv_npe_err);
    output_tree.Branch("onsum_npe_val", &onsum_npe_val);
    output_tree.Branch("onsum_npe_err", &onsum_npe_err);
    output_tree.Branch("offsum_npe_val", &offsum_npe_val);
    output_tree.Branch("offsum_npe_err", &offsum_npe_err);

    for (Int_t i = 0, n_run_num = ana_run_num.size(); i < n_run_num; i++) {
        tmp_run_num = ana_run_num[i];
        std::pair<Int_t, Int_t> position = ana_helper::get_scan_position(ana_run_num[i]);
        pos_x = position.first;
        pos_y = position.second;
        
        // -- analyze -----
        Int_t pdf_save_mode = 0;
        if (i == 0) pdf_save_mode = 1;
        else if (i == n_run_num-1) pdf_save_mode = 2;
        std::unordered_map<std::string, std::vector<FitResult>> result_container = analyze(ana_run_num[i], pdf_save_mode);

        // -- initialize -----
        linear_a.clear(); linear_b.clear();
        indiv_npe_val.clear(); indiv_npe_err.clear();
        onsum_npe_val.clear(); onsum_npe_err.clear();
        offsum_npe_val.clear(); offsum_npe_err.clear();

        // -- efficiency -----
        n_trig = result_container["eff"][0].additional[0];
        n_hit  = result_container["eff"][0].additional[1];
        
        // -- linear -----
        for (const auto &result : result_container["linear"]) {
            linear_a.push_back(result.par.size() > 0 ? result.par[0] : 0.0);  // par[0] が存在しない場合は 0.0 を詰める
            linear_b.push_back(result.par.size() > 1 ? result.par[1] : 0.0);  // par[1] が存在しない場合も 0.0 を詰める
        }

        // -- indiv -----
        for (const auto &result : result_container["indiv_npe"]) {
            indiv_npe_val.push_back(result.par.size() > 1 ? result.par[1] : 0.0);  // par[1] がなければ 0.0
            indiv_npe_err.push_back(result.err.size() > 1 ? result.err[1] : 0.0);  // err[1] がなければ 0.0
        }

        // -- onsum -----
        for (const auto &result : result_container["onsum_npe"]) {
            onsum_npe_val.push_back(result.par.size() > 1 ? result.par[1] : 0.0);  // par[1] がなければ 0.0
            onsum_npe_err.push_back(result.err.size() > 1 ? result.err[1] : 0.0);  // err[1] がなければ 0.0
        }

        for (const auto &result : result_container["offsum_npe"]) {
            offsum_npe_val.push_back(result.par.size() > 1 ? result.par[1] : 0.0);  // par[1] がなければ 0.0
            offsum_npe_err.push_back(result.err.size() > 1 ? result.err[1] : 0.0);  // err[1] がなければ 0.0
        }


        output_tree.Fill();
    }

    // +------------+
    // | Write data |
    // +------------+
    fout.cd(); // 明示的にカレントディレクトリを設定
    output_tree.Write();
    fout.Close(); 


    return 0;
}