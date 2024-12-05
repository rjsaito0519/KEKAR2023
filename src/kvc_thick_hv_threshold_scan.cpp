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

static const TString pdf_name =  OUTPUT_DIR + "/img/kvc_2cm_hv_threshold_scan.pdf";

std::unordered_map<std::string, std::vector<FitResult>> analyze(Int_t run_num_beam, Int_t run_num_ped, Int_t hv, Int_t start_or_end = 0) // 0: mid_page, 1:start_page, 2: end_page, 3: both
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
    // -- w/ beam data ------
    TString root_file_path = Form("%s/kekar_run%05d.root", DATA_DIR.Data(), run_num_beam);
    auto *f_beam = new TFile( root_file_path.Data() );
    if (!f_beam || f_beam->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return {};
    }
    TTreeReader reader_beam("tree", f_beam);
    TTreeReaderArray<Double_t> t1a(reader_beam, "t1a");
    TTreeReaderArray<Double_t> t1t(reader_beam, "t1t");
    TTreeReaderArray<Double_t> t2a(reader_beam, "t2a");
    TTreeReaderArray<Double_t> t2t(reader_beam, "t2t");
    TTreeReaderArray<Double_t> t3a(reader_beam, "t3a");
    TTreeReaderArray<Double_t> t3t(reader_beam, "t3t");
    TTreeReaderArray<Double_t> t4a(reader_beam, "t4a");
    TTreeReaderArray<Double_t> t4t(reader_beam, "t4t");
    TTreeReaderArray<Double_t> kvca_beam(reader_beam, "kvca");
    TTreeReaderArray<Double_t> kvcsuma_beam(reader_beam, "kvcsuma");
    TTreeReaderArray<Double_t> kvcsumt_beam(reader_beam, "kvcsumt");

    // -- w/o beam data ------
    root_file_path = Form("%s/kekar_run%05d.root", DATA_DIR.Data(), run_num_ped);
    auto *f_ped = new TFile( root_file_path.Data() );
    if (!f_ped || f_ped->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return {};
    }
    TTreeReader reader_ped("tree", f_ped);
    TTreeReaderArray<Double_t> kvca_ped(reader_ped, "kvca");
    TTreeReaderArray<Double_t> kvcsuma_ped(reader_ped, "kvcsuma");
    TTreeReaderArray<Double_t> kvcsumt_ped(reader_ped, "kvcsumt");


    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    auto *h_t1a = new TH1D(Form("T1a_%d", run_num_beam), Form("run%05d T1(ADC);ADC;", run_num_beam), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t1t = new TH1D(Form("T1t_%d", run_num_beam), Form("run%05d T1(TDC);TDC;", run_num_beam), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    auto *h_t2a = new TH1D(Form("T2a_%d", run_num_beam), Form("run%05d T2(ADC);ADC;", run_num_beam), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t2t = new TH1D(Form("T2t_%d", run_num_beam), Form("run%05d T2(TDC);TDC;", run_num_beam), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    auto *h_t3a = new TH1D(Form("T3a_%d", run_num_beam), Form("run%05d T3(ADC);ADC;", run_num_beam), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t3t = new TH1D(Form("T3t_%d", run_num_beam), Form("run%05d T3(TDC);TDC;", run_num_beam), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    auto *h_t4a = new TH1D(Form("T4a_%d", run_num_beam), Form("run%05d T4(ADC);ADC;", run_num_beam), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_t4t = new TH1D(Form("T4t_%d", run_num_beam), Form("run%05d T4(TDC);TDC;", run_num_beam), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);

    std::vector<HistPair> h_kvca;
    std::vector<HistPair> h_onsum_adc;
    std::vector<HistPair> h_offsum_adc;
    TH1D *h_kvca_ped[conf.max_kvc_ch];
    TH1D *h_onsum_ped[conf.max_kvc_ch];
    TH1D *h_kvcsumt[conf.max_kvc_ch];
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        TString name  = Form("KVCa_%d_%d", run_num_beam, ch + 1);
        TString title = Form("run%05d KVC(ADC) seg2 %s %d;ADC;", run_num_beam, ch/2==0 ? "UP" : "DOWN", ch%2+1);
        h_kvca.emplace_back(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);

        name  = Form("KVConSUMa_%d_%d", run_num_beam, ch+1);
        title = Form("run%05d KVC online sum(ADC) seg%d;ADC;", run_num_beam, ch+1);
        h_onsum_adc.emplace_back(name, title, conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);

        name  = Form("KVCoffSUMa_%d_%d", run_num_beam, ch+1);
        title = Form("run%05d KVC offline sum(ADC) seg%d;ADC;", run_num_beam, ch+1);
        h_offsum_adc.emplace_back(name, title, conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);

        name  = Form("KVCa_%d_%d", run_num_ped, ch + 1);
        title = Form("run%05d KVC(pedestal) seg2 %s %d;ADC;", run_num_ped, ch/2==0 ? "UP" : "DOWN", ch%2+1);
        h_kvca_ped[ch] = new TH1D(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);

        name  = Form("KVCSUMa_%d_%d", run_num_ped, ch + 1);
        title = Form("run%05d KVC SUM(pedestal) ch%d;ADC;", run_num_ped, ch + 1);
        h_onsum_ped[ch] = new TH1D(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);

        h_kvcsumt[ch] = new TH1D(Form("KVCSUMt_%d_%d", run_num_beam, ch+1), Form("run%05d KVCSUM(TDC) seg%d;TDC;", run_num_beam, ch+1), conf.adjust_tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }
    
    // -- NPE -----
    std::vector<HistPair> h_npe;
    std::vector<HistPair> h_onsum_npe;
    std::vector<HistPair> h_offsum_npe;
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        TString name  = Form("KVCnpe_%d_%d", run_num_beam, ch + 1);
        TString title = Form("run%05d KVC(NPE) seg2 %s %d;NPE;", run_num_beam, ch/2==0 ? "UP" : "DOWN", ch%2+1);
        h_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);

        name  = Form("KVConsumnpe_%d_%d", run_num_beam, ch+1);
        title = Form("run%05d KVC online SUM(NPE) seg%d;NPE;", run_num_beam, ch+1);
        h_onsum_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);

        name  = Form("KVCoffsumnpe_%d_%d", run_num_beam, ch+1);
        title = Form("run%05d KVC offline SUM(NPE) seg%d;NPE;", run_num_beam, ch+1);
        h_offsum_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);
    }

    // -- 2d histogram -----
    auto *h_correlation = new TH2D("on_off_correlation", ";online sum[adc];offline sum[npe]", conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max, conf.npe_bin_num, conf.npe_min, conf.npe_max);

    // -- threshold -----
    auto *h_onsum_ratio = new TH1D("onsum_ratio", Form("run%05d online sum ratio;NPE;ratio", run_num_beam), conf.npe_bin_num, conf.npe_min, conf.npe_max);
    auto *h_offsum_ratio = new TH1D("offsum_ratio", Form("run%05d offline sum ratio;NPE;ratio", run_num_beam), conf.npe_bin_num, conf.npe_min, conf.npe_max);


    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    // -- for TDC -----
    reader_beam.Restart();
    while (reader_beam.Next()){
        h_t1a->Fill(t1a[0]);
        h_t2a->Fill(t2a[0]);
        h_t3a->Fill(t3a[0]);
        h_t4a->Fill(t4a[0]);
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            h_t1t->Fill(t1t[n_hit]);
            h_t2t->Fill(t2t[n_hit]);
            h_t3t->Fill(t3t[n_hit]);
            h_t4t->Fill(t4t[n_hit]);
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) h_kvcsumt[ch]->Fill(kvcsumt_beam[conf.max_nhit_tdc*ch+n_hit]);
        }
    }
    // -- for pedestal -----
    reader_ped.Restart();
    while (reader_ped.Next()){
        for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
            h_kvca_ped[ch]->Fill(kvca_ped[ch]);
            h_onsum_ped[ch]->Fill(kvcsuma_ped[ch]);
        }
    }


    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- make canvas and draw -----
    Int_t nth_pad = 1;
    Int_t rows = 2;
    Int_t cols = 2;
    Int_t max_pads = rows * cols;
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
    nth_pad++;


    // -- KVC pedestal -----
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        FitResult result = ana_helper::pedestal_fit_with_gauss(h_kvca_ped[ch], c, nth_pad);
        result_container["KVC_ped"].push_back(result);
        nth_pad++;
    }
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        FitResult result = ana_helper::pedestal_fit_with_gauss(h_onsum_ped[ch], c, nth_pad);
        result_container["KVCSUM_ped"].push_back(result);
        nth_pad++;
    }

    // -- kvc sum TDC -----
    Double_t kvc_tdc_min[conf.max_kvc_ch];
    Double_t kvc_tdc_max[conf.max_kvc_ch];
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        tmp_fit_result = ana_helper::cherenkov_tdc_fit(h_kvcsumt[ch], c, nth_pad);
        kvc_tdc_min[ch] = tmp_fit_result.additional[0];
        kvc_tdc_max[ch] = tmp_fit_result.additional[1];
        nth_pad++;
    }


    // +------------------+
    // | Fill event (2nd) |
    // +------------------+
    // -- for w/ beam data -----
    std::vector<Int_t> n_hitkvc_beam(conf.max_kvc_ch, 0);
    std::unordered_map<std::string, std::vector<std::vector<Double_t>>> online_sum_container{ {"raw", {{}, {}, {}, {}}}, {"trig", {{}, {}, {}, {}}} };
    reader_beam.Restart();
    while (reader_beam.Next()){
        // -- set up flag -----
        Bool_t do_hit_t1a = false, do_hit_t2a = false, do_hit_t3a = false, do_hit_t4a = false;
        Bool_t do_hit_t1t = false, do_hit_t2t = false, do_hit_t3t = false, do_hit_t4t = false, do_hit_kvc = false;
        if ( t1_adc_min < t1a[0] && t1a[0] < t1_adc_max ) do_hit_t1a = true;
        if ( t2_adc_min < t2a[0] && t2a[0] < t2_adc_max ) do_hit_t2a = true;
        if ( t3_adc_min < t3a[0] && t3a[0] < t3_adc_max ) do_hit_t3a = true;
        if ( t4_adc_min < t4a[0] && t4a[0] < t4_adc_max ) do_hit_t4a = true;

        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            if ( t1_tdc_min < t1t[n_hit] && t1t[n_hit] < t1_tdc_max ) do_hit_t1t = true;
            if ( t2_tdc_min < t2t[n_hit] && t2t[n_hit] < t2_tdc_max ) do_hit_t2t = true;
            if ( t3_tdc_min < t3t[n_hit] && t3t[n_hit] < t3_tdc_max ) do_hit_t3t = true;
            if ( t4_tdc_min < t4t[n_hit] && t4t[n_hit] < t4_tdc_max ) do_hit_t4t = true;
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
                // if ( kvc_tdc_min[ch] < kvcsumt_beam[conf.max_nhit_tdc*ch+n_hit] && kvcsumt_beam[conf.max_nhit_tdc*ch+n_hit] < kvc_tdc_max[ch] ) do_hit_kvc = true;
                
                if ( conf.hv_th_scan_kvc_tdc_gate_min < kvcsumt_beam[conf.max_nhit_tdc*ch+n_hit] && kvcsumt_beam[conf.max_nhit_tdc*ch+n_hit] < conf.hv_th_scan_kvc_tdc_gate_min ) {
                    do_hit_kvc = true;
                    n_hitkvc_beam[ch]++;
                }
            }
        }
        Bool_t trig_flag_adc = do_hit_t4a && do_hit_t2a && do_hit_t3a && do_hit_t4a;
        Bool_t trig_flag_tdc = do_hit_t4t && do_hit_t2t && do_hit_t3t && do_hit_t4t;


        // -- event selection and fill data -----
        Double_t offline_sum_a   = 0.;
        Double_t offline_sum_npe = 0.;
        for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
            Int_t index_opg = conf.max_kvc_ch * ch + 1; // seg2
            h_kvca[ch].raw->Fill(kvca_beam[ch]);
            h_npe[ch].raw->Fill( (kvca_beam[ch] - result_container["KVC_ped"][ch].par[1]) / conf.kvc_thick_opg[hv][index_opg].first );
            offline_sum_a   +=  kvca_beam[ch] - result_container["KVC_ped"][ch].par[1];
            offline_sum_npe += (kvca_beam[ch] - result_container["KVC_ped"][ch].par[1]) / conf.kvc_thick_opg[hv][index_opg].first;
            h_onsum_adc[ch].raw->Fill(kvcsuma_beam[ch] - result_container["KVCSUM_ped"][ch].par[1]);
            online_sum_container["raw"][ch].push_back(kvcsuma_beam[ch]);
        }
        h_offsum_adc[1].raw->Fill(offline_sum_a);
        h_offsum_npe[1].raw->Fill(offline_sum_npe);
        
        if (do_hit_kvc) {
            offline_sum_a   = 0.;
            offline_sum_npe = 0.;
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
                Int_t index_opg = conf.max_kvc_ch * ch + 1; // seg2
                h_kvca[ch].trig->Fill(kvca_beam[ch]);
                h_npe[ch].trig->Fill( (kvca_beam[ch] - result_container["KVC_ped"][ch].par[1]) / conf.kvc_thick_opg[hv][index_opg].first );
                offline_sum_a   +=  kvca_beam[ch] - result_container["KVC_ped"][ch].par[1];
                offline_sum_npe += (kvca_beam[ch] - result_container["KVC_ped"][ch].par[1]) / conf.kvc_thick_opg[hv][index_opg].first;
                h_onsum_adc[ch].trig->Fill(kvcsuma_beam[ch] - result_container["KVCSUM_ped"][ch].par[1]);
                online_sum_container["trig"][ch].push_back(kvcsuma_beam[ch]);
            }
            h_offsum_adc[1].trig->Fill(offline_sum_a);
            h_offsum_npe[1].trig->Fill(offline_sum_npe);

            // 最初にOnline sumのone photon gainをoffline sumとの相関より推測する
            if (trig_flag_adc && trig_flag_tdc) h_correlation->Fill( kvcsuma_beam[1] - result_container["KVCSUM_ped"][1].par[1], offline_sum_npe );
        }
    }
    FitResult eff_result_beam;
    eff_result_beam.additional.push_back( static_cast<Double_t>( reader_beam.GetEntries() ) );
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) eff_result_beam.additional.push_back( static_cast<Double_t>(n_hitkvc_beam[ch]) );
    result_container["eff"].push_back(eff_result_beam);
    std::cout << (Double_t) n_hitkvc_beam[1] / reader_beam.GetEntries() << std::endl;

    // -- for w/0 beam data (check noise level) -----
    std::vector<Int_t> n_hitkvc_ped(conf.max_kvc_ch, 0);
    reader_ped.Restart();
    while (reader_ped.Next()){
        // -- set up flag -----
        Bool_t do_hit_kvc = false;        
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
                // if ( kvc_tdc_min[ch] < kvcsumt_ped[conf.max_nhit_tdc*ch+n_hit] && kvcsumt_ped[conf.max_nhit_tdc*ch+n_hit] < kvc_tdc_max[ch] ) do_hit_kvc = true;

                if ( conf.hv_th_scan_kvc_tdc_gate_min < kvcsumt_ped[conf.max_nhit_tdc*ch+n_hit] && kvcsumt_ped[conf.max_nhit_tdc*ch+n_hit] < conf.hv_th_scan_kvc_tdc_gate_min ) {
                    do_hit_kvc = true;
                    n_hitkvc_ped[ch]++;
                }
            }
        }
        // if (do_hit_kvc) n_hitkvc_ped++;
    }
    FitResult eff_result_ped;
    eff_result_ped.additional.push_back( static_cast<Double_t>( reader_ped.GetEntries() ) );
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) eff_result_ped.additional.push_back( static_cast<Double_t>(n_hitkvc_ped[ch]) );
    result_container["eff"].push_back(eff_result_ped);
    std::cout << (Double_t) n_hitkvc_ped[ch] / reader_ped.GetEntries() << std::endl;


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
            //           * (conf.kvc_thick_opg[hv][1].first + conf.kvc_thick_opg[hv][5].first + conf.kvc_thick_opg[hv][9].first + conf.kvc_thick_opg[hv][13].first)
            //           / (conf.kvc_thick_opg[hv][ch].first + conf.kvc_thick_opg[hv][ch+conf.max_kvc_ch].first + conf.kvc_thick_opg[hv][ch+2*conf.max_kvc_ch].first + conf.kvc_thick_opg[hv][ch+3*conf.max_kvc_ch].first);

            // coeff_b = linear_fit_result.par[1] 
            //           * (conf.kvc_thick_opg[hv][1].first + conf.kvc_thick_opg[hv][5].first + conf.kvc_thick_opg[hv][9].first + conf.kvc_thick_opg[hv][13].first)
            //           / (conf.kvc_thick_opg[hv][ch].first + conf.kvc_thick_opg[hv][ch+conf.max_kvc_ch].first + conf.kvc_thick_opg[hv][ch+2*conf.max_kvc_ch].first + conf.kvc_thick_opg[hv][ch+3*conf.max_kvc_ch].first);

            coeff_a = 4.0/(conf.kvc_thick_opg[hv][ch].first + conf.kvc_thick_opg[hv][ch+conf.max_kvc_ch].first + conf.kvc_thick_opg[hv][ch+2*conf.max_kvc_ch].first + conf.kvc_thick_opg[hv][ch+3*conf.max_kvc_ch].first);
            coeff_b = 0.0;

            // coeff_a = (1.0/conf.kvc_thick_opg[hv][ch].first + 1.0/conf.kvc_thick_opg[hv][ch+conf.max_kvc_ch].first + 1.0/conf.kvc_thick_opg[hv][ch+2*conf.max_kvc_ch].first + 1.0/conf.kvc_thick_opg[hv][ch+3*conf.max_kvc_ch].first)/4.0;
            // coeff_b = 0.0;

            // std::cout << linear_fit_result.par[0] << ", " << linear_fit_result.par[1] << std::endl;
            // std::cout << 4.0/(conf.kvc_thick_opg[hv][1].first + conf.kvc_thick_opg[hv][5].first + conf.kvc_thick_opg[hv][9].first + conf.kvc_thick_opg[hv][13].first) << std::endl;
            // std::cout << coeff_a << ", " << coeff_b << std::endl;
        }

        for (const auto& it : online_sum_container["raw"][ch]) h_onsum_npe[ch].raw->Fill( coeff_a*(it - result_container["KVCSUM_ped"][0].par[1]) + coeff_b );        
        for (const auto& it : online_sum_container["trig"][ch]) h_onsum_npe[ch].trig->Fill( coeff_a*(it - result_container["KVCSUM_ped"][0].par[1]) + coeff_b );
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
        h_kvca[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_kvca[ch].raw->SetLineColor(kBlue);
        h_kvca[ch].raw->Draw();
        h_kvca[ch].trig->SetLineColor(kRed);
        h_kvca[ch].trig->SetFillColor(kRed);
        h_kvca[ch].trig->SetFillStyle(3003);
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
        h_npe[ch].raw->GetXaxis()->SetRangeUser(-5, 100.0);
        h_npe[ch].raw->SetLineColor(kBlue);
        h_npe[ch].raw->Draw();
        h_npe[ch].trig->SetLineColor(kRed);
        h_npe[ch].trig->SetFillColor(kRed);
        h_npe[ch].trig->SetFillStyle(3003);
        h_npe[ch].trig->Draw("same");
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
        h_offsum_adc[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_offsum_adc[ch].raw->SetLineColor(kBlue);
        h_offsum_adc[ch].raw->Draw();
        h_offsum_adc[ch].trig->SetLineColor(kRed);
        h_offsum_adc[ch].trig->SetFillColor(kRed);
        h_offsum_adc[ch].trig->SetFillStyle(3003);
        h_offsum_adc[ch].trig->Draw("same");
        nth_pad++;

        c->cd(nth_pad);
        h_onsum_adc[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_onsum_adc[ch].raw->SetLineColor(kBlue);
        h_onsum_adc[ch].raw->Draw();
        h_onsum_adc[ch].trig->SetLineColor(kRed);
        h_onsum_adc[ch].trig->SetFillColor(kRed);
        h_onsum_adc[ch].trig->SetFillStyle(3003);
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

        c->cd(nth_pad);
        Double_t peak_pos = h_offsum_npe[ch].raw->GetMean();
        Double_t stdev = h_offsum_npe[ch].raw->GetStdDev();
        h_offsum_npe[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_offsum_npe[ch].raw->SetLineColor(kBlue);
        h_offsum_npe[ch].raw->Draw();
        h_offsum_npe[ch].trig->SetLineColor(kRed);
        h_offsum_npe[ch].trig->SetFillColor(kRed);
        h_offsum_npe[ch].trig->SetFillStyle(3003);
        h_offsum_npe[ch].trig->Draw("same");
        nth_pad++;

        c->cd(nth_pad);
        h_onsum_npe[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_onsum_npe[ch].raw->SetLineColor(kBlue);
        h_onsum_npe[ch].raw->Draw();
        h_onsum_npe[ch].trig->SetLineColor(kRed);
        h_onsum_npe[ch].trig->SetFillColor(kRed);
        h_onsum_npe[ch].trig->SetFillStyle(3003);
        h_onsum_npe[ch].trig->Draw("same");
        nth_pad++;
    }

    // -- check threshold -----
    {
        c->Print(pdf_name);
        c->Clear();
        c->Divide(2, 1);
        nth_pad = 1;

        c->cd(nth_pad);
        h_offsum_ratio->Divide( h_offsum_npe[1].trig, h_offsum_npe[1].raw, 1, 1 );
        // h_offsum_ratio->GetXaxis()->SetRangeUser(conf.threshold_fit_range_min - 10, conf.threshold_fit_range_max + 10);
        // h_offsum_ratio->Draw();
        FitResult offsum_result = ana_helper::threshold_erf_fit(h_offsum_ratio, c, nth_pad);
        result_container["offsum_thre"].push_back(offsum_result);
        nth_pad++;

        c->cd(nth_pad);
        h_onsum_ratio->Divide( h_onsum_npe[1].trig, h_onsum_npe[1].raw, 1, 1 );
        // h_onsum_ratio->GetXaxis()->SetRangeUser(conf.threshold_fit_range_min - 10, conf.threshold_fit_range_max + 10);
        // h_onsum_ratio->Draw();
        FitResult onsum_result = ana_helper::threshold_erf_fit(h_onsum_ratio, c, nth_pad);
        result_container["onsum_thre"].push_back(onsum_result);
        nth_pad++;
    }


    c->Print(pdf_name);
    if (start_or_end == 2 || start_or_end == 3) c->Print(pdf_name + "]"); // end
    gROOT->Reset();

    return result_container;
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();
    conf.kvc_thick_initialize(3000.0);

    // +-------------+
    // | pro version |
    // +-------------+
    std::vector<std::vector<Double_t>> run_info{
    //   beam ped  HV  Vth
    
        // Condition 2
        {412, 429, 58,  50},
        {414, 430, 58,  75},
        {415, 431, 58, 100},
        {416, 432, 58, 125},
        {417, 433, 58, 150},
        
        {418, 434, 57,  50},
        {419, 435, 57,  75},
        {420, 436, 57, 100},
        {421, 437, 57, 125},
        {422, 438, 57, 150},
        
        {423, 439, 56,  50},
        {424, 440, 56,  75},
        {425, 441, 56, 100},
        {426, 442, 56, 125},
        {427, 443, 56, 150},
    };

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = OUTPUT_DIR + "/root/kvc_2cm_hv_threshold_scan.root";
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", ""); 

    // -- prepare root file branch -----
    Int_t run_num_beam, hv, threshold;
    Double_t n_entry_beam, n_hit_beam, n_entry_ped, n_hit_ped;
    std::vector<Double_t> linear_a, linear_b;
    std::vector<Double_t> onsum_thre_val, onsum_thre_err, offsum_thre_val, offsum_thre_err;

    output_tree.Branch("run_num", &run_num_beam, "run_num/I");
    output_tree.Branch("hv", &hv, "hv/I");
    output_tree.Branch("threshold", &threshold, "threshold/I");
    output_tree.Branch("n_entry_beam", &n_entry_beam, "n_entry_beam/D");
    output_tree.Branch("n_hit_beam", &n_hit_beam, "n_hit_beam/D");
    output_tree.Branch("n_entry_ped", &n_entry_ped, "n_entry_ped/D");
    output_tree.Branch("n_hit_ped", &n_hit_ped, "n_hit_ped/D");    
    output_tree.Branch("linear_a", &linear_a);
    output_tree.Branch("linear_b", &linear_b);
    output_tree.Branch("onsum_thre_val", &onsum_thre_val);
    output_tree.Branch("onsum_thre_err", &onsum_thre_err);
    output_tree.Branch("offsum_thre_val", &offsum_thre_val);
    output_tree.Branch("offsum_thre_err", &offsum_thre_err);

    for (Int_t i = 0, n_run_info = run_info.size(); i < n_run_info; i++) {
        run_num_beam = run_info[i][0];
        hv = run_info[i][2];
        threshold = run_info[i][3];
        
        // -- analyze -----
        Int_t pdf_save_mode = 0;
        if (i == 0) pdf_save_mode = 1;
        else if (i == n_run_info-1) pdf_save_mode = 2;
        std::unordered_map<std::string, std::vector<FitResult>> result_container = analyze(run_info[i][0], run_info[i][1], run_info[i][2], pdf_save_mode);

        // -- initialize -----
        linear_a.clear(); linear_b.clear();
        onsum_thre_val.clear(); onsum_thre_err.clear();
        offsum_thre_val.clear(); offsum_thre_err.clear();

        // -- efficiency -----
        n_entry_beam = result_container["eff"][0].additional[0];
        n_hit_beam   = result_container["eff"][0].additional[1];
        n_entry_ped  = result_container["eff"][1].additional[0];
        n_hit_ped    = result_container["eff"][1].additional[1];
        

        // -- linear -----
        for (const auto &result : result_container["linear"]) {
            linear_a.push_back(  result.par[0] );
            linear_b.push_back(  result.par[1] );
        }

        // -- onsum -----
        for (const auto &result : result_container["onsum_thre"]) {
            onsum_thre_val.push_back( result.par[1] );
            onsum_thre_err.push_back( result.err[1] );
        }
        for (const auto &result : result_container["offsum_thre"]) {
            offsum_thre_val.push_back( result.par[1] );
            offsum_thre_err.push_back( result.err[1] );
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