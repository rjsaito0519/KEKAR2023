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

static const TString pdf_name =  OUTPUT_DIR + "/img/bac_hv_threshold_scan.pdf";

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
    TTreeReaderArray<Double_t> baca_beam(reader_beam, "baca");
    TTreeReaderArray<Double_t> bacsuma_beam(reader_beam, "bacsuma");
    TTreeReaderArray<Double_t> bacsumt_beam(reader_beam, "bacsumt");

    // -- w/o beam data ------
    root_file_path = Form("%s/kekar_run%05d.root", DATA_DIR.Data(), run_num_ped);
    auto *f_ped = new TFile( root_file_path.Data() );
    if (!f_ped || f_ped->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return {};
    }
    TTreeReader reader_ped("tree", f_ped);
    TTreeReaderArray<Double_t> baca_ped(reader_ped, "baca");
    TTreeReaderArray<Double_t> bacsuma_ped(reader_ped, "bacsuma");
    TTreeReaderArray<Double_t> bacsumt_ped(reader_ped, "bacsumt");


    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    std::vector<HistPair> h_baca;
    TH1D *h_baca_ped[conf.max_bac_ch]; 
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        TString name  = Form("BACa_%d_%d", run_num_beam, ch + 1);
        TString title = Form("run%05d BAC(ADC) ch%d;ADC;", run_num_beam, ch + 1);
        h_baca.emplace_back(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);

        name  = Form("BACa_%d_%d", run_num_ped, ch + 1);
        title = Form("run%05d BAC(pedestal) ch%d;ADC;", run_num_ped, ch + 1);
        h_baca_ped[ch] = new TH1D(name, title, conf.adjust_adc_bin_num, conf.adc_min, conf.adc_max);
    }
    auto *h_onsum_ped = new TH1D(Form("BAConSUMa_%d", run_num_beam), Form("run%05d BAC online sum(ADC);ADC;", run_num_beam), conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);
    HistPair h_onsum_adc(Form("BAConSUMa_%d", run_num_beam), Form("run%05d BAC online sum(ADC);ADC;", run_num_beam), conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);
    HistPair h_offsum_adc(Form("BACoffSUMa_%d", run_num_beam), Form("run%05d BAC offline sum(ADC);ADC;", run_num_beam), conf.sumadc_bin_num, conf.sumadc_min, conf.sumadc_max);
    auto *h_bacsumt = new TH1D(Form("BACSUMt_%d", run_num_beam), Form("run%05d BACSUM(TDC);TDC;", run_num_beam), conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    
    // -- NPE -----
    std::vector<HistPair> h_npe;
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        TString name  = Form("BACnpe_%d_%d", run_num_beam, ch + 1);
        TString title = Form("run%05d BAC(NPE) ch%d;NPE;", run_num_beam, ch + 1);
        h_npe.emplace_back(name, title, conf.npe_bin_num, conf.npe_min, conf.npe_max);
    }
    HistPair h_onsum_npe(Form("BAConsumnpe_%d", run_num_beam), Form("run%05d BAC online SUM(NPE);NPE;", run_num_beam), conf.npe_bin_num, conf.npe_min, conf.npe_max);
    HistPair h_offsum_npe(Form("BACoffsumnpe_%d", run_num_beam), Form("run%05d BAC offline sum(NPE);NPE;", run_num_beam), conf.npe_bin_num, conf.npe_min, conf.npe_max);
   
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
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            h_bacsumt->Fill(bacsumt_beam[n_hit]);
        }
    }
    // -- for pedestal -----
    reader_ped.Restart();
    while (reader_ped.Next()){
        for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) h_baca_ped[ch]->Fill(baca_ped[ch]);
        h_onsum_ped->Fill(bacsuma_ped[0]);
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

    // -- BAC pedestal -----
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        FitResult result = ana_helper::pedestal_fit_with_gauss(h_baca_ped[ch], c, nth_pad);
        result_container["BAC_ped"].push_back(result);
        nth_pad++;
    }
    {
        c->Print(pdf_name);
        c->Clear();
        c->Divide(1, 1);
        nth_pad = 1;
        FitResult result = ana_helper::pedestal_fit_with_gauss(h_onsum_ped, c, nth_pad);
        result_container["BACSUM_ped"].push_back(result);
    }

    // -- bac sum TDC -----
    c->Print(pdf_name);
    c->Clear();
    c->Divide(1, 1);
    nth_pad = 1;
    FitResult tmp_fit_result = ana_helper::cherenkov_tdc_fit(h_bacsumt, c, nth_pad);
    Double_t bac_tdc_min = tmp_fit_result.additional[0];
    Double_t bac_tdc_max = tmp_fit_result.additional[1];



    // +------------------+
    // | Fill event (2nd) |
    // +------------------+
    // -- for w/ beam data -----
    Int_t n_hitbac_beam = 0;
    std::unordered_map<std::string, std::vector<Double_t>> online_sum_container{ {"raw", {}}, {"trig", {}} };
    reader_beam.Restart();
    while (reader_beam.Next()){
        // -- set up flag -----
        Bool_t do_hit_bac = false;
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            // if ( bac_tdc_min < bacsumt_beam[n_hit] && bacsumt_beam[n_hit] < bac_tdc_max ) do_hit_bac = true;
            if ( conf.hv_th_scan_bac_tdc_gate_min < bacsumt_beam[n_hit] && bacsumt_beam[n_hit] < conf.hv_th_scan_bac_tdc_gate_max) {
                do_hit_bac = true;
                n_hitbac_beam++;
            }
        }

        // -- event selection and fill data -----
        Double_t offline_sum_a   = 0.;
        Double_t offline_sum_npe = 0.;
        for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
            h_baca[ch].raw->Fill(baca_beam[ch]);
            h_npe[ch].raw->Fill( (baca_beam[ch] - result_container["BAC_ped"][ch].par[1]) / conf.bac_opg[hv][ch].first );
            offline_sum_a   +=  baca_beam[ch] - result_container["BAC_ped"][ch].par[1];
            offline_sum_npe += (baca_beam[ch] - result_container["BAC_ped"][ch].par[1]) / conf.bac_opg[hv][ch].first;
        }
        h_offsum_adc.raw->Fill(offline_sum_a);
        h_onsum_adc.raw->Fill(bacsuma_beam[0] - result_container["BACSUM_ped"][0].par[1]);
        h_offsum_npe.raw->Fill(offline_sum_npe);
        online_sum_container["raw"].push_back(bacsuma_beam[0]);
        
        if (do_hit_bac) {
            offline_sum_a   = 0.;
            offline_sum_npe = 0.;
            for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
                h_baca[ch].trig->Fill(baca_beam[ch]);
                h_npe[ch].trig->Fill( (baca_beam[ch] - result_container["BAC_ped"][ch].par[1]) / conf.bac_opg[hv][ch].first );
                offline_sum_a   +=  baca_beam[ch] - result_container["BAC_ped"][ch].par[1];
                offline_sum_npe += (baca_beam[ch] - result_container["BAC_ped"][ch].par[1]) / conf.bac_opg[hv][ch].first;
            }
            h_offsum_adc.trig->Fill(offline_sum_a);
            h_onsum_adc.trig->Fill(bacsuma_beam[0] - result_container["BACSUM_ped"][0].par[1]);
            h_offsum_npe.trig->Fill(offline_sum_npe);
            online_sum_container["trig"].push_back(bacsuma_beam[0]);

            // 最初にOnline sumのone photon gainをoffline sumとの相関より推測する
            h_correlation->Fill( bacsuma_beam[0] - result_container["BACSUM_ped"][0].par[1], offline_sum_npe );
        }
    }
    FitResult eff_result_beam;
    eff_result_beam.additional.push_back( static_cast<Double_t>( reader_beam.GetEntries() ) );
    eff_result_beam.additional.push_back( static_cast<Double_t>(n_hitbac_beam) );
    result_container["eff"].push_back(eff_result_beam);
    std::cout << (Double_t) n_hitbac_beam / reader_beam.GetEntries() << std::endl;

    // -- for w/0 beam data (check noise level) -----
    Int_t n_hitbac_ped = 0;
    reader_ped.Restart();
    while (reader_ped.Next()){
        // -- set up flag -----
        Bool_t do_hit_bac = false;
        for (Int_t n_hit = 0; n_hit < conf.max_nhit_tdc; n_hit++) {
            // if ( bac_tdc_min < bacsumt_ped[n_hit] && bacsumt_ped[n_hit] < bac_tdc_max ) do_hit_bac = true;
            if ( conf.hv_th_scan_bac_tdc_gate_min < bacsumt_ped[n_hit] && bacsumt_ped[n_hit] < conf.hv_th_scan_bac_tdc_gate_max) {
                n_hitbac_ped++;
            }
        }
    }
    FitResult eff_result_ped;
    eff_result_ped.additional.push_back( static_cast<Double_t>( reader_ped.GetEntries() ) );
    eff_result_ped.additional.push_back( static_cast<Double_t>(n_hitbac_ped) );
    result_container["eff"].push_back(eff_result_ped);
    std::cout << (Double_t) n_hitbac_ped / reader_ped.GetEntries() << std::endl;


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
    for (const auto& it : online_sum_container["raw"]) h_onsum_npe.raw->Fill( linear_fit_result.par[0]*(it - result_container["BACSUM_ped"][0].par[1]) + linear_fit_result.par[1] );
    for (const auto& it : online_sum_container["trig"]) h_onsum_npe.trig->Fill( linear_fit_result.par[0]*(it - result_container["BACSUM_ped"][0].par[1]) + linear_fit_result.par[1] );

    //  +-----------------+
    //  | Draw histograms |
    //  +-----------------+
    // -- indiv adc -----
    c->Print(pdf_name);
    c->Clear();
    c->Divide(cols, rows);
    nth_pad = 1;

    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        c->cd(nth_pad);
        // gPad->SetLogy(1);
        Double_t peak_pos = h_baca[ch].raw->GetMean();
        Double_t stdev = h_baca[ch].raw->GetStdDev();
        h_baca[ch].raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_baca[ch].raw->SetLineColor(kBlue);
        h_baca[ch].raw->Draw();
        h_baca[ch].trig->SetLineColor(kRed);
        h_baca[ch].trig->SetFillColor(kRed);
        h_baca[ch].trig->SetFillStyle(3003);
        h_baca[ch].trig->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_baca[ch].trig->Draw("same");
        nth_pad++;
    }

    // -- indiv NPE -----
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) {
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
    {
        c->Print(pdf_name);
        c->Clear();
        c->Divide(2, 1);
        nth_pad = 1;

        c->cd(nth_pad);
        Double_t peak_pos = h_offsum_adc.raw->GetMean();
        Double_t stdev = h_offsum_adc.raw->GetStdDev();
        h_offsum_adc.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_offsum_adc.raw->SetLineColor(kBlue);
        h_offsum_adc.raw->Draw();
        h_offsum_adc.trig->SetLineColor(kRed);
        h_offsum_adc.trig->SetFillColor(kRed);
        h_offsum_adc.trig->SetFillStyle(3003);
        h_offsum_adc.trig->Draw("same");
        nth_pad++;

        c->cd(nth_pad);
        h_onsum_adc.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_onsum_adc.raw->SetLineColor(kBlue);
        h_onsum_adc.raw->Draw();
        h_onsum_adc.trig->SetLineColor(kRed);
        h_onsum_adc.trig->SetFillColor(kRed);
        h_onsum_adc.trig->SetFillStyle(3003);
        h_onsum_adc.trig->Draw("same");
        nth_pad++;
    }

    // -- sum npe -----
    {
        c->Print(pdf_name);
        c->Clear();
        c->Divide(2, 1);
        nth_pad = 1;
        
        c->cd(nth_pad);
        Double_t peak_pos = h_offsum_npe.raw->GetMean();
        Double_t stdev = h_offsum_npe.raw->GetStdDev();
        h_offsum_npe.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_offsum_npe.raw->SetLineColor(kBlue);
        h_offsum_npe.raw->Draw();
        h_offsum_npe.trig->SetLineColor(kRed);
        h_offsum_npe.trig->SetFillColor(kRed);
        h_offsum_npe.trig->SetFillStyle(3003);
        h_offsum_npe.trig->Draw("same");
        nth_pad++;

        c->cd(nth_pad);
        h_onsum_npe.raw->GetXaxis()->SetRangeUser(peak_pos - 5.0*stdev, peak_pos + 5.0*stdev);
        h_onsum_npe.raw->SetLineColor(kBlue);
        h_onsum_npe.raw->Draw();
        h_onsum_npe.trig->SetLineColor(kRed);
        h_onsum_npe.trig->SetFillColor(kRed);
        h_onsum_npe.trig->SetFillStyle(3003);
        h_onsum_npe.trig->Draw("same");
        nth_pad++;
    }

    // -- check threshold -----
    {
        c->Print(pdf_name);
        c->Clear();
        c->Divide(2, 1);
        nth_pad = 1;

        c->cd(nth_pad);
        h_offsum_ratio->Divide( h_offsum_npe.trig, h_offsum_npe.raw, 1, 1 );
        // h_offsum_ratio->GetXaxis()->SetRangeUser(conf.threshold_fit_range_min - 10, conf.threshold_fit_range_max + 10);
        // h_offsum_ratio->Draw();
        FitResult offsum_result = ana_helper::threshold_erf_fit(h_offsum_ratio, c, nth_pad);
        result_container["offsum_thre"].push_back(offsum_result);
        nth_pad++;

        c->cd(nth_pad);
        h_onsum_ratio->Divide( h_onsum_npe.trig, h_onsum_npe.raw, 1, 1 );
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
    conf.bac_initialize(525);

    // +-------------+
    // | pro version |
    // +-------------+
    std::vector<std::vector<Double_t>> run_info{
    //   beam ped  HV  Vth
    
        // Condition 1
        {269, 295, 58,  50},
        {270, 296, 58,  75},
        {271, 297, 58, 100},
        {272, 298, 58, 125},
        {273, 299, 58, 150},
        
        {274, 290, 57,  50},
        {275, 291, 57,  75},
        {276, 292, 57, 100},
        {277, 293, 57, 125},
        {278, 294, 57, 150},
        
        {280, 285, 56,  50},
        {281, 286, 56,  75},
        {282, 287, 56, 100},
        {283, 288, 56, 125},
        {284, 289, 56, 150},

        // Condition 2
        {412, 429, 58,  50},
        {445, 446, 58,  60},
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
    TString output_path = OUTPUT_DIR + "/root/bac_hv_threshold_scan.root";
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", ""); 

    // -- prepare root file branch -----
    Int_t run_num_beam, hv, threshold;
    Double_t n_entry_beam, n_entry_ped;
    std::vector<Double_t> linear_a, linear_b, n_hit_beam, n_hit_ped;
    std::vector<Double_t> onsum_thre_val, onsum_thre_err, offsum_thre_val, offsum_thre_err;

    output_tree.Branch("run_num", &run_num_beam, "run_num/I");
    output_tree.Branch("hv", &hv, "hv/I");
    output_tree.Branch("threshold", &threshold, "threshold/I");
    output_tree.Branch("n_entry_beam", &n_entry_beam, "n_entry_beam/D");
    output_tree.Branch("n_hit_beam", &n_hit_beam);
    output_tree.Branch("n_entry_ped", &n_entry_ped, "n_entry_ped/D");
    output_tree.Branch("n_hit_ped", &n_hit_ped);    
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
        linear_a.clear(); linear_b.clear(); n_hit_beam.clear(); n_hit_ped.clear();
        onsum_thre_val.clear(); onsum_thre_err.clear();
        offsum_thre_val.clear(); offsum_thre_err.clear();

        // -- efficiency -----
        n_entry_beam = result_container["eff"][0].additional[0];
        n_entry_ped  = result_container["eff"][1].additional[0];
        for (Int_t j = 1, n = result_container["eff"][0].additional.size(); j < n; j++) {
            n_hit_beam.push_back(result_container["eff"][0].additional[j]);
            n_hit_ped.push_back(result_container["eff"][1].additional[j]);        
        }

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