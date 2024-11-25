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

static const TString pdf_name =  OUTPUT_DIR + "/img/pedestal.pdf";

std::unordered_map<std::string, std::vector<FitResult>> analyze(Int_t run_num, Int_t start_or_end = 0) // 0: mid_page, 1:start_page, 2: end_page, 3: both
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
    TString root_file_path = Form("%s/kekar_run%05d.root", DATA_DIR.Data(), run_num);
    auto *f = new TFile( root_file_path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return {};
    }
    TTreeReader reader("tree", f);
    TTreeReaderArray<Double_t> saca(reader, "saca");
    TTreeReaderArray<Double_t> sacsuma(reader, "sacsuma");
    TTreeReaderArray<Double_t> baca(reader, "baca");
    TTreeReaderArray<Double_t> bacsuma(reader, "bacsuma");
    TTreeReaderArray<Double_t> kvca(reader, "kvca");
    TTreeReaderArray<Double_t> kvcsuma(reader, "kvcsuma");
    

    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    TH1D *h_saca[conf.max_sac_ch];
    for (Int_t ch = 0; ch < conf.max_sac_ch; ch++) h_saca[ch] = new TH1D(Form("SACa_%d_%d", run_num, ch+1), Form("run%05d SAC(ADC) ch%d;ADC;", run_num, ch+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_sacsuma = new TH1D(Form("SACSUMa_%d", run_num), Form("run%05d SACSUM(ADC);ADC;", run_num), conf.adc_bin_num, conf.adc_min, conf.adc_max);

    TH1D *h_baca[conf.max_bac_ch];
    for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) h_baca[ch] = new TH1D(Form("BACa_%d_%d", run_num, ch+1), Form("run%05d BAC(ADC) ch%d;ADC;", run_num, ch+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    auto *h_bacsuma = new TH1D(Form("BACSUMa_%d", run_num), Form("run%05d BACSUM(ADC);ADC;", run_num), conf.adc_bin_num, conf.adc_min, conf.adc_max);

    TH1D *h_kvca[conf.max_kvc_ch];
    TH1D *h_kvcsuma[conf.max_kvc_ch];
    for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
        if (run_num < 412) h_kvca[ch] = new TH1D(Form("KVCa_%d_%d", run_num, ch+1), Form("run%05d KVC(ADC) seg%d %s;ADC;", run_num, ch/2+2, ch%2==0 ? "UP" : "DOWN"), conf.adc_bin_num, conf.adc_min, conf.adc_max);
        else h_kvca[ch] = new TH1D(Form("KVCa_%d_%d", run_num, ch+1), Form("run%05d KVC(ADC) seg2 %s %d;ADC;", run_num, ch<2 ? "UP" : "DOWN", ch%2+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);
        h_kvcsuma[ch] = new TH1D(Form("KVCSUMa_%d_%d", run_num, ch+1), Form("run%05d KVCSUM(ADC) ch%d;ADC;", run_num, ch+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);
    }

    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    reader.Restart();
    while (reader.Next()){
        for (Int_t ch = 0; ch < conf.max_sac_ch; ch++) h_saca[ch]->Fill(saca[ch]);
        h_sacsuma->Fill(sacsuma[0]);
        
        for (Int_t ch = 0; ch < conf.max_bac_ch; ch++) h_baca[ch]->Fill(baca[ch]);
        h_bacsuma->Fill(bacsuma[0]);

        for (Int_t ch = 0; ch < conf.max_kvc_ch; ch++) {
            h_kvca[ch]->Fill(kvca[ch]);
            h_kvcsuma[ch]->Fill(kvcsuma[ch]);
        }
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    
    // -- prepare container -----
    std::unordered_map<std::string, std::vector<FitResult>> result_container;

    // -- make canvas and draw -----
    Int_t nth_pad = 1;
    Int_t rows = 2;
    Int_t cols = 2;
    Int_t max_pads = rows * cols;
    auto *c = new TCanvas("", "", 1500, 1200);
    c->Divide(cols, rows);
    if (start_or_end == 1 || start_or_end == 3) c->Print(pdf_name + "["); // start

    // -- BAC -----
    for (Int_t ch = 0; ch < conf.max_bac_ch+1; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        if (ch < conf.max_bac_ch) {
            FitResult result = ana_helper::pedestal_fit_with_gauss(h_baca[ch], c, nth_pad);
            result_container["BAC"].push_back(result);
        } else {
            FitResult result = ana_helper::pedestal_fit_with_gauss(h_bacsuma, c, nth_pad);
            result_container["BACSUM"].push_back(result);
        }
        nth_pad++;
    }

    // -- SAC -----
    for (Int_t ch = 0; ch < conf.max_sac_ch+1; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        if (ch < conf.max_sac_ch) {
            FitResult result = ana_helper::pedestal_fit_with_gauss(h_saca[ch], c, nth_pad, 3.0);
            result_container["SAC"].push_back(result);
        } else { 
            FitResult result = ana_helper::pedestal_fit_with_gauss(h_sacsuma, c, nth_pad, 3.0);
            result_container["SACSUM"].push_back(result);
        }
        nth_pad++;
    }

    // -- KVC -----
    for (Int_t ch = 0; ch < 2*conf.max_kvc_ch; ch++) {
        if (nth_pad > max_pads) {
            c->Print(pdf_name);
            c->Clear();
            c->Divide(cols, rows);
            nth_pad = 1;
        }
        if (ch < conf.max_kvc_ch) {
            FitResult result = ana_helper::pedestal_fit_with_gauss(h_kvca[ch], c, nth_pad);
            result_container["KVC"].push_back(result);
        } else {
            FitResult result = ana_helper::pedestal_fit_with_gauss(h_kvcsuma[ch%conf.max_kvc_ch], c, nth_pad);
            result_container["KVCSUM"].push_back(result);
        }
        nth_pad++;
    }

    c->Print(pdf_name);
    if (start_or_end == 2 || start_or_end == 3) c->Print(pdf_name + "]"); // end
    gROOT->Reset();

    return result_container;
}

Int_t main(int argc, char** argv) {
    Config& conf = Config::getInstance();

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
    std::vector<Int_t> ana_run_num{ 309, 320, 332, 343, 355, 377, 387, 456, 466, 476, 486, 496, 506, 516, 520, 524, 528, 532 };

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = WORK_DIR + "/data/pedestal.root";
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile fout(output_path.Data(), "create");
    TTree output_tree("tree", ""); 

    // -- prepare root file branch -----
    std::vector<Double_t> bac_ped_pos_val, bac_ped_pos_err, bacsum_ped_pos_val, bacsum_ped_pos_err;
    std::vector<Double_t> sac_ped_pos_val, sac_ped_pos_err, sacsum_ped_pos_val, sacsum_ped_pos_err;
    std::vector<Double_t> kvc_ped_pos_val, kvc_ped_pos_err, kvcsum_ped_pos_val, kvcsum_ped_pos_err;
    std::vector<Double_t> bac_ped_sig_val, bac_ped_sig_err, bacsum_ped_sig_val, bacsum_ped_sig_err;
    std::vector<Double_t> sac_ped_sig_val, sac_ped_sig_err, sacsum_ped_sig_val, sacsum_ped_sig_err;
    std::vector<Double_t> kvc_ped_sig_val, kvc_ped_sig_err, kvcsum_ped_sig_val, kvcsum_ped_sig_err;
    Int_t tmp_run_num;

    output_tree.Branch("run_num", &tmp_run_num, "run_num/I");
    output_tree.Branch("bac_ped_pos_val", &bac_ped_pos_val);
    output_tree.Branch("bac_ped_pos_err", &bac_ped_pos_err);
    output_tree.Branch("bacsum_ped_pos_val", &bacsum_ped_pos_val);
    output_tree.Branch("bacsum_ped_pos_err", &bacsum_ped_pos_err);
    output_tree.Branch("bac_ped_sig_val", &bac_ped_sig_val);
    output_tree.Branch("bac_ped_sig_err", &bac_ped_sig_err);
    output_tree.Branch("bacsum_ped_sig_val", &bacsum_ped_sig_val);
    output_tree.Branch("bacsum_ped_sig_err", &bacsum_ped_sig_err);
    output_tree.Branch("sac_ped_pos_val", &sac_ped_pos_val);
    output_tree.Branch("sac_ped_pos_err", &sac_ped_pos_err);
    output_tree.Branch("sacsum_ped_pos_val", &sacsum_ped_pos_val);
    output_tree.Branch("sacsum_ped_pos_err", &sacsum_ped_pos_err);
    output_tree.Branch("sac_ped_sig_val", &sac_ped_sig_val);
    output_tree.Branch("sac_ped_sig_err", &sac_ped_sig_err);
    output_tree.Branch("sacsum_ped_sig_val", &sacsum_ped_sig_val);
    output_tree.Branch("sacsum_ped_sig_err", &sacsum_ped_sig_err);
    output_tree.Branch("kvc_ped_pos_val", &kvc_ped_pos_val);
    output_tree.Branch("kvc_ped_pos_err", &kvc_ped_pos_err);
    output_tree.Branch("kvcsum_ped_pos_val", &kvcsum_ped_pos_val);
    output_tree.Branch("kvcsum_ped_pos_err", &kvcsum_ped_pos_err);
    output_tree.Branch("kvc_ped_sig_val", &kvc_ped_sig_val);
    output_tree.Branch("kvc_ped_sig_err", &kvc_ped_sig_err);
    output_tree.Branch("kvcsum_ped_sig_val", &kvcsum_ped_sig_val);
    output_tree.Branch("kvcsum_ped_sig_err", &kvcsum_ped_sig_err);

    for (Int_t i = 0, n_run_num = ana_run_num.size(); i < n_run_num; i++) {
        tmp_run_num = ana_run_num[i];
        
        // -- analyze -----
        Int_t pdf_save_mode = 0;
        if (i == 0) pdf_save_mode = 1;
        else if (i == n_run_num-1) pdf_save_mode = 2;
        std::unordered_map<std::string, std::vector<FitResult>> result_container = analyze(ana_run_num[i], pdf_save_mode);

        // -- initialize -----
        bac_ped_pos_val.clear(); bac_ped_pos_err.clear(); bacsum_ped_pos_val.clear(); bacsum_ped_pos_err.clear();
        sac_ped_pos_val.clear(); sac_ped_pos_err.clear(); sacsum_ped_pos_val.clear(); sacsum_ped_pos_err.clear();
        kvc_ped_pos_val.clear(); kvc_ped_pos_err.clear(); kvcsum_ped_pos_val.clear(); kvcsum_ped_pos_err.clear();
        bac_ped_sig_val.clear(); bac_ped_sig_err.clear(); bacsum_ped_sig_val.clear(); bacsum_ped_sig_err.clear();
        sac_ped_sig_val.clear(); sac_ped_sig_err.clear(); sacsum_ped_sig_val.clear(); sacsum_ped_sig_err.clear();
        kvc_ped_sig_val.clear(); kvc_ped_sig_err.clear(); kvcsum_ped_sig_val.clear(); kvcsum_ped_sig_err.clear();

        // -- bac -----
        for (const auto &result : result_container["BAC"]) {
            bac_ped_pos_val.push_back( result.par[1] );
            bac_ped_pos_err.push_back( result.err[1] );
            bac_ped_sig_val.push_back( result.par[2] );
            bac_ped_sig_err.push_back( result.err[2] );
        }
        for (const auto &result : result_container["BACSUM"]) {
            bacsum_ped_pos_val.push_back( result.par[1] );
            bacsum_ped_pos_err.push_back( result.err[1] );
            bacsum_ped_sig_val.push_back( result.par[2] );
            bacsum_ped_sig_err.push_back( result.err[2] );
        }

        // -- sac -----
        for (const auto &result : result_container["SAC"]) {
            sac_ped_pos_val.push_back( result.par[1] );
            sac_ped_pos_err.push_back( result.err[1] );
            sac_ped_sig_val.push_back( result.par[2] );
            sac_ped_sig_err.push_back( result.err[2] );
        }
        for (const auto &result : result_container["SACSUM"]) {
            sacsum_ped_pos_val.push_back( result.par[1] );
            sacsum_ped_pos_err.push_back( result.err[1] );
            sacsum_ped_sig_val.push_back( result.par[2] );
            sacsum_ped_sig_err.push_back( result.err[2] );
        }

        // -- kvc -----
        for (const auto &result : result_container["KVC"]) {
            kvc_ped_pos_val.push_back( result.par[1] );
            kvc_ped_pos_err.push_back( result.err[1] );
            kvc_ped_sig_val.push_back( result.par[2] );
            kvc_ped_sig_err.push_back( result.err[2] );
        }
        for (const auto &result : result_container["KVCSUM"]) {
            kvcsum_ped_pos_val.push_back( result.par[1] );
            kvcsum_ped_pos_err.push_back( result.err[1] );
            kvcsum_ped_sig_val.push_back( result.par[2] );
            kvcsum_ped_sig_err.push_back( result.err[2] );
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