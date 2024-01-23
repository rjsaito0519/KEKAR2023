#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>

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
#include "TTreeReaderArray.h"

#include "include/macro.hh"
#include "include/nagao_macro.hh"
#include "include/variable.hh"


void analyze(Int_t run_num, Bool_t isSave = false) {
    
    // いろんな設定
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.05, "XY");
    gStyle->SetTitleSize(1., "XY");
    gStyle->SetTitleFontSize(0.08);
    gROOT->GetColor(0)->SetAlpha(0.01);
    gStyle->SetPadBottomMargin(0.15);

    // -- load file ------------------------------------------------------------------
    auto *f = new TFile( Form( "../root/kekar_run%05d.root", run_num ) );
    TTreeReader reader("tree", f);
    Int_t tot_num = reader.GetEntries();

//  +------+
//  |  T1  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT1 = 0, MaxAdcT1 = adc_max;
    TTreeReaderValue<Double_t> T1a(reader, "t1a");
    TH1D *hT1a = new TH1D("T1a", Form("run%05d T1(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT1 = TdcGateT1[0], MaxTdcT1 = TdcGateT1[1];
    TTreeReaderArray<Double_t> T1t(reader, "t1t");
    TH1D *hT1t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT1t[n_hit] = new TH1D(Form("T1%d", n_hit+1), Form("run%05d T1(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +------+
//  |  T2  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT2 = 0, MaxAdcT2 = adc_max;
    TTreeReaderValue<Double_t> T2a(reader, "t2a");
    TH1D *hT2a = new TH1D("T2a", Form("run%05d T2(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT2 = TdcGateT2[0], MaxTdcT2 = TdcGateT2[1];
    TTreeReaderArray<Double_t> T2t(reader, "t2t");
    TH1D *hT2t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT2t[n_hit] = new TH1D(Form("T2%d", n_hit+1), Form("run%05d T2(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +------+
//  |  T3  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT3 = 0, MaxAdcT3 = adc_max;
    TTreeReaderValue<Double_t> T3a(reader, "t3a");
    TH1D *hT3a = new TH1D("T3a", Form("run%05d T3(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT3 = TdcGateT3[0], MaxTdcT3 = TdcGateT3[1];
    TTreeReaderArray<Double_t> T3t(reader, "t3t");
    TH1D *hT3t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT3t[n_hit] = new TH1D(Form("T3%d", n_hit+1), Form("run%05d T3(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +------+
//  |  T4  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT4 = 0, MaxAdcT4 = adc_max;
    TTreeReaderValue<Double_t> T4a(reader, "t4a");
    TH1D *hT4a = new TH1D("T4a", Form("run%05d T4(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT4 = TdcGateT4[0], MaxTdcT4 = TdcGateT4[1];
    TTreeReaderArray<Double_t> T4t(reader, "t4t");
    TH1D *hT4t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT4t[n_hit] = new TH1D(Form("T4%d", n_hit+1), Form("run%05d T4(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +-------+
//  |  SAC  |
//  +-------+
// --------------------------------------------------------------------------------------
    // SUMTDC (indivのTDCは取っていない)
    Double_t MinTdcSACSUM = TdcGateSACSUM[0], MaxTdcSACSUM = TdcGateSACSUM[1];
    TTreeReaderArray<Double_t> SACSUMt(reader, "sacsumt");
    TH1D *hSACSUMt[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hSACSUMt[n_hit] = new TH1D(Form("SACSUMt%d", n_hit+1), Form("run%05d SACSUM(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +-------+
//  |  BAC  |
//  +-------+
// --------------------------------------------------------------------------------------    
    // SUMTDC (indivのTDCは取っていない)
    Double_t MinTdcBACSUM = TdcGateBACSUM[0], MaxTdcBACSUM = TdcGateBACSUM[1];
    TTreeReaderArray<Double_t> BACSUMt(reader, "bacsumt");
    TH1D *hBACSUMt[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hBACSUMt[n_hit] = new TH1D(Form("BACSUMt%d", n_hit+1), Form("run%05d BACSUM(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +-------+
//  |  KVC  |
//  +-------+
// --------------------------------------------------------------------------------------    
    // SUMTDC (indivのTDCは取っていない)
    Double_t MinTdcKVCSUM[maxKVCSUMch], MaxTdcKVCSUM[maxKVCSUMch];
    for (Int_t ch = 0; ch < maxKVCSUMch; ch++) {
        MinTdcKVCSUM[ch] = TdcGateKVCSUM[0];
        MaxTdcKVCSUM[ch] = TdcGateKVCSUM[1];
    }
    TTreeReaderArray<Double_t> KVCSUMt(reader, "kvcsumt");
    TH1D *hKVCSUMt[maxKVCSUMch][maxTDChit +1];
    for (Int_t ch = 0; ch < maxKVCSUMch; ch++ ) for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hKVCSUMt[ch][n_hit] = new TH1D(Form("KVCSUMt_ch%d_%d", ch+1, n_hit+1), Form("run%05d KVCSUM ch%d(TDC)", run_num, ch+1), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------


// --- fill event -----------------------------------------------------------------------------------
    reader.Restart();
    while (reader.Next() ){

        // -- trigger ADC ----------
        hT1a->Fill(*T1a);
        hT2a->Fill(*T2a);
        hT3a->Fill(*T3a);
        hT4a->Fill(*T4a);

        // -- TDC ---------
        for (Int_t n_hit = 0; n_hit < maxTDChit; n_hit++) {
            // -- trigger ----------
            hT1t[n_hit]    ->Fill(T1t[n_hit]);
            hT1t[maxTDChit]->Fill(T1t[n_hit]);
            hT2t[n_hit]    ->Fill(T2t[n_hit]);
            hT2t[maxTDChit]->Fill(T2t[n_hit]);
            hT3t[n_hit]    ->Fill(T3t[n_hit]);
            hT3t[maxTDChit]->Fill(T3t[n_hit]);
            hT4t[n_hit]    ->Fill(T4t[n_hit]);
            hT4t[maxTDChit]->Fill(T4t[n_hit]);
            // -- SAC ----------
            if (TdcCutThreshold < SACSUMt[n_hit]) {
                hSACSUMt[n_hit]    ->Fill(SACSUMt[n_hit]);
                hSACSUMt[maxTDChit]->Fill(SACSUMt[n_hit]);
            }
            // -- BAC ----------
            hBACSUMt[n_hit]    ->Fill(BACSUMt[n_hit]);
            hBACSUMt[maxTDChit]->Fill(BACSUMt[n_hit]);
            // -- KVC ----------
            for (Int_t ch = 0; ch < maxKVCSUMch; ch++) {
                hKVCSUMt[ch][n_hit]    ->Fill(KVCSUMt[n_hit + ch*maxTDChit]);
                hKVCSUMt[ch][maxTDChit]->Fill(KVCSUMt[n_hit + ch*maxTDChit]);
            }
        }
    }
// --------------------------------------------------------------------------------------


// --- fitting -----------------------------------------------------------------------------------
    //  +-----------+
    //  |  Trigger  |
    //  +-----------+
    TCanvas *c_trigger_gate = new TCanvas("", "", 1200, 600);
    c_trigger_gate->Divide(4,2);
    Double_t par[3];

    // --- T1 ----------------------------------------------------
    fit_TDC( hT1t[maxTDChit], par, c_trigger_gate, 1 );
    MinTdcT1 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT1 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT1a, par, c_trigger_gate, 5 );
    MinAdcT1 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT1 = par[1] + adc_n_sigma*par[2];
    
    // --- T2 ----------------------------------------------------
    fit_TDC( hT2t[maxTDChit], par, c_trigger_gate, 2 );
    MinTdcT2 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT2 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT2a, par, c_trigger_gate, 6 );
    MinAdcT2 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT2 = par[1] + adc_n_sigma*par[2];
    
    // --- T3 ----------------------------------------------------
    fit_TDC( hT3t[maxTDChit], par, c_trigger_gate, 3 );
    MinTdcT3 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT3 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT3a, par, c_trigger_gate, 7 );
    MinAdcT3 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT3 = par[1] + adc_n_sigma*par[2];
    
    // --- T4 ----------------------------------------------------
    fit_TDC( hT4t[maxTDChit], par, c_trigger_gate, 4 );
    MinTdcT4 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT4 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT4a, par, c_trigger_gate, 8 );
    MinAdcT4 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT4 = par[1] + adc_n_sigma*par[2];
    
    // if (isSave) c_trigger_gate->Print( Form("../img/SAC_run%05d_TriggerandSAC.svg", run_num) );

    //  +---------+
    //  | Counter |
    //  +---------+
    TCanvas *c_counter_gate = new TCanvas("", "", 1200, 600);
    c_counter_gate->Divide(4,2);
    // --- SAC ----------------------------------------------------
    fit_TDC( hSACSUMt[maxTDChit], par, c_counter_gate, 1 );
    MinTdcSACSUM = par[1] - tdc_n_sigma*par[2];
    MaxTdcSACSUM = par[1] + tdc_n_sigma*par[2];

    // --- BAC ----------------------------------------------------
    fit_TDC( hBACSUMt[maxTDChit], par, c_counter_gate, 2 );
    MinTdcBACSUM = par[1] - tdc_n_sigma*par[2];
    MaxTdcBACSUM = par[1] + tdc_n_sigma*par[2];
    
    // --- KVC ----------------------------------------------------
    for (Int_t ch = 0; ch < maxKVCSUMch; ch++) {
        fit_TDC( hKVCSUMt[ch][maxTDChit], par, c_counter_gate, 5+ch );
        MinTdcKVCSUM[ch] = par[1] - tdc_n_sigma*par[2];
        MaxTdcKVCSUM[ch] = par[1] + tdc_n_sigma*par[2];
    }

    // if (isSave) c_trigger_gate->Print( Form("../img/SAC_run%05d_TriggerandSAC.svg", run_num) );
// --------------------------------------------------------------------------------------


// --- check efficiency -----------------------------------------------------------------------------------
    Int_t n_trig = 0, n_SAC = 0, n_BAC = 0, n_KVC = 0;
    Bool_t isHitT1t = false, isHitT2t = false, isHitT3t = false, isHitT4t = false, isHitT1a = false, isHitT2a = false, isHitT3a = false, isHitT4a = false;
    Bool_t isHitSACSUMt = false, isHitBACSUMt = false, isHitKVCSUMt = false;
    reader.Restart();
    while (reader.Next() ){

        // --- make flag ----------------------------------------------------
        isHitT1t = false, isHitT2t = false, isHitT3t = false, isHitT4t = false, isHitT1a = false, isHitT2a = false, isHitT3a = false, isHitT4a = false;
        isHitSACSUMt = false, isHitBACSUMt = false, isHitKVCSUMt = false;
        for (Int_t n_hit = 0; n_hit < maxTDChit; n_hit++) {
            if ( MinTdcT1 < T1t[n_hit] && T1t[n_hit] < MaxTdcT1 ) isHitT1t = true;
            if ( MinTdcT2 < T2t[n_hit] && T2t[n_hit] < MaxTdcT2 ) isHitT2t = true;
            if ( MinTdcT3 < T3t[n_hit] && T3t[n_hit] < MaxTdcT3 ) isHitT3t = true;
            if ( MinTdcT4 < T4t[n_hit] && T4t[n_hit] < MaxTdcT4 ) isHitT4t = true;
            if ( MinAdcT1 < *T1a ) isHitT1a = true;
            if ( MinAdcT2 < *T2a ) isHitT2a = true;
            if ( MinAdcT3 < *T3a ) isHitT3a = true;
            if ( MinAdcT4 < *T4a ) isHitT4a = true;

            if ( MinTdcSACSUM < SACSUMt[n_hit] && SACSUMt[n_hit] < MaxTdcSACSUM ) isHitSACSUMt = true;
            if ( MinTdcBACSUM < BACSUMt[n_hit] && BACSUMt[n_hit] < MaxTdcBACSUM ) isHitBACSUMt = true;
            for (Int_t ch = 0; ch < maxKVCSUMch; ch++ ) if ( MinTdcKVCSUM[ch] < KVCSUMt[n_hit + ch*maxTDChit] && KVCSUMt[n_hit + ch*maxTDChit] < MaxTdcKVCSUM[ch] ) isHitKVCSUMt = true;
        }

        if (isHitT1t && isHitT2t && isHitT3t && isHitT4t && isHitT1a && isHitT2a && isHitT3a && isHitT4a ) {
            n_trig++;
            if (isHitSACSUMt) n_SAC++;
            if (isHitBACSUMt) n_BAC++;
            if (isHitKVCSUMt) n_KVC++;
        }
    }

    std::cout << (Double_t) n_SAC/n_trig << std::endl;
    std::cout << (Double_t) n_BAC/n_trig << std::endl;
    std::cout << (Double_t) n_KVC/n_trig << std::endl;
    
    
// --- Write data -----------------------------------------------------------------------------------
    system(Form("rm tmp_data.csv"));
    std::ofstream ofs("./tmp_data.csv", std::ios::app);
    ofs << "run,tot_num,n_trig,n_SAC,n_BAC,n_KVC,effSAC,effBAC,effKVC\n";
    ofs << run_num << "," << tot_num << "," << n_trig << "," << n_SAC << "," << n_BAC << "," << n_KVC << "," << (Double_t) n_SAC/n_trig << "," << (Double_t) n_BAC/n_trig << "," << (Double_t) n_KVC/n_trig << "\n";
    ofs.close();
    system(Form("python3 write.py efficiency"));
}

void efficiency(Int_t input_run_num){

    analyze( input_run_num );
    
}

void efficiency_all(Int_t init, Int_t end){

    std::vector<Int_t> except_data1{ 309, 311, 320, 321, 326, 332, 337, 343, 350, 355, 356, 358, 359, 364, 367, 369, 370, 375, 377, 387, 389 };
    // 300 - 392
    for (Int_t i = init; i <= end; i++) {
        Bool_t ana_flag = true;
        for (Int_t j = 0; j < except_data1.size(); j++) if (i == except_data1[j]) ana_flag = false;
        if (ana_flag) analyze( i );
    }
    
    // std::vector<Int_t> except_data2{ 456, 458, 466, 476, 486, 496, 506, 516, 520, 524, 528 };
    // // 447 - 531
    // for (Int_t i = init; i <= end; i++) {
    //     Bool_t ana_flag = true;
    //     for (Int_t j = 0; j < except_data2.size(); j++) if (i == except_data2[j]) ana_flag = false;
    //     if (ana_flag) analyze( i );
    // }
    
}

Int_t main(int argc, char** argv) {
    TApplication *theApp = new TApplication("App", &argc, argv);
    if (argc == 2) {
        efficiency(atoi(argv[1]));
    } else {
        efficiency_all(atoi(argv[1]), atoi(argv[2]));
    }
    theApp->Run();

    return 0;
}
