#include <iostream>
#include <sys/stat.h>
#include <fstream>

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

    // fileの読み込み
    TFile *f = new TFile( Form( "../root/kekar_run%05d.root", run_num ) );
    TTree *t = (TTree*)f->Get("tree"); 

//  +------+
//  |  T1  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT1 = 0, MaxAdcT1 = adc_max;
    Double_t T1a[1];
    t->SetBranchAddress("t1a", &T1a);
    TH1D *hT1a = new TH1D("T1a", Form("run%05d T1(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT1 = TdcGateT1[0], MaxTdcT1 = TdcGateT1[1];
    Double_t T1t[1][maxTDChit];
    t->SetBranchAddress("t1t", &T1t);
    TH1D *hT1t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT1t[n_hit] = new TH1D(Form("T1%d", n_hit+1), Form("run%05d T1(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +------+
//  |  T2  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT2 = 0, MaxAdcT2 = adc_max;
    Double_t T2a[1];
    t->SetBranchAddress("t2a", &T2a);
    TH1D *hT2a = new TH1D("T2a", Form("run%05d T2(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT2 = TdcGateT2[0], MaxTdcT2 = TdcGateT2[1];
    Double_t T2t[1][maxTDChit];
    t->SetBranchAddress("t2t", &T2t);
    TH1D *hT2t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT2t[n_hit] = new TH1D(Form("T2%d", n_hit+1), Form("run%05d T2(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +------+
//  |  T3  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT3 = 0, MaxAdcT3 = adc_max;
    Double_t T3a[1];
    t->SetBranchAddress("t3a", &T3a);
    TH1D *hT3a = new TH1D("T3a", Form("run%05d T3(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT3 = TdcGateT3[0], MaxTdcT3 = TdcGateT3[1];
    Double_t T3t[1][maxTDChit];
    t->SetBranchAddress("t3t", &T3t);
    TH1D *hT3t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT3t[n_hit] = new TH1D(Form("T3%d", n_hit+1), Form("run%05d T3(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +------+
//  |  T4  |
//  +------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t MinAdcT4 = 0, MaxAdcT4 = adc_max;
    Double_t T4a[1];
    t->SetBranchAddress("t4a", &T4a);
    TH1D *hT4a = new TH1D("T4a", Form("run%05d T4(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // TDC
    Double_t MinTdcT4 = TdcGateT4[0], MaxTdcT4 = TdcGateT4[1];
    Double_t T4t[1][maxTDChit];
    t->SetBranchAddress("t4t", &T4t);
    TH1D *hT4t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT4t[n_hit] = new TH1D(Form("T4%d", n_hit+1), Form("run%05d T4(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +-------+
//  |  SAC  |
//  +-------+
// --------------------------------------------------------------------------------------
    // SUMADC
    Double_t SACSUMa[1];
    t->SetBranchAddress("sacsuma", &SACSUMa);
    TH1D *hSACSUMa     = new TH1D(    "SACOnlineSUMa", Form("run%05d SAC Online SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // SUMTDC (indivのTDCは取っていない)
    Double_t MinTdcSACSUM = TdcGateSACSUM[0], MaxTdcSACSUM = TdcGateSACSUM[1];
    Double_t SACSUMt[1][maxTDChit];
    t->SetBranchAddress("sacsumt", &SACSUMt);
    TH1D *hSACSUMt[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hSACSUMt[n_hit] = new TH1D(Form("SACSUMt%d", n_hit+1), Form("run%05d SACSUM(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +-------+
//  |  BAC  |
//  +-------+
// --------------------------------------------------------------------------------------
    // SUMADC
    Double_t BACSUMa[1];
    t->SetBranchAddress("bacsuma", &BACSUMa);
    TH1D *hBACSUMa     = new TH1D(    "BACOnlineSUMa", Form("run%05d BAC Online SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // SUMTDC (indivのTDCは取っていない)
    Double_t MinTdcBACSUM = TdcGateBACSUM[0], MaxTdcBACSUM = TdcGateBACSUM[1];
    Double_t BACSUMt[1][maxTDChit];
    t->SetBranchAddress("bacsumt", &BACSUMt);
    TH1D *hBACSUMt[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hBACSUMt[n_hit] = new TH1D(Form("BACSUMt%d", n_hit+1), Form("run%05d BACSUM(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +-------+
//  |  KVC  |
//  +-------+
// --------------------------------------------------------------------------------------
    // SUMADC
    Double_t KVCSUMa[maxKVCSUMch];
    t->SetBranchAddress("kvcsuma", &KVCSUMa);
    TH1D *hKVCSUMa[maxKVCSUMch];
    for (Int_t ch = 0; ch < maxKVCSUMch; ch++ ) hKVCSUMa[ch] = new TH1D( Form("KVCOnlineSUMa_ch%d", ch+1), Form("run%05d KVC Online SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // SUMTDC (indivのTDCは取っていない)
    Double_t MinTdcKVCSUM = TdcGateKVCSUM[0], MaxTdcKVCSUM = TdcGateKVCSUM[1];
    Double_t KVCSUMt[maxKVCSUMch][maxTDChit];
    t->SetBranchAddress("kvcsumt", &KVCSUMt);
    TH1D *hKVCSUMt[maxKVCSUMch][maxTDChit +1];
    for (Int_t ch = 0; ch < maxKVCSUMch; ch++ ) for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hKVCSUMt[ch][n_hit] = new TH1D(Form("KVCSUMt_ch%d_%d", ch+1, n_hit+1), Form("run%05d KVCSUM(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------


// --- fill event -----------------------------------------------------------------------------------
    //TEventListを生成する
    t->Draw(">>elist");
    // t->Draw(">>elist");
    TEventList *elist = (TEventList*)gROOT->FindObject("elist");
    //tree->GetEventList()でもいいと思う
    Int_t tot_num = elist->GetN(); //GetNはその条件をみたすEventの数を返す

    Int_t n_trig = 0, n_SAC = 0, n_BAC = 0, n_KVC = 0;
    Bool_t isHitT1t = false, isHitT2t = false, isHitT3t = false, isHitT4t = false;
    Bool_t isHitSACSUMt = false, isHitBACSUMt = false, isHitKVCSUMt = false;
    for (Int_t i = 0; i < tot_num; i++ ){
        t->GetEntry(elist->GetEntry(i));
        
        // --- set flag ----------------------------------------------------
        isHitT1t = false, isHitT2t = false, isHitT3t = false, isHitT4t = false;
        isHitSACSUMt = false, isHitBACSUMt = false, isHitKVCSUMt = false;
        for (Int_t n_hit = 0; n_hit < maxTDChit; n_hit++) {
            if ( MinTdcT1 < T1t[0][n_hit] && T1t[0][n_hit] < MaxTdcT1 ) isHitT1t = true;
            if ( MinTdcT2 < T2t[0][n_hit] && T2t[0][n_hit] < MaxTdcT2 ) isHitT2t = true;
            if ( MinTdcT3 < T3t[0][n_hit] && T3t[0][n_hit] < MaxTdcT3 ) isHitT3t = true;
            if ( MinTdcT4 < T4t[0][n_hit] && T4t[0][n_hit] < MaxTdcT4 ) isHitT4t = true;

            if ( MinTdcSACSUM < SACSUMt[0][n_hit] && SACSUMt[0][n_hit] < MaxTdcSACSUM ) isHitSACSUMt = true;
            if ( MinTdcBACSUM < BACSUMt[0][n_hit] && BACSUMt[0][n_hit] < MaxTdcBACSUM ) isHitBACSUMt = true;
            for (Int_t ch = 0; ch < maxKVCSUMch; ch++ ) if ( MinTdcKVCSUM < KVCSUMt[ch][n_hit] && KVCSUMt[ch][n_hit] < MaxTdcKVCSUM ) isHitKVCSUMt = true;
        }

        if (isHitT1t && isHitT2t && isHitT3t && isHitT4t ) {
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
    Char_t cmd[100];
    sprintf(cmd, "rm tmp_data.csv");
    system(cmd);
    std::ofstream ofs("./tmp_data.csv", std::ios::app);
    ofs << "run,tot_num,n_trig,n_SAC,n_BAC,n_KVC,effSAC,effBAC,effKVC\n";
    ofs << run_num << "," << tot_num << "," << n_trig << "," << n_SAC << "," << n_BAC << "," << n_KVC << "," << (Double_t) n_SAC/n_trig << "," << (Double_t) n_BAC/n_trig << "," << (Double_t) n_KVC/n_trig << "\n";
    ofs.close();
    sprintf(cmd, "python3 write.py efficiency");
    system(cmd);
}

void efficiency(Int_t input_run_num){

    analyze( input_run_num );
    
}

void efficiency_all(){

    // 300 - 392
    // for (Int_t i = 372; i <= 392; i++) {
    //     analyze( i );
    // }
    
    // 445 - 531
    for (Int_t i = 526; i <= 531; i++) {
        analyze( i );
    }
    
}

Int_t main(int argc, char** argv) {
    TApplication *theApp = new TApplication("App", &argc, argv);
    if (argc == 2) {
        efficiency(atoi(argv[1]));
    } else {
        efficiency_all();
    }
    theApp->Run();

    return 0;
}
