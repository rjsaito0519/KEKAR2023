#include <iostream>
#include <sys/stat.h>
#include <ncurses.h>
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
    gStyle->SetTitleSize(1, "XY");
    gStyle->SetTitleFontSize(0.08);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // fileの読み込み
    TFile *f = new TFile( Form( "../root/kekar_run%05d.root", run_num ) );
    TTree *t = (TTree*)f->Get("tree"); 
    
    // --- SAC ----------------------------------------------------
    Double_t SACa[maxSACch];
    t->SetBranchAddress("saca", &SACa);
    TH1D *hSACa[maxSACch];
    for (Int_t ch = 0; ch < maxSACch; ch++) hSACa[ch] = new TH1D( Form("SACa_%d", ch), Form("run%05d SAC ch%d (ADC)", run_num, ch+1), adc_bin_num, adc_min, adc_max);
    // --------------------------------------------------------------------------------------

    // // --- BAC ----------------------------------------------------
    // Double_t BACa[maxBACch];
    // t->SetBranchAddress("E72BACa", &BACa);
    // TH1D *hBACa[maxBACch];
    // for (Int_t ch = 0; ch < maxBACch; ch++) hBACa[ch] = new TH1D( Form("BACa_%d", ch), Form("run%05d BAC ch%d SUM(ADC)", run_num, ch+1), adc_bin_num, adc_min, adc_max);
    // // --------------------------------------------------------------------------------------

    // // --- KVC ----------------------------------------------------
    // Double_t KVCSUMa[1];
    // t->SetBranchAddress("KVCSUMa", &KVCSUMa);
    // TH1D *hKVCSUMa = new TH1D( "KVCSUMa", Form("run%05d KVC SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    // // --------------------------------------------------------------------------------------


// --- fill event -----------------------------------------------------------------------------------
    //TEventListを生成する
    t->Draw(">>elist");
    TEventList *elist = (TEventList*)gROOT->FindObject("elist");
    //tree->GetEventList()でもいいと思う
    Int_t tot_num = elist->GetN(); //GetNはその条件をみたすEventの数を返す

    // ヒストグラムにfill
    for (Int_t i = 0; i < tot_num; i++ ){
        t->GetEntry(elist->GetEntry(i));
        for (Int_t ch = 0; ch < maxSACch; ch++ ) hSACa[ch]->Fill(SACa[ch]);        

        // hSACSUMa->Fill(SACSUMa[0]);
        // for (Int_t ch = 0; ch < maxBACch; ch++ ) hBACa[ch]->Fill(BACa[ch]);
        // hKVCSUMa->Fill(KVCSUMa[0]);     
    }
// --------------------------------------------------------------------------------------


// --- fitting -----------------------------------------------------------------------------------
    TCanvas *c_pedestal = new TCanvas("", "", 1200, 600);
    c_pedestal->Divide(maxSACch,1);
    Double_t par[6];
    // エラーを受け取る配列は用意しているがmacro側の設定がめんどいので入れてない。どうせ2桁小さいので気にしなくて良い気もする
    Double_t SAC_pedestal_par[maxSACch][3], SAC_pedestal_par_err[maxSACch][3];
    Double_t BAC_pedestal_par[maxBACch][3], BAC_pedestal_par_err[maxBACch][3];
    Double_t KVC_pedestal_par[3], KVC_pedestal_par_err[3];
    
    // --- SAC ----------------------------------------------------
    // fit_pedestal( hSACSUMa, par, c_pedestal, 1, 5, 30 );
    // for (int i = 0; i < 3; i++) SAC_pedestal_par[i]     = par[i];
    // // for (int i = 0; i < 3; i++) SAC_pedestal_par_err[i] = par[i+3];
    for (Int_t ch = 0; ch < maxSACch; ch++) {
        fit_pedestal( hSACa[ch], par, c_pedestal, ch+1, 100, 3 );
        for (int i = 0; i < 3; i++) SAC_pedestal_par[ch][i] = par[i];
        // for (int i = 0; i < 3; i++) BAC_pedestal_par_err[i] = par[i+3];
    }
    

    // // --- BAC ----------------------------------------------------
    // for (Int_t ch = 0; ch < maxBACch; ch++) {
    //     fit_pedestal( hBACa[ch], par, c_pedestal, ch+1, 5, 40 );
    //     for (int i = 0; i < 3; i++) BAC_pedestal_par[ch][i] = par[i];
    //     // for (int i = 0; i < 3; i++) BAC_pedestal_par_err[i] = par[i+3];
    // }
    
    // // --- KVC ----------------------------------------------------
    // fit_pedestal( hKVCSUMa, par, c_pedestal, 3, 5, 10 );
    // for (int i = 0; i < 3; i++) KVC_pedestal_par[i]     = par[i];
    // // for (int i = 0; i < 3; i++) KVC_pedestal_par_err[i] = par[i+3];
    
    // if (isSave) c_pedestal->Print( Form("../img/run%05d_pedestal.svg", run_num) );

// --- Write data -----------------------------------------------------------------------------------
    Char_t cmd[100];
    sprintf(cmd, "rm tmp_data.csv");
    system(cmd);
    
    std::ofstream ofs("./tmp_data.csv", std::ios::app);
    ofs << "run,counter,ch,ped_amp,ped_center,ped_sigma";
    // ofs << run_num << "," << "SAC";
    // for (Int_t i = 0; i < 3; i++) ofs << "," << SAC_pedestal_par[i];
    for (Int_t ch = 0; ch < maxSACch; ch++) {
        ofs << "\n" << run_num << "," << "SAC," << ch+1;
        for (Int_t i = 0; i < 3; i++) ofs << "," << SAC_pedestal_par[ch][i];
    }
    // ofs << "\n"  << run_num << "," << "KVC";
    // for (Int_t i = 0; i < 3; i++) ofs << "," << KVC_pedestal_par[i];
    ofs << "\n";
    ofs.close();

    sprintf(cmd, "python3 write.py pedestal_indiv");
    system(cmd);

}

void pedestal(Int_t input_run_num, Bool_t isSave = false){

    analyze( input_run_num, isSave );
    
}

Int_t main(int argc, char** argv) {
    TApplication *theApp = new TApplication("App", &argc, argv);
    if (argc > 2) {
        pedestal(atoi(argv[1]), argv[2]);
    } else {
        pedestal(atoi(argv[1]));
    }
    theApp->Run();
    
    // initscr();
    // while (true) {
    //     int ch = getch();
    //     if (ch == 'q') break;
    // }
    // endwin();

    return 0;
}