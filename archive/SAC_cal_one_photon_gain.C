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
#include "TGraphErrors.h"

#include "include/macro.hh"
#include "include/nagao_macro.hh"
#include "include/variable.hh"


void analyze(Int_t run_num, Int_t use_ch, Double_t *cal_result, TVirtualPad *c, Int_t n_c) {

    // いろんな設定
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.05, "XY");
    gStyle->SetTitleSize(1., "XY");
    gStyle->SetTitleFontSize(0.08);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // fileの読み込み
    TFile *f = new TFile( Form( "../root/kekar_run%05d.root", run_num ) );
    TTree *t = (TTree*)f->Get("tree");

    // fill用の変数用意
    Double_t SACa[maxSACch];
    t->SetBranchAddress("saca", &SACa);
    
    // ヒストグラムの準備    
    TH1D *hSACa[maxSACch];
    
    for (Int_t ch = 0; ch < maxSACch; ch++) hSACa[ch] = new TH1D(Form("SACa%d", ch+1), Form("run%05d SAC(ADC) ch%d", run_num, ch+1), 4096, 0, 4096);

    //TEventListを生成する
    t->Draw(">>elist");
    TEventList *elist = (TEventList*)gROOT->FindObject("elist");
    //tree->GetEventList()でもいいと思う
    Int_t tot_num = elist->GetN(); //GetNはその条件をみたすEventの数を返す

    // ヒストグラムにfill
    for (Int_t i = 0; i < tot_num; i++ ){
        t->GetEntry(elist->GetEntry(i));
        for (Int_t ch = 0; ch < maxSACch; ch++) hSACa[ch]->Fill(SACa[ch]); 
    }

    // ペデスタルのfit
    std::cout << "----------------------------------------------------" << std::endl;
    
    Double_t peak_pos = hSACa[use_ch]->GetMaximumBin();

    Int_t bin_info[4];
    searchRange(hSACa[use_ch], bin_info, 5, 10000);
    Double_t pedestal_par[3];
    Double_t pedestal_par_err[3];
    TF1 *pedestal_f = new TF1( "pedestal", "gausn", bin_info[1], peak_pos+2 );
    hSACa[use_ch]->Fit(pedestal_f, "0", "", bin_info[1], peak_pos+2  );
    for (int i = 0; i < 3; i++) {
        pedestal_par[i] = pedestal_f->GetParameter(i);
        pedestal_par_err[i] = pedestal_f->GetParError(i);
    }
    std::cout << "----------------------------------------------------" << std::endl;

    // one photon gainの計算
    Double_t lambda_p0, sigma_lambda_p0, lambda_mean, sigma_lambda_mean, ped_pos, sigma_ped_pos, gain, sigma_gain;
    lambda_p0 = -log( pedestal_par[0]/tot_num );
    sigma_lambda_p0 = pedestal_par_err[0] / pedestal_par[0];

    ped_pos = pedestal_par[1];
    sigma_ped_pos = pedestal_par_err[1];

    lambda_mean = hSACa[use_ch]->GetMean();
    sigma_lambda_mean = hSACa[use_ch]->GetMeanError();

    gain = ( lambda_mean - ped_pos ) / lambda_p0;
    sigma_gain = sqrt( pow( sigma_lambda_mean, 2 )/pow( lambda_p0, 2 ) + pow( sigma_ped_pos, 2 )/pow( lambda_p0, 2 ) + pow( lambda_mean - ped_pos, 2 ) * pow( sigma_lambda_p0, 2 ) / pow( lambda_p0, 4 ) );
    cal_result[0] = gain;
    cal_result[1] = sigma_gain;

    // グラフに描画
    c->cd(n_c);
    gPad->SetLogy(1);
    hSACa[use_ch]->GetXaxis()->SetRangeUser(bin_info[0]-15, bin_info[3]);
    hSACa[use_ch]->GetXaxis()->SetNdivisions(508);
    SetTH1(hSACa[use_ch], Form("Run%05d SAC", run_num), "ADC", "");
    hSACa[use_ch]->Draw();
    
    pedestal_f->SetNpx(1000);
    pedestal_f->SetLineColor(kOrange);
    pedestal_f->SetLineWidth( 2 ); // 線の太さ変更
    pedestal_f->Draw("same");

}

void wrapper(Int_t *run_num_list){

    TColor *color = gROOT->GetColor(0);
    color->SetAlpha(0.01);

    TCanvas *c[maxSACch];
    TVirtualPad *c_gain[maxSACch];
    for (Int_t i = 0; i < maxSACch; i++) {
        c[i] = new TCanvas("", Form("LED %f ch%d", (Double_t)run_num_list[0]/10, i+1), 900, 600);
        c[i]->Divide(1, 2);
        c_gain[i] = c[i]->cd(1);
        c_gain[i]->Divide(3, 1, 0.00001, 0.01);
    }

    Double_t x[maxSACch][3] = {
        {2200, 2150, 2100},
        {1881, 1858, 1836},
        {1986, 1957, 1928},
        {2018, 1985, 1951},
        {2082, 2045, 2007},
        {2038, 2003, 1967},
        {1826, 1813, 1799},
        {2109, 2060, 2011}
    };
    Double_t x_err[3] = {0, 0, 0};
    Double_t y[maxSACch][3], y_err[maxSACch][3];
    TGraph *g[maxSACch];

    Double_t cal_result[2];
    for (Int_t ch = 0; ch < maxSACch; ch++) {
        for (Int_t cond = 0; cond < 3; cond++) {
            analyze(run_num_list[cond+1], ch, cal_result, c_gain[ch], cond+1);
            y[ch][cond] = cal_result[0];
            y_err[ch][cond] = cal_result[1];
        }
        g[ch] = new TGraphErrors(3, &x[ch][0], &y[ch][0], &x_err[0], &y_err[ch][0]);
        // g[ch]->SetTitle( Form("LED %.1f V, ch%d", (Double_t)run_num_list[0]/10, ch+1) );
        // g[ch]->GetXaxis()->SetTitle("PMT HV [V]");
        // g[ch]->GetYaxis()->SetTitle("Gain [ch]");
        SetGr(g[ch], Form("LED %.1f V, ch%d", (Double_t)run_num_list[0]/10, ch+1), "PMT HV [V]", "Gain [ch]");
        g[ch]->SetMarkerStyle(24); 
        c[ch]->cd(2);
        g[ch]->Draw("AP");
        c[ch]->Print( Form("../img/gain_LED%d_ch%d.svg", run_num_list[0], ch+1) );
    }

    // gain出力
    std::cout << "-----------------------------------\n" << run_num_list[0] << "\n";
    for (Int_t i = 0; i < maxSACch; i++) {
        std::cout << i+1 << " , " << y[i][0] << " , " << y[i][1] << " , " << y[i][2] << std::endl; 
    }
    std::cout << "-----------------------------------" << std::endl;

    // データ書き出し
    Char_t cmd[100];
    sprintf(cmd, "rm tmp_data.csv");
    system(cmd);
    std::ofstream ofs("./tmp_data.csv", std::ios::app);
    ofs << "ch,1,1,2,2,3,3,LED_vol\n";
    for (Int_t ch = 0; ch < maxSACch; ch++) {
        ofs << ch+1;
        for (Int_t i = 0; i < 3; i++) ofs << "," << y[ch][i] << "," << y_err[ch][i];
        ofs <<"," << (Double_t)run_num_list[0]/10 << "\n"; 
    }
    ofs.close();

    sprintf(cmd, "python3 write.py SAC_LED");
    system(cmd);

}

void SAC_cal_one_photon_gain(){

    //   LED_vol*10, condition1, condition2, condition3
    // Int_t   veto[4] = { 15, 231, 230, 229 };
    Int_t woveto[4] = { 15, 212, 214, 216 };

    // wrapper(veto);
    wrapper(woveto);

}

Int_t main(int argc, char** argv) {
    TApplication *theApp = new TApplication("App", &argc, argv);
    SAC_cal_one_photon_gain();
    theApp->Run();
    
    return 0;
}
