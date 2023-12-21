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
    Double_t MinTdcT1 = 0, MaxTdcT1 = tdc_max;
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
    Double_t MinTdcT2 = 0, MaxTdcT2 = tdc_max;
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
    Double_t MinTdcT3 = 0, MaxTdcT3 = tdc_max;
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
    Double_t MinTdcT4 = 0, MaxTdcT4 = tdc_max;
    Double_t T4t[1][maxTDChit];
    t->SetBranchAddress("t4t", &T4t);
    TH1D *hT4t[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hT4t[n_hit] = new TH1D(Form("T4%d", n_hit+1), Form("run%05d T4(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);
// --------------------------------------------------------------------------------------

//  +-------+
//  |  SAC  |
//  +-------+
// --------------------------------------------------------------------------------------
    // ADC
    Double_t SACa[maxSACch];
    t->SetBranchAddress("saca", &SACa);
    TH1D *hSACa[maxSACch];
    TH1D *hTrigSACa[maxSACch];
    for (Int_t ch = 0; ch < maxSACch; ch++){
        hSACa[ch]     = new TH1D(Form(    "SACa_%d", ch+1), Form("run%05d SAC(ADC) ch%d", run_num, ch+1), adc_bin_num, adc_min, adc_max);
        hTrigSACa[ch] = new TH1D(Form("trigSACa_%d", ch+1), Form("run%05d SAC(ADC) ch%d", run_num, ch+1), adc_bin_num, adc_min, adc_max);
    }

    // NPE
    TH1D *hSACnpe[maxSACch];
    TH1D *hTrigSACnpe[maxSACch];
    for (Int_t ch = 0; ch < maxSACch; ch++){
        hSACnpe[ch]     = new TH1D(Form(    "SACnpe_%d", ch+1), Form("run%05d SAC(NPE) ch%d", run_num, ch+1), npe_bin_num, npe_min, npe_max);
        hTrigSACnpe[ch] = new TH1D(Form("trigSACnpe_%d", ch+1), Form("run%05d SAC(NPE) ch%d", run_num, ch+1), npe_bin_num, npe_min, npe_max);
    }

    // SUMADC
    Double_t SACSUMa[1];
    t->SetBranchAddress("sacsuma", &SACSUMa);
    TH1D *hSACOnlineSUMa     = new TH1D(    "SACOnlineSUMa", Form("run%05d SAC Online SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    TH1D *hTrigSACOnlineSUMa = new TH1D("trigSACOnlineSUMa", Form("run%05d SAC Online SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    TH1D *hSACOfflineSUMa     = new TH1D(    "SACOfflineSUMa", Form("run%05d SAC Offline SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    TH1D *hTrigSACOfflineSUMa = new TH1D("trigSACOfflineSUMa", Form("run%05d SAC Offline SUM(ADC)", run_num), adc_bin_num, adc_min, adc_max);
    
    // SUMNPE
    TH1D *hSACOnlineSUMnpe     = new TH1D(    "SACOnlineSUMnpe", Form("run%05d SAC Online SUM(NPE)", run_num), npe_bin_num, npe_min, npe_max);
    TH1D *hTrigSACOnlineSUMnpe = new TH1D("trigSACOnlineSUMnpe", Form("run%05d SAC Online SUM(NPE)", run_num), npe_bin_num, npe_min, npe_max);
    TH1D *hSACOfflineSUMnpe     = new TH1D(    "SACOfflineSUMNPE", Form("run%05d SAC Offline SUM(NPE)", run_num), npe_bin_num, npe_min, npe_max);
    TH1D *hTrigSACOfflineSUMnpe = new TH1D("trigSACOfflineSUMNPE", Form("run%05d SAC Offline SUM(NPE)", run_num), npe_bin_num, npe_min, npe_max);

    // SUMTDC (indivのTDCは取っていない)
    Double_t MinTdcSACSUM = 0., MaxTdcSACSUM = tdc_max;
    Double_t SACSUMt[1][maxTDChit];
    t->SetBranchAddress("sacsumt", &SACSUMt);
    TH1D *hSACSUMt[maxTDChit +1];
    for (Int_t n_hit = 0; n_hit < maxTDChit+1; n_hit++) hSACSUMt[n_hit] = new TH1D(Form("SACSUMt%d", n_hit+1), Form("run%05d SACSUM(TDC)", run_num), tdc_bin_num, tdc_min, tdc_max);

    // 2d histogram
    TH2D *hSACSUMcorrelation = new TH2D("SACSUMcorrelation","", npe_bin_num, npe_min, npe_max, adc_bin_num*5, adc_min, adc_max);
    TH2D *hSACSUMcorrelation2 = new TH2D("SACSUMcorrelation2","", npe_bin_num, npe_min, npe_max, npe_bin_num, npe_min, npe_max);

    // divided histogram for checking threshold
    TH1D *hSACdivided[2];
    hSACdivided[0] = new TH1D("dividedOnline", Form("run%05d", run_num), npe_bin_num, npe_min, npe_max);
    hSACdivided[1] = new TH1D("dividedOffline", Form("run%05d", run_num), npe_bin_num, npe_min, npe_max);
// --------------------------------------------------------------------------------------


// --- fill event -----------------------------------------------------------------------------------
    //TEventListを生成する
    t->Draw(">>elist");
    // t->Draw(">>elist");
    TEventList *elist = (TEventList*)gROOT->FindObject("elist");
    //tree->GetEventList()でもいいと思う
    Int_t tot_num = elist->GetN(); //GetNはその条件をみたすEventの数を返す

    // ヒストグラムにfill
    for (Int_t i = 0; i < tot_num; i++ ){
        t->GetEntry(elist->GetEntry(i));

        hT1a->Fill(T1a[0]);
        hT2a->Fill(T2a[0]);
        hT3a->Fill(T3a[0]);
        hT4a->Fill(T4a[0]);
        for (Int_t n_hit = 0; n_hit < maxTDChit; n_hit++) {
            hT1t[n_hit]    ->Fill(T1t[0][n_hit]);
            hT1t[maxTDChit]->Fill(T1t[0][n_hit]);

            hT2t[n_hit]    ->Fill(T2t[0][n_hit]);
            hT2t[maxTDChit]->Fill(T2t[0][n_hit]);

            hT3t[n_hit]    ->Fill(T3t[0][n_hit]);
            hT3t[maxTDChit]->Fill(T3t[0][n_hit]);

            hT4t[n_hit]    ->Fill(T4t[0][n_hit]);
            hT4t[maxTDChit]->Fill(T4t[0][n_hit]);

            hSACSUMt[n_hit]    ->Fill(SACSUMt[0][n_hit]);
            hSACSUMt[maxTDChit]->Fill(SACSUMt[0][n_hit]);
        }
        for (Int_t ch = 0; ch < maxSACch; ch++) hSACa[ch]->Fill(SACa[ch]);
        hSACOnlineSUMa->Fill(SACSUMa[0]);
    }
// --------------------------------------------------------------------------------------


// --- fitting -----------------------------------------------------------------------------------
    // TDCのfitting及びヒスト表示
    TCanvas *c_trigger_gate = new TCanvas("", "", 1200, 600);
    c_trigger_gate->Divide(5,2);
    Double_t par[3];

    // --- T1 ----------------------------------------------------
    fit_TDC( hT1t[maxTDChit], par, c_trigger_gate, 1 );
    MinTdcT1 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT1 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT1a, par, c_trigger_gate, 6 );
    MinAdcT1 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT1 = par[1] + adc_n_sigma*par[2];
    
    // --- T2 ----------------------------------------------------
    fit_TDC( hT2t[maxTDChit], par, c_trigger_gate, 2 );
    MinTdcT2 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT2 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT2a, par, c_trigger_gate, 7 );
    MinAdcT2 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT2 = par[1] + adc_n_sigma*par[2];
    
    // --- T3 ----------------------------------------------------
    fit_TDC( hT3t[maxTDChit], par, c_trigger_gate, 3 );
    MinTdcT3 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT3 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT3a, par, c_trigger_gate, 8 );
    MinAdcT3 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT3 = par[1] + adc_n_sigma*par[2];
    
    // --- T4 ----------------------------------------------------
    fit_TDC( hT4t[maxTDChit], par, c_trigger_gate, 4 );
    MinTdcT4 = par[1] - tdc_n_sigma*par[2];
    MaxTdcT4 = par[1] + tdc_n_sigma*par[2];
    fit_ADC( hT4a, par, c_trigger_gate, 9 );
    MinAdcT4 = par[1] - adc_n_sigma*par[2];
    // MaxAdcT4 = par[1] + adc_n_sigma*par[2];
    
    // --- SAC ----------------------------------------------------
    fit_TDC( hSACSUMt[maxTDChit], par, c_trigger_gate, 5 );
    MinTdcSACSUM = par[1] - tdc_n_sigma*par[2];
    MaxTdcSACSUM = par[1] + tdc_n_sigma*par[2];

    if (isSave) c_trigger_gate->Print( Form("../img/SAC_run%05d_TriggerandSAC.svg", run_num) );

    // ペデスタルのfit
    TCanvas *c_pedestal = new TCanvas("", "", 1200, 600);
    c_pedestal->Divide(4, 2);
    Double_t pedestal_indiv_par[maxSACch][3];
    for (Int_t ch = 0; ch < maxSACch; ch++) {
        fit_pedestal( hSACa[ch], par, c_pedestal, ch+1 );
        for (int i = 0; i < 3; i++) pedestal_indiv_par[ch][i] = par[i];
    }
    if (isSave) c_pedestal->Print( Form("../img/SAC_run%05d_pedestal.svg", run_num) );
// --------------------------------------------------------------------------------------


// --- Fill trigged event -----------------------------------------------------------------------------------
    // TDCで制限かけてヒストグラムにfill ついでにN_p.e.のヒストもfill
    Double_t tmp_sum_a = 0., tmp_sum_npe = 0.;
    Bool_t isHitT1t = false, isHitT2t = false, isHitT3t = false, isHitT4t = false, isHitSACSUM = false, isHitT1a = false, isHitT2a = false, isHitT3a = false, isHitT4a = false;
    Int_t n_trig = 0, n_SAC = 0;
    Double_t sum_pedestal_pos = SACsum_pedestal_pos[ pedestal_index(run_num) ];
    for (Int_t i = 0; i < tot_num; i++ ){
        t->GetEntry(elist->GetEntry(i));
        
        // --- set up flag ----------------------------------------------------
        isHitT1t = false, isHitT2t = false, isHitT3t = false, isHitT4t = false, isHitSACSUM = false, isHitT1a = false, isHitT2a = false, isHitT3a = false, isHitT4a = false;
        if ( MinAdcT1 < T1a[0] && T1a[0] < MaxAdcT1 ) isHitT1a = true;
        if ( MinAdcT2 < T2a[0] && T1a[0] < MaxAdcT2 ) isHitT2a = true;
        if ( MinAdcT3 < T3a[0] && T1a[0] < MaxAdcT3 ) isHitT3a = true;
        if ( MinAdcT4 < T4a[0] && T1a[0] < MaxAdcT4 ) isHitT4a = true;
        for (Int_t n_hit = 0; n_hit < maxTDChit; n_hit++) {
            if ( MinTdcT1 < T1t[0][n_hit] && T1t[0][n_hit] < MaxTdcT1 ) isHitT1t = true;
            if ( MinTdcT2 < T2t[0][n_hit] && T2t[0][n_hit] < MaxTdcT2 ) isHitT2t = true;
            if ( MinTdcT3 < T3t[0][n_hit] && T3t[0][n_hit] < MaxTdcT3 ) isHitT3t = true;
            if ( MinTdcT4 < T4t[0][n_hit] && T4t[0][n_hit] < MaxTdcT4 ) isHitT4t = true;

            if ( MinTdcSACSUM < SACSUMt[0][n_hit] && SACSUMt[0][n_hit] < MaxTdcSACSUM ) isHitSACSUM = true;
        }

        if (isHitT1t && isHitT2t && isHitT3t && isHitT4t && isHitT1a && isHitT2a && isHitT3a && isHitT4a ) {
            n_trig++;
            // --- only trigger ----------------------------------------------------
            tmp_sum_a = 0.;
            tmp_sum_npe = 0.;
            for (Int_t ch = 0; ch < maxSACch; ch++) {
                tmp_sum_a   += SACa[ch];
                tmp_sum_npe += (SACa[ch] - pedestal_indiv_par[ch][1]) / SAC_one_photon_gain[ch][2];
                hSACnpe[ch] -> Fill( (SACa[ch] - pedestal_indiv_par[ch][1]) / SAC_one_photon_gain[ch][2] );
            } 
            hSACOfflineSUMa->Fill( tmp_sum_a ); 
            hSACOfflineSUMnpe->Fill( tmp_sum_npe );

            // --- with SACSUM ----------------------------------------------------
            if (isHitSACSUM) {
                n_SAC++;
                tmp_sum_a = 0.;
                tmp_sum_npe = 0.;
                for (Int_t ch = 0; ch < maxSACch; ch++) {
                    hTrigSACa[ch]         ->Fill(  SACa[ch] );
                    hTrigSACnpe[ch]       ->Fill( (SACa[ch] - pedestal_indiv_par[ch][1]) / SAC_one_photon_gain[ch][2] );
                    tmp_sum_a   += SACa[ch];
                    tmp_sum_npe += (SACa[ch] - pedestal_indiv_par[ch][1]) / SAC_one_photon_gain[ch][2];
                }
                hTrigSACOfflineSUMa->Fill(tmp_sum_a);
                hTrigSACOfflineSUMnpe->Fill( tmp_sum_npe );
                hTrigSACOnlineSUMa->Fill( SACSUMa[0] );
                // 最初にOnline sumのone photon gainをoffline sumとの相関より推測する
                hSACSUMcorrelation->Fill( tmp_sum_npe, SACSUMa[0] - sum_pedestal_pos );
            }
        }
    }

    // one photon gainの推測
    TCanvas *c_correlation = new TCanvas("", "", 800, 800);
    c_correlation->cd(1);
    TF1 *linear_f = new TF1("linear","[0]*x+[1]", 0, 50);
    Double_t par_linear[2];
    linear_f->SetParameters(1., 0.);
    linear_f->SetLineColor(kRed);
    hSACSUMcorrelation->Fit(linear_f, "", "", 0, 50);
    for (Int_t i = 0; i < 2; i++) par_linear[i] = linear_f->GetParameter(i);
    SetTH2(hSACSUMcorrelation, "", "Offline SUM (NPE)", "Online SUM (ADC)");
    hSACSUMcorrelation->GetXaxis()->SetRangeUser(-5, 50);
    hSACSUMcorrelation->Draw("colz");
    linear_f->Draw("same");
    if (isSave) c_correlation->Print( Form("../img/SAC_run%05d_correlation.svg", run_num) );

    for (Int_t i = 0; i < tot_num; i++ ){
        t->GetEntry(elist->GetEntry(i));
        
        // --- set up flag ----------------------------------------------------
        isHitT1t = false, isHitT2t = false, isHitT3t = false, isHitT4t = false, isHitSACSUM = false, isHitT1a = false, isHitT2a = false, isHitT3a = false, isHitT4a = false;
        if ( MinAdcT1 < T1a[0] && T1a[0] < MaxAdcT1 ) isHitT1a = true;
        if ( MinAdcT2 < T2a[0] && T1a[0] < MaxAdcT2 ) isHitT2a = true;
        if ( MinAdcT3 < T3a[0] && T1a[0] < MaxAdcT3 ) isHitT3a = true;
        if ( MinAdcT4 < T4a[0] && T1a[0] < MaxAdcT4 ) isHitT4a = true;
        for (Int_t n_hit = 0; n_hit < maxTDChit; n_hit++) {
            if ( MinTdcT1 < T1t[0][n_hit] && T1t[0][n_hit] < MaxTdcT1 ) isHitT1t = true;
            if ( MinTdcT2 < T2t[0][n_hit] && T2t[0][n_hit] < MaxTdcT2 ) isHitT2t = true;
            if ( MinTdcT3 < T3t[0][n_hit] && T3t[0][n_hit] < MaxTdcT3 ) isHitT3t = true;
            if ( MinTdcT4 < T4t[0][n_hit] && T4t[0][n_hit] < MaxTdcT4 ) isHitT4t = true;

            if ( MinTdcSACSUM < SACSUMt[0][n_hit] && SACSUMt[0][n_hit] < MaxTdcSACSUM ) isHitSACSUM = true;
        }

        if (isHitT1t && isHitT2t && isHitT3t && isHitT4t && isHitT1a && isHitT2a && isHitT3a && isHitT4a ) {
            // --- only trigger ----------------------------------------------------
            hSACOnlineSUMnpe->Fill( (SACSUMa[0] - sum_pedestal_pos - par_linear[1]) / par_linear[0] );            
            // --- with SACSUM ----------------------------------------------------
            if (isHitSACSUM) {
                tmp_sum_npe = 0.;
                for (Int_t ch = 0; ch < maxSACch; ch++) tmp_sum_npe += (SACa[ch] - pedestal_indiv_par[ch][1]) / SAC_one_photon_gain[ch][2];
                hTrigSACOnlineSUMnpe->Fill( (SACSUMa[0] - sum_pedestal_pos - par_linear[1]) / par_linear[0] );
                hSACSUMcorrelation2->Fill( tmp_sum_npe, (SACSUMa[0] - sum_pedestal_pos - par_linear[1]) / par_linear[0]  );
            }
        }
    }

// --- Draw figure -----------------------------------------------------------------------------------

    //  +--------------+
    //  |  individual  |
    //  +--------------+
    Int_t bin_info[4];
    Double_t par_indiv_poisson[maxSACch][2], par_indiv_poisson_err[maxSACch][2];
    TF1 *indiv_poisson = new TF1("indiv_poisson", "[0] * TMath::Poisson(x, [1])", 2, 30);    

    TCanvas *c_SAC = new TCanvas("", "", 1200, 800);
    c_SAC->Divide(maxSACch,2);
    for (Int_t ch = 0; ch < maxSACch; ch++) {
        c_SAC->cd(ch+1);
        gPad->SetLogy(1);
        searchRange(hTrigSACa[ch], bin_info, 5, 1000);
        hTrigSACa[ch]->GetXaxis()->SetRangeUser(bin_info[0]-15, bin_info[2]+15);
        SetTH1(hTrigSACa[ch], Form("run%05d SAC ch%d", run_num, ch+1), "ADC", "");
        hTrigSACa[ch]->Draw();

        c_SAC->cd(ch+1+maxSACch);
        gPad->SetLogy(1);
        SetTH1(hTrigSACnpe[ch], Form("NPE %d", ch+1), "N_{p.e.}", "");
        hTrigSACnpe[ch]->GetXaxis()->SetRangeUser(-5, 30);
        hTrigSACnpe[ch]->Draw();    

        indiv_poisson->SetParameter(0, 2000);
        indiv_poisson->SetParameter(1, hTrigSACnpe[ch]->GetMean());
        indiv_poisson->SetLineColor(kOrange);
        hTrigSACnpe[ch]->Fit(indiv_poisson, "", "", 2, 30);
        for (Int_t i = 0; i < 2; i++) {
            par_indiv_poisson[ch][i]     = indiv_poisson->GetParameter(i);
            par_indiv_poisson_err[ch][i] = indiv_poisson->GetParError(i);
        }
    }
    if (isSave) c_SAC->Print( Form("../img/SAC_run%05d_SAC.svg", run_num) );

    //  +-------+
    //  |  SUM  |
    //  +-------+
    Double_t par_sum_poisson[2][2], par_sum_poisson_err[2][2];
    TF1 *sum_poisson[2];
    for (Int_t i = 0; i < 2; i++) {
        sum_poisson[i] = new TF1(Form("sum_poisson%d", i), "[0] * TMath::Poisson(x, [1])", 2, 8);
        sum_poisson[i]->SetParameter(0, 2000);
        sum_poisson[i]->SetLineColor(kBlue);
        sum_poisson[i]->SetLineStyle( 2 );
    }

    TCanvas *c_SACSUM = new TCanvas("", "", 900, 600);
    c_SACSUM->Divide(2, 1);
    c_SACSUM->cd(1);
    SetTH1(hSACOnlineSUMnpe, Form("run%05d SAC Online SUM", run_num), "NPE", "");
    hSACOnlineSUMnpe->GetXaxis()->SetRangeUser(-5, 50);
    hSACOnlineSUMnpe->Draw(); 
    hTrigSACOnlineSUMnpe->SetLineColor(kRed);
    hTrigSACOnlineSUMnpe->Draw("same");

    sum_poisson[0]->SetParameter(1, hTrigSACOnlineSUMnpe->GetMean());
    hTrigSACOnlineSUMnpe->Fit(sum_poisson[0], "0", "", 2, 8);
    for (Int_t i = 0; i < 2; i++) {
        par_sum_poisson[0][i]     = sum_poisson[0]->GetParameter(i);
        par_sum_poisson_err[0][i] = sum_poisson[0]->GetParError(i);
    }

    c_SACSUM->cd(2);
    SetTH1(hSACOfflineSUMnpe, Form("run%05d SAC Offline SUM", run_num), "NPE", "");
    hSACOfflineSUMnpe->GetXaxis()->SetRangeUser(-5, 50);    
    hSACOfflineSUMnpe->Draw();  
    hTrigSACOfflineSUMnpe->SetLineColor(kRed);    
    hTrigSACOfflineSUMnpe->Draw("same");

    sum_poisson[1]->SetParameter(1, hTrigSACOfflineSUMnpe->GetMean());
    hTrigSACOfflineSUMnpe->Fit(sum_poisson[1], "0", "", 2, 8);
    for (Int_t i = 0; i < 2; i++) {
        par_sum_poisson[1][i]     = sum_poisson[1]->GetParameter(i);
        par_sum_poisson_err[1][i] = sum_poisson[1]->GetParError(i);
    }
    if (isSave) c_SACSUM->Print( Form("../img/SAC_run%05d_SACSUM.svg", run_num) );


    // ----------------------------
    TCanvas *c_SACSUM_test = new TCanvas("", "", 600, 600);

    c_SACSUM_test->cd(1);
    // gPad->SetLogy(1);
    SetTH1(hSACOfflineSUMnpe, "", "N_{p.e.}", "");
    hSACOfflineSUMnpe->GetXaxis()->SetRangeUser(-5, 80);    
    hSACOfflineSUMnpe->Draw();  
    // hTrigSACOfflineSUMnpe->SetFillColor(kRed);
    // hTrigSACOfflineSUMnpe->SetFillStyle(1001);    
    hTrigSACOfflineSUMnpe->Draw("same");

    if (isSave) c_SACSUM_test->Print( Form("../img/SAC_run%05d_test.svg", run_num) );

    //  +---------------+
    //  |  correlation  |
    //  +---------------+
    TCanvas *c_correlation2 = new TCanvas("", "", 800, 800);
    c_correlation2->cd(1);
    SetTH2(hSACSUMcorrelation2, "", "Offline SUM (NPE)", "Online SUM (NPE)");
    hSACSUMcorrelation2->GetXaxis()->SetRangeUser(-2, 30);
    hSACSUMcorrelation2->GetYaxis()->SetRangeUser(-2, 30);
    hSACSUMcorrelation2->Draw("colz");
    TF1 *XequalY = new TF1("XequalY", "x", -1, 25);
    linear_f->SetLineColor(kRed);
    XequalY->Draw("same");
    if (isSave) c_correlation->Print( Form("../img/SAC_run%05d_correlation2.svg", run_num) );

    //  +-------------+
    //  |  threshold  |
    //  +-------------+
    Double_t par_threshold[2][2], par_threshold_err[2][2];
    TF1 *threshold_f[2];
    for (Int_t i = 0; i < 2; i++) {
        threshold_f[i] = new TF1(Form("sum_threshold%d", i), "1 / ( 1+ TMath::Exp( -[0]*(x-[1]) ) )", -5, 25);
        threshold_f[i]->SetParameter(0, 0.1);
        threshold_f[i] ->SetParameter(1, 1.);
        threshold_f[i] ->SetLineColor(kOrange);
    }
    auto vline = new TGraph(2);
    vline->SetLineColor(kBlue);
    vline->SetLineWidth( 2 );
    vline->SetLineStyle( 2 );
    
    TCanvas *c_threshold = new TCanvas("", "", 1200, 600);
    c_threshold->Divide(2, 1);

    c_threshold->cd(1);
    hSACdivided[0]->Divide( hTrigSACOnlineSUMnpe, hSACOnlineSUMnpe, 1, 1 );
    SetTH1(hSACdivided[0], Form("run%05d SAC Online SUM", run_num), "NPE", "");
    hSACdivided[0]->GetXaxis()->SetRangeUser(-5, 25);    
    hSACdivided[0]->Draw();
    hSACdivided[0]->Fit(threshold_f[0], "", "", -5, 25);
    for (Int_t i = 0; i < 2; i++) {
        par_threshold[0][i]     = threshold_f[0]->GetParameter(i);
        par_threshold_err[0][i] = threshold_f[0]->GetParError(i);
    }
    vline->SetPoint(0, par_threshold[0][1], 0);
	vline->SetPoint(1, par_threshold[0][1], 1);
    vline->Draw("same");

    c_threshold->cd(2);
    hSACdivided[1]->Divide( hTrigSACOfflineSUMnpe, hSACOfflineSUMnpe, 1, 1 );
    SetTH1(hSACdivided[1], Form("run%05d SAC Offline SUM", run_num), "NPE", "");
    hSACdivided[1]->GetXaxis()->SetRangeUser(-5, 25);    
    hSACdivided[1]->Draw();
    hSACdivided[1]->Fit(threshold_f[1], "", "", -5, 25);
    for (Int_t i = 0; i < 2; i++) {
        par_threshold[1][i]     = threshold_f[1]->GetParameter(i);
        par_threshold_err[1][i] = threshold_f[1]->GetParError(i);
    }
    vline->SetPoint(0, par_threshold[1][1], 0);
	vline->SetPoint(1, par_threshold[1][1], 1);
    vline->Draw("same");

    if (isSave) c_threshold->Print( Form("../img/SAC_run%05d_threshold.svg", run_num) );
 
// --- Write data -----------------------------------------------------------------------------------
    Char_t cmd[100];
    sprintf(cmd, "rm tmp_data.csv");
    system(cmd);
    std::ofstream ofs("./tmp_data.csv", std::ios::app);
    ofs << "run,tot_num,w/osac,w/sac,linear_a,linear_b,OnSumA,OnSumAerr,OnSumLamb,OnSumLamberr,OffSumA,OffSumAErr,OffSumLamb,OffSumLambErr,OnSumAlpha,OnSumAlphaErr,OnSumN,OnSumNErr,OffSumAlpha,OffSumAlphaErr,OffSumN,OffSumNErr\n";
    ofs << run_num << "," << tot_num << "," << n_trig << "," << n_SAC << "," << par_linear[0] << "," << par_linear[1];
    for (Int_t i = 0; i < 2; i++) ofs << "," << par_sum_poisson[i][0] << "," << par_sum_poisson_err[i][0] << "," << par_sum_poisson[i][1] << "," << par_sum_poisson_err[i][1];
    for (Int_t i = 0; i < 2; i++) ofs << "," << par_threshold[i][0] << "," << par_threshold_err[i][0] << "," << par_threshold[i][1] << "," << par_threshold_err[i][1];
    ofs << "\n";
    ofs.close();
    sprintf(cmd, "python3 write.py SAC_summary");
    system(cmd);
}

void SAC_analysis(Int_t input_run_num, Bool_t isSave = false){

    analyze( input_run_num, isSave );
    
}

void SAC_analysis_all(){

    // threshold
    analyze( 39, true );
    analyze( 40, true );
    analyze( 41, true );
    analyze( 42, true );
    analyze( 43, true );
    analyze( 45, true );

    for (Int_t i = 59; i <= 78; i++) {
        analyze( i, true );
    }
    
}

Int_t main(int argc, char** argv) {
    TApplication *theApp = new TApplication("App", &argc, argv);
    if (argc > 2) {
        SAC_analysis(atoi(argv[1]), argv[2]);
    } else if (argc == 2) {
        SAC_analysis(atoi(argv[1]));
    } else {
        SAC_analysis_all();
    }
    theApp->Run();

    return 0;
}
