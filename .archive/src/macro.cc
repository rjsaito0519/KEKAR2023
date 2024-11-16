#include <iostream>

#include "TMath.h"
#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TColor.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"

#include "../include/variable.hh"

void searchRange(TH1D *h, Int_t *bin_info, Int_t n = 5, Double_t n_threshold = 1000) {
    // 描画範囲の調整
    Int_t max_bin_pos = h->GetMaximumBin(), min_bin_pos = h->GetMaximumBin();
    Double_t tmp_value = 10000000;
    Double_t buf[n];
    Double_t threshold = h->GetBinContent( h->GetMaximumBin() )/n_threshold;
    // 最大値を調べる
    while (tmp_value > threshold) {
        max_bin_pos += n;
        for (Int_t i = 0; i < n; i++) buf[i] = h->GetBinContent(max_bin_pos+i);
        tmp_value = *std::max_element(buf, buf + n);
    }
    Int_t max_fit_range = max_bin_pos;
    while (tmp_value > threshold/10) {
        max_bin_pos += n;
        for (Int_t i = 0; i < n; i++) buf[i] = h->GetBinContent(max_bin_pos+i);
        tmp_value = *std::max_element(buf, buf + n);
    }
    // 最小値を調べる
    tmp_value = 10000000;
    while (tmp_value > threshold) {
        min_bin_pos -= n;
        for (Int_t i = 0; i < n; i++) buf[i] = h->GetBinContent(min_bin_pos-i);
        tmp_value = *std::max_element(buf, buf + n);
    }
    Int_t min_fit_range = min_bin_pos;
    while (tmp_value > threshold/10) {
        min_bin_pos -= 5;
        for (Int_t i = 0; i < n; i++) buf[i] = h->GetBinContent(min_bin_pos-i);
        tmp_value = *std::max_element(buf, buf + n);
    }
    bin_info[0] = min_bin_pos;
    bin_info[1] = min_fit_range;
    bin_info[2] = max_bin_pos;
    bin_info[3] = max_fit_range;
}

void fit_TDC(TH1D *h, Double_t *par, TCanvas *c, Int_t n_c) {
    Int_t bin_info[4];
    searchRange(h, bin_info, 5, 10);
    for (Int_t i = 0; i < 4; i++) bin_info[i] *= tdc_ticks_size;
    TF1 *fit_f = new TF1( "gauss", "gausn", bin_info[1], bin_info[3]  );
    fit_f->SetParameter(1, h->GetMaximumBin()*tdc_ticks_size);
    fit_f->SetParameter(2, 5);
    std::cout << h->GetMaximumBin()*tdc_ticks_size << "," << bin_info[1] << ", " << bin_info[3] <<  "\n-------------------" << std::endl;
    TString fit_option = '0';
    if ( h->GetBinContent( h->GetMaximumBin() ) < 100 ) fit_option += 'L';
    h->Fit(fit_f, fit_option.Data(), "", bin_info[1], bin_info[3] );
    for (Int_t i = 0; i < 3; i++) par[i] = fit_f->GetParameter(i);
    c->cd(n_c);
    h->GetXaxis()->SetRangeUser(par[1] - (tdc_n_sigma+2)*par[2], par[1] + (tdc_n_sigma+2)*par[2]);
    h->Draw();
    fit_f->SetNpx(1000);
    fit_f->SetLineColor(kOrange);
    fit_f->SetLineWidth( 2 ); // 線の太さ変更
    fit_f->Draw("same");

    auto gr = new TGraph(5);
	Double_t x1 = par[1] - tdc_n_sigma*par[2];
	Double_t x2 = par[1] + tdc_n_sigma*par[2];
	Double_t y1 = 0;
	Double_t y2 = h->GetBinContent( h->GetMaximumBin() );
	gr->SetPoint(0,x1, y1);
	gr->SetPoint(1,x2, y1);
	gr->SetPoint(2,x2, y2);
	gr->SetPoint(3,x1, y2);
	gr->SetPoint(4,x1, y1);
	gr->Draw("F");
	gr->SetFillColor(kBlue);
    gr->SetFillStyle(3353);
	c->Update();

    h->GetXaxis()->SetNdivisions(505);
}

void fit_ADC(TH1D *h, Double_t *par, TCanvas *c, Int_t n_c, Int_t peak_index = 0) {
    c->cd(n_c);

    TSpectrum *s = new TSpectrum(2);
    s->Search(h, 5, "new", 0.005);
    Double_t *peak_pos = s->GetPositionX();

    // 三角マーカーを消す
    TList *functions = h->GetListOfFunctions();
    TPolyMarker *pm = ( TPolyMarker *)functions-> FindObject ( "TPolyMarker" );
    functions->Remove(pm);
    delete pm;

    TF1 *fit_f = new TF1( "landau", "landaun", peak_pos[peak_index]-25, peak_pos[peak_index]+75  );
    fit_f->SetParameter(1, peak_pos[peak_index]);
    fit_f->SetParameter(2, 1.);
    h->Fit(fit_f, "0", "",  peak_pos[peak_index]-25, peak_pos[peak_index]+75 );
    for (int i = 0; i < 3; i++) par[i] = fit_f->GetParameter(i);

    gPad->SetLogy(1);
    h->GetXaxis()->SetRangeUser(par[1] - adc_n_sigma*par[2] - 50, peak_pos[peak_index]+200);
    h->Draw();
    fit_f->SetNpx(1000);
    fit_f->SetLineColor(kOrange);
    fit_f->SetLineWidth( 2 ); // 線の太さ変更
    fit_f->Draw("same");

    auto gr = new TGraph(5);
	Double_t x1 = par[1] - adc_n_sigma*par[2];
	Double_t x2 = peak_pos[peak_index]+500;
	Double_t y1 = 0;
	Double_t y2 = h->GetBinContent( h->GetMaximumBin() );
	gr->SetPoint(0,x1, y1);
	gr->SetPoint(1,x2, y1);
	gr->SetPoint(2,x2, y2);
	gr->SetPoint(3,x1, y2);
	gr->SetPoint(4,x1, y1);
	gr->Draw("F");
	gr->SetFillColor(kBlue);
    gr->SetFillStyle(3353);
	c->Update();

    h->GetXaxis()->SetNdivisions(505);
}

void fit_pedestal(TH1D *h, Double_t *par, TCanvas *c, Int_t n_c, Double_t n_threshold = 1000, Int_t range_extend = 2, Bool_t isFullRange = false) {

    c->cd(n_c);

    Int_t bin_info[4];
    searchRange(h, bin_info, 5, n_threshold);
    Int_t range_max = h->GetMaximumBin()+range_extend;
    if (isFullRange) range_max = bin_info[3];
    TF1 *pedestal_f = new TF1( "pedestal", "gausn", bin_info[1], range_max  );
    pedestal_f->SetParameter(1, h->GetMaximumBin());
    pedestal_f->SetParameter(2, 5);
    h->Fit(pedestal_f, "0", "", bin_info[1], range_max );
    for (int i = 0; i < 3; i++) par[i] = pedestal_f->GetParameter(i);

    gPad->SetLogy(1);
    h->GetXaxis()->SetRangeUser(bin_info[0]-10, bin_info[2]);
    h->Draw();
    pedestal_f->SetNpx(1000);
    pedestal_f->SetLineColor(kOrange);
    pedestal_f->SetLineWidth( 2 ); // 線の太さ変更
    pedestal_f->Draw("same");

}

Int_t pedestal_index(Int_t run_num) {
    Int_t index = 0;
    if (run_num < 44) {
        index = 0;
    } else if ( 44 < run_num && run_num < 57 ) {
        index = 1;
    } else if ( 57 < run_num && run_num < 79 ) {
        index = 2;
    } else if ( 79 < run_num ) {
        index = 3;
    }
    return index;
}

Double_t peak_search(TH1D *h, Int_t peak_n = 1, Int_t peak_index = 0) {

    TSpectrum *s = new TSpectrum(peak_n);
    s->Search(h, 5, "new");
    Double_t *peak_pos = s->GetPositionX();

    // // 三角マーカーを消す
    // TList *functions = h->GetListOfFunctions();
    // TPolyMarker *pm = ( TPolyMarker *)functions-> FindObject ( "TPolyMarker" );
    // functions->Remove(pm);
    // delete pm;

    return peak_pos[peak_index];
}