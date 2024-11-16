#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>

#include "TMath.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TColor.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include <TGClient.h>
#include <TGFrame.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TLine.h>
#include <TProfile.h>
#include <TBox.h>
#include <Math/RootFinder.h>
#include <TLatex.h>
#include <TF1Convolution.h>
#include <Math/RootFinder.h>

#include "variable.hh"
#include "param.hh"
#include "fit_functions.hh"

// タブを追加し、埋め込まれたキャンバスを返す関数
TCanvas* add_tab(TGTab *tab, const char* tabName) {
    // タブを作成し、キャンバスを埋め込む
    TGCompositeFrame *tf = tab->AddTab(tabName);
    TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas(tabName, tf, 1000, 800);
    tf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    return embeddedCanvas->GetCanvas();
}

void write_result(TString file_path, Int_t run_num, const std::vector<std::vector<Double_t>>& result, Int_t header_mode) {
        // ファイルを削除
        if (std::ifstream(file_path.Data())) std::remove(file_path.Data());
        // 新しいファイルを作成し、追記モードで開く
        std::ofstream ofs(file_path.Data(), std::ios::app);
        if (header_mode == 0) {
            // gaussian
            ofs << "run_num,ch,amp,center,sigma,amp_err,center_err,sigma_err\n";
        } else if (header_mode == 1) {
            // poisson
            ofs << "run_num,ch,amp,mu,amp_err,mu_err,conv_sigma\n";
        } else if (header_mode == 2) {
            // threshold
            ofs << "run_num,ch,p0,p1,p0_err,p1_err\n";
        } else if (header_mode == 3) {
            // sac one photon gain
            ofs << "run_num,ch,";
            for (Int_t i = 0, n = result[0].size(); i < n/2-1; i++) ofs << Form("p%d,", i);
            for (Int_t i = 0, n = result[0].size(); i < n/2-1; i++) ofs << Form("p%d_err,", i);
            ofs << "opg,opg_err\n";
        }
        for (Int_t ch = 0, len_data = result.size(); ch < len_data; ch++) {
            ofs << run_num << "," << ch;
            for (Int_t i = 0, n = result[ch].size(); i < n; i++) ofs << "," << result[ch][i];
            ofs << "\n";
        }
        ofs.close();
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

Double_t peak_search(TH1D *h, Int_t n_peak = 1, Int_t nth_peak = 0) {

    TSpectrum *s = new TSpectrum(n_peak);
    s->Search(h, 5, "new");
    Double_t *peak_pos_container = s->GetPositionX();
    Double_t peak_pos = peak_pos_container[nth_peak];
    delete s;

    // 三角マーカーを消す
    TList *functions = h->GetListOfFunctions();
    TPolyMarker *pm = ( TPolyMarker *)functions-> FindObject ( "TPolyMarker" );
    functions->Remove(pm);
    delete pm;

    return peak_pos;
}

std::vector<Double_t> search_range(TH1D *h, Int_t n_cluster = 5, Double_t ratio = 0.01) {
    Int_t max_bin_pos = h->GetMaximumBin(), min_bin_pos = h->GetMaximumBin();
    Double_t tmp_value = 10000000;
    Double_t buf[n_cluster];
    Double_t threshold = h->GetBinContent( h->GetMaximumBin() )*ratio;

    // left side
    while (tmp_value > threshold) {
        min_bin_pos -= n_cluster;
        for (Int_t i = 0; i < n_cluster; i++) buf[i] = h->GetBinContent(min_bin_pos-i);
        tmp_value = *std::max_element(buf, buf + n_cluster);
    }
    // right side
    tmp_value = 10000000;
    while (tmp_value > threshold) {
        max_bin_pos += n_cluster;
        for (Int_t i = 0; i < n_cluster; i++) buf[i] = h->GetBinContent(max_bin_pos+i);
        tmp_value = *std::max_element(buf, buf + n_cluster);
    }

    std::vector<Double_t> range_container{ h->GetBinCenter(min_bin_pos), h->GetBinCenter(max_bin_pos) };
    return range_container;
}

std::vector<Double_t> fit_tdc(TH1D *h, TCanvas *c, Int_t n_c) {
    c->cd(n_c);

    std::vector<Double_t> result{};
    Double_t peak_pos = h->GetBinCenter( h->GetMaximumBin() );
    Double_t half_width = 10.0;
    Double_t n_sigma = 2;

    // -- first fit -----
    TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-half_width, peak_pos+half_width);
    f_prefit->SetParameter(1, peak_pos);
    f_prefit->SetParameter(2, half_width/2);
    TString fit_option = "0Q";
    if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
    h->Fit(f_prefit, fit_option.Data(), "", peak_pos-half_width, peak_pos+half_width );
    for (Int_t i = 0; i < 3; i++) result.push_back(f_prefit->GetParameter(i));
    delete f_prefit;

    // -- second fit -----
    TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), "gausn", result[1]-n_sigma*result[2], result[1]+n_sigma*result[2]);
    f_fit->SetParameter(1, result[1]);
    f_fit->SetParameter(2, result[2]*0.9);
    f_fit->SetLineColor(kOrange);
    f_fit->SetLineWidth(2);
    h->Fit(f_fit, fit_option.Data(), "", result[1]-n_sigma*result[2], result[1]+n_sigma*result[2]);
    result.clear();
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParameter(i));
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParError(i));
    Double_t chi2 = f_fit->GetChisquare();
    Double_t ndf  = f_fit->GetNDF();
    Double_t p_value = TMath::Prob(chi2, ndf);
    // result.push_back(p_value);

    // -- draw -----
    h->GetXaxis()->SetRangeUser(result[1]-(tdc_gate_n_sigma+5)*result[2], result[1]+(tdc_gate_n_sigma+3)*result[2]);
    h->Draw();
    f_fit->Draw("same");

    Double_t x1 = result[1] - tdc_gate_n_sigma * result[2];
    Double_t x2 = result[1] + tdc_gate_n_sigma * result[2];
    Double_t y1 = 0;
    Double_t y2 = h->GetBinContent(h->GetMaximumBin());

    TBox* box = new TBox(x1, y1, x2, y2);
    box->SetFillColor(kBlue);
    box->SetFillStyle(3353);
    box->Draw("same");

    c->Update();

    return result;
}

std::vector<Double_t> fit_trig_adc(TH1D *h, TCanvas *c, Int_t n_c, Int_t ch, TString key) {
    c->cd(n_c);

    std::vector<Double_t> result{};
    std::vector<Double_t> fit_par = param::trigger_adc.at("default");
    auto it = param::trigger_adc.find(key.Data());
    if (it != param::trigger_adc.end()) fit_par = it->second;
    Double_t peak_pos  = fit_par[3*ch];
    Double_t range_min = fit_par[3*ch+1];
    Double_t range_max = fit_par[3*ch+2];
    Double_t n_sigma_left  = 1.5;
    Double_t n_sigma_right = 2.0;

    // -- first fit -----
    TF1 *f_prefit = new TF1("pre_fit_gauss", "landaun", range_min, range_max);
    f_prefit->SetParameter(1, peak_pos);
    f_prefit->SetParameter(2, (range_max-range_min)/2);
    TString fit_option = "0Q";
    if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
    h->Fit(f_prefit, fit_option.Data(), "", range_min, range_max );
    for (Int_t i = 0; i < 3; i++) result.push_back(f_prefit->GetParameter(i));
    delete f_prefit;

    // -- second fit -----
    TF1 *f_fit = new TF1( Form("landau_%s", h->GetName()), "landaun", result[1]-n_sigma_left*result[2], result[1]+n_sigma_right*result[2]);
    f_fit->SetParameter(1, result[1]);
    f_fit->SetParameter(2, result[2]*0.9);
    f_fit->SetLineColor(kOrange);
    f_fit->SetLineWidth(2);
    h->Fit(f_fit, fit_option.Data(), "", result[1]-n_sigma_left*result[2], result[1]+n_sigma_right*result[2]);
    result.clear();
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParameter(i));
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParError(i));
    Double_t chi2 = f_fit->GetChisquare();
    Double_t ndf  = f_fit->GetNDF();
    Double_t p_value = TMath::Prob(chi2, ndf);
    // result.push_back(p_value);

    // -- draw -----
    h->GetXaxis()->SetRangeUser(result[1]-(n_sigma_left+5)*result[2], result[1]+(n_sigma_right+10)*result[2]);
    h->Draw();
    f_fit->Draw("same");

    Double_t x1 = result[1] - adc_gate_n_sigma * result[2];
    Double_t x2 = result[1] + (n_sigma_right+10) * result[2];
    Double_t y1 = 0;
    Double_t y2 = f_fit->GetMaximum();

    TBox* box = new TBox(x1, y1, x2, y2);
    box->SetFillColor(kBlue);
    box->SetFillStyle(3353);
    box->Draw("same");

    c->Update();

    return result;
}

std::vector<Double_t> fit_pedestal(TH1D *h, TCanvas *c, Int_t n_c, Double_t n_sigma_left  = 2.0, Double_t n_sigma_right = 0.5) {
    c->cd(n_c);

    std::vector<Double_t> result{};
    Double_t peak_pos = h->GetBinCenter( h->GetMaximumBin() );
    Double_t half_width = 10.0;

    // -- first fit -----
    TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-half_width, peak_pos+half_width);
    f_prefit->SetParameter(1, peak_pos);
    f_prefit->SetParameter(2, half_width/2);
    TString fit_option = "0Q";
    if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
    h->Fit(f_prefit, fit_option.Data(), "", peak_pos-half_width, peak_pos+half_width );
    for (Int_t i = 0; i < 3; i++) result.push_back(f_prefit->GetParameter(i));
    delete f_prefit;

    // -- second fit -----
    TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), "gausn", result[1]-n_sigma_left*result[2], result[1]+n_sigma_right*result[2]);
    f_fit->SetParameter(1, result[1]);
    f_fit->SetParameter(2, result[2]*0.9);
    f_fit->SetLineColor(kOrange);
    f_fit->SetLineWidth(2);
    h->Fit(f_fit, fit_option.Data(), "", result[1]-n_sigma_left*result[2], result[1]+n_sigma_right*result[2]);
    result.clear();
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParameter(i));
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParError(i));
    Double_t chi2 = f_fit->GetChisquare();
    Double_t ndf  = f_fit->GetNDF();
    Double_t p_value = TMath::Prob(chi2, ndf);
    // result.push_back(p_value);

    // -- draw -----
    h->GetXaxis()->SetRangeUser(result[1]-(n_sigma_left+5)*result[2], result[1]+(n_sigma_right+5)*result[2]);
    h->Draw();
    f_fit->Draw("same");

    TLine *line = new TLine(result[1], 0, result[1], h->GetMaximum());
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw("same");

    c->Update();

    return result;

}

std::vector<Double_t> fit_poisson(TH1D *h, TCanvas *c, Int_t n_c, Double_t fit_range_min = 2, Double_t fit_range_max = 10, Bool_t do_convolution = false) {
    c->cd(n_c);
    // gPad->SetLogy(1);
    std::vector<Double_t> result{};
    Double_t conv_sigma = 0.4;

    TF1 *f_poisson;
    if (do_convolution) {
        auto poisson_conv = fit_functions::f_poisson_conv(conv_sigma, fit_range_min, fit_range_max);
        f_poisson = new TF1(Form("poisson_%s", h->GetName()), poisson_conv, fit_range_min, fit_range_max, 2);
    } else {
        f_poisson = new TF1(Form("poisson_%s", h->GetName()), "[0] * TMath::Poisson(x, [1])", fit_range_min, fit_range_max);
    }

    f_poisson->SetParameter(1, h->GetMean());
    TString fit_option = "0Q";
    if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
    h->Fit(f_poisson, fit_option.Data(), "", fit_range_min, fit_range_max);
    for (Int_t i = 0; i < 2; i++) result.push_back(f_poisson->GetParameter(i));
    for (Int_t i = 0; i < 2; i++) result.push_back(f_poisson->GetParError(i));
    Double_t chi2 = f_poisson->GetChisquare();
    Double_t ndf  = f_poisson->GetNDF();
    Double_t p_value = TMath::Prob(chi2, ndf);
    // result.push_back(p_value);
    result.push_back( do_convolution ? conv_sigma : 0.0 );

    std::vector<Double_t> draw_range = search_range(h);
    h->GetXaxis()->SetRangeUser(-3, draw_range[1]+3);
    h->Draw();

    f_poisson->SetLineColor(kOrange);
    f_poisson->Draw("same");

    auto *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.1);
    text->DrawLatex(0.35, 0.8, Form("#mu = %.2f", result[1]));
    text->Draw("same");

    return result;
}

std::vector<Double_t> fit_threshold(TH1D *h, TCanvas *c, Int_t n_c, Double_t fit_range_min = -5, Double_t fit_range_max = 15) {
    c->cd(n_c);
    std::vector<Double_t> result{};

    auto *f_threshold = new TF1(Form("threshold_%s", h->GetName()), "1 / ( 1+TMath::Exp(-[0]*(x-[1])) )", fit_range_min, fit_range_max);
    f_threshold->SetParameters(0.1, 1.0);
    h->Fit(f_threshold, "0Q", "", fit_range_min, fit_range_max);
    for (Int_t i = 0; i < 2; i++) result.push_back(f_threshold->GetParameter(i));
    for (Int_t i = 0; i < 2; i++) result.push_back(f_threshold->GetParError(i));
    Double_t chi2 = f_threshold->GetChisquare();
    Double_t ndf  = f_threshold->GetNDF();
    Double_t p_value = TMath::Prob(chi2, ndf);
    // result.push_back(p_value);

    // 数値的に解を求めるためのRootFinderの設定, ほぼ1になる(0.999)の値を調べる
    Double_t target_value = 0.999;
    ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
    ROOT::Math::Functor1D step_func([=](Double_t x) { return 1 / ( 1+TMath::Exp(-result[0]*(x-result[1])) ) - target_value; });
    rootFinder.SetFunction(step_func, result[1], fit_range_max); // 探索範囲が第2, 3引数
    Bool_t is_success = rootFinder.Solve();
    Double_t higher_thre = is_success ? rootFinder.Root() : result[1];

    h->GetXaxis()->SetRangeUser(fit_range_min, fit_range_max);
    h->Draw();

    f_threshold->SetLineColor(kOrange);
    f_threshold->Draw("same");

    auto *text_mid = new TLatex();
    text_mid->SetNDC();
    text_mid->SetTextSize(0.05);
    text_mid->DrawLatex(0.35, 0.5, Form("NPE_{ mid th} = %.2f", result[1]));
    text_mid->Draw("same");

    auto *text_high = new TLatex();
    text_high->SetNDC();
    text_high->SetTextSize(0.05);
    text_high->DrawLatex(0.35, 0.45, Form("NPE_{high th} = %.2f", higher_thre));
    text_high->Draw("same");


    auto *line_mid = new TLine(result[1], 0, result[1], 1.);
    line_mid->SetLineStyle(2);
    line_mid->SetLineWidth(2);
    line_mid->SetLineColor(kRed);
    line_mid->Draw("same");

    auto *line_high = new TLine(higher_thre, 0, higher_thre, 1.);
    line_high->SetLineStyle(2);
    line_high->SetLineWidth(2);
    line_high->SetLineColor(kRed);
    line_high->Draw("same");

    return result;
}

std::vector<Double_t> fit_adc(TH1D *h, TCanvas *c, Int_t n_c, Double_t peak_pos, Double_t half_width) {
    c->cd(n_c);

    std::vector<Double_t> result{};
    Double_t n_sigma_left  = 2.0;
    Double_t n_sigma_right = 1.5;

    // -- first fit -----
    TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-half_width, peak_pos+half_width);
    f_prefit->SetParameter(1, peak_pos);
    f_prefit->SetParameter(2, half_width/2);
    TString fit_option = "0Q";
    if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
    h->Fit(f_prefit, fit_option.Data(), "", peak_pos-half_width, peak_pos+half_width );
    for (Int_t i = 0; i < 3; i++) result.push_back(f_prefit->GetParameter(i));
    delete f_prefit;

    // -- second fit -----
    TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), "gausn", result[1]-n_sigma_left*result[2], result[1]+n_sigma_right*result[2]);
    f_fit->SetParameter(1, result[1]);
    f_fit->SetParameter(2, result[2]*0.9);
    f_fit->SetLineColor(kOrange);
    f_fit->SetLineWidth(2);
    h->Fit(f_fit, fit_option.Data(), "", result[1]-n_sigma_left*result[2], result[1]+n_sigma_right*result[2]);
    result.clear();
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParameter(i));
    for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParError(i));
    Double_t chi2 = f_fit->GetChisquare();
    Double_t ndf  = f_fit->GetNDF();
    Double_t p_value = TMath::Prob(chi2, ndf);
    // result.push_back(p_value);

    // -- draw -----
    h->GetXaxis()->SetRangeUser(result[1]-(n_sigma_left+5)*result[2], result[1]+(n_sigma_right+5)*result[2]);
    h->Draw();
    f_fit->Draw("same");

    TLine *line = new TLine(result[1], 0, result[1], h->GetMaximum());
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw("same");

    c->Update();

    return result;

}

