#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    FitResult cherenkov_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        // -- almost same as trigger_counter_tdc_fit -----
        Config& conf = Config::getInstance();
        c->cd(n_c);

        std::vector<Double_t> par, err;
        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        Double_t n_sigma = 5.0;

        // -- first fit -----
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-stdev, peak_pos+stdev);
        f_prefit->SetParameter(1, peak_pos);
        f_prefit->SetParameter(2, stdev/2);
        h->Fit(f_prefit, "0Q", "", peak_pos-stdev, peak_pos+stdev );
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
        delete f_prefit;

        // -- second fit -----
        TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), "gausn", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
        f_fit->SetParameter(0, par[0]);
        f_fit->SetParameter(1, par[1]);
        f_fit->SetParameter(2, par[2]*0.9);
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(2);
        h->Fit(f_fit, "0Q", "", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);

        FitResult result;
        par.clear();
        Int_t n_par = f_fit->GetNpar();
        for (Int_t i = 0; i < n_par; i++) {
            par.push_back(f_fit->GetParameter(i));
            err.push_back(f_fit->GetParError(i));
        }
        Double_t chi2 = f_fit->GetChisquare();
        Double_t ndf  = f_fit->GetNDF();
        result.par = par;
        result.err = err;
        result.reduced_chi2 = (Double_t) chi2/ndf;
        result.additional.push_back(result.par[1] - conf.cherenkov_tdc_n_sigma*result.par[2]);
        result.additional.push_back(result.par[1] + conf.cherenkov_tdc_n_sigma*result.par[2]);

        // -- draw -----
        h->GetXaxis()->SetRangeUser(
            result.par[1]-(conf.cherenkov_tdc_n_sigma+5.0)*result.par[2], 
            result.par[1]+(conf.cherenkov_tdc_n_sigma+5.0)*result.par[2]
        );
        h->Draw();
        f_fit->Draw("same");

        Double_t x1 = result.par[1] - conf.cherenkov_tdc_n_sigma  * result.par[2];
        Double_t x2 = result.par[1] + conf.cherenkov_tdc_n_sigma  * result.par[2];
        Double_t y1 = 0;
        Double_t y2 = h->GetBinContent(h->GetMaximumBin());

        TBox* box = new TBox(x1, y1, x2, y2);
        box->SetFillColor(kBlue);
        box->SetFillStyle(3353);
        box->Draw("same");

        c->Update();

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult correlation_fit(TH2D *h, TCanvas *c, Int_t n_c) {
        // -- almost same as trigger_counter_tdc_fit -----
        Config& conf = Config::getInstance();
        c->cd(n_c);
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.15);

        std::vector<Double_t> par, err;

        Double_t fit_range_min = 0.0;
        Double_t fit_range_max = 3500.0;

        // -- linear fit -----
        auto *f_linear = new TF1("linear", "[0]*x+[1]", fit_range_min, fit_range_max);
        f_linear->SetParameters(1., 0.);
        f_linear->SetLineColor(kRed);
        f_linear->SetLineWidth(2.0);
        h->Fit(f_linear, "0Q", "", fit_range_min, fit_range_max);

        FitResult result;
        par.clear();
        Int_t n_par = f_linear->GetNpar();
        for (Int_t i = 0; i < n_par; i++) {
            par.push_back(f_linear->GetParameter(i));
            err.push_back(f_linear->GetParError(i));
        }
        Double_t chi2 = f_linear->GetChisquare();
        Double_t ndf  = f_linear->GetNDF();
        result.par = par;
        result.err = err;
        result.reduced_chi2 = (Double_t) chi2/ndf;

        // -- draw -----
        h->GetXaxis()->SetRangeUser(0.0, 4096.0);
        h->Draw("colz");
        f_linear->Draw("same");
        c->Update();

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult poisson_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        c->cd(n_c);
        std::vector<Double_t> par, err;

        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        Double_t n_sigma = 1.5;

        // -- first fit -----
        Int_t n_iter = 3;
        par.insert(par.end(), {0.0, peak_pos, 2.0*stdev});
        for (Int_t dummy = 0; dummy < n_iter; dummy++) {
            TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
            f_prefit->SetParameter(1, par[1]);
            f_prefit->SetParameter(2, par[2]*0.5);
            h->Fit(f_prefit, "0Q", "", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
            par.clear();
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;            
        }

        Double_t fit_range_min = par[1]-n_sigma*par[2];
        Double_t fit_range_max = par[1]+n_sigma*par[2];
        TF1 *f_poisson = new TF1(Form("poisson_%s", h->GetName()), "[0]*TMath::Poisson(x, [1])", fit_range_min, fit_range_max);
        f_poisson->SetParameter(1, par[1]);
        TString fit_option = "0Q";
        if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
        h->Fit(f_poisson, fit_option.Data(), "", fit_range_min, fit_range_max);

        FitResult result;
        par.clear();
        Int_t n_par = f_poisson->GetNpar();
        for (Int_t i = 0; i < n_par; i++) {
            par.push_back(f_poisson->GetParameter(i));
            err.push_back(f_poisson->GetParError(i));
        }
        Double_t chi2 = f_poisson->GetChisquare();
        Double_t ndf  = f_poisson->GetNDF();
        result.par = par;
        result.err = err;
        result.reduced_chi2 = (Double_t) chi2/ndf;

        // -- draw -----
        h->GetXaxis()->SetRangeUser(-3.0, 3.0*result.par[1]);
        h->Draw();
        f_poisson->SetLineColor(kOrange);
        f_poisson->SetLineWidth(2);
        f_poisson->Draw("same");

        auto *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.1);
        text->DrawLatex(0.5, 0.8, Form("#lambda = %.2f", result.par[1]));
        text->Draw("same");

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult conv_poisson_fit(TH1D *h, TCanvas *c, Int_t n_c, Double_t pedestal_sigma) {
        c->cd(n_c);
        std::vector<Double_t> par, err;

        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        Double_t n_sigma = 1.5;

        // -- first fit -----
        Int_t n_iter = 3;
        par.insert(par.end(), {0.0, peak_pos, 2.0*stdev});
        for (Int_t dummy = 0; dummy < n_iter; dummy++) {
            TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
            f_prefit->SetParameter(1, par[1]);
            f_prefit->SetParameter(2, par[2]*0.5);
            h->Fit(f_prefit, "0Q", "", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
            par.clear();
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;            
        }

        Double_t fit_range_min = par[1]-n_sigma*par[2];
        Double_t fit_range_max = par[1]+n_sigma*par[2];
        Double_t conv_range_min = (fit_range_min - 2.0*stdev > 0) ? fit_range_min - 2.0*stdev : 0.0;
        Double_t conv_range_max = fit_range_max + 2.0*stdev;
        
        TF1Convolution* poisson_conv = fit_functions::poisson_conv(pedestal_sigma, conv_range_min, conv_range_max);
        TF1 *f_poisson = new TF1(Form("poisson_%s", h->GetName()), poisson_conv, fit_range_min, fit_range_max, 2);
        f_poisson->SetParameter(1, par[1]);
        TString fit_option = "0Q";
        if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
        h->Fit(f_poisson, fit_option.Data(), "", fit_range_min, fit_range_max);

        FitResult result;
        par.clear();
        Int_t n_par = f_poisson->GetNpar();
        for (Int_t i = 0; i < n_par; i++) {
            par.push_back(f_poisson->GetParameter(i));
            err.push_back(f_poisson->GetParError(i));
        }
        Double_t chi2 = f_poisson->GetChisquare();
        Double_t ndf  = f_poisson->GetNDF();
        result.par = par;
        result.err = err;
        result.reduced_chi2 = (Double_t) chi2/ndf;

        // -- draw -----
        h->GetXaxis()->SetRangeUser(-3.0, 3.0*result.par[1]);
        h->Draw();
        f_poisson->SetLineColor(kOrange);
        f_poisson->SetLineWidth(2);
        f_poisson->Draw("same");

        auto *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.1);
        text->DrawLatex(0.5, 0.8, Form("#lambda = %.2f", result.par[1]));
        text->Draw("same");

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult npe_gauss_fit(TH1D *h, TCanvas *c, Int_t n_c, Double_t n_sigma) {
        c->cd(n_c);
        std::vector<Double_t> par, err;

        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        
        // -- first fit -----
        Int_t n_iter = 3;
        par.insert(par.end(), {0.0, peak_pos, 2.0*stdev});
        for (Int_t dummy = 0; dummy < n_iter; dummy++) {
            TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
            f_prefit->SetParameter(1, par[1]);
            f_prefit->SetParameter(2, par[2]*0.5);
            h->Fit(f_prefit, "0Q", "", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
            par.clear();
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;            
        }

        Double_t fit_range_min = par[1]-n_sigma*par[2];
        Double_t fit_range_max = par[1]+n_sigma*par[2];
        TF1 *f_fit = new TF1(Form("gauss_%s", h->GetName()), "gausn", fit_range_min, fit_range_max);
        f_fit->SetParameters(par[0], par[1], par[2]);
        TString fit_option = "0Q";
        if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
        h->Fit(f_fit, fit_option.Data(), "", fit_range_min, fit_range_max);

        FitResult result;
        par.clear();
        Int_t n_par = f_fit->GetNpar();
        for (Int_t i = 0; i < n_par; i++) {
            par.push_back(f_fit->GetParameter(i));
            err.push_back(f_fit->GetParError(i));
        }
        Double_t chi2 = f_fit->GetChisquare();
        Double_t ndf  = f_fit->GetNDF();
        result.par = par;
        result.err = err;
        result.reduced_chi2 = (Double_t) chi2/ndf;

        // -- draw -----
        h->GetXaxis()->SetRangeUser(-3.0, 3.0*result.par[1]);
        h->Draw();
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(2);
        f_fit->Draw("same");

        auto *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.1);
        text->DrawLatex(0.5, 0.8, Form("#lambda = %.2f", result.par[1]));
        text->Draw("same");

        return result;
    }


}