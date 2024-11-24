#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    FitResult pedestal_fit_with_landau_conv(TH1D *h, TCanvas *c, Int_t n_c) {
        c->cd(n_c);
        std::vector<Double_t> par, err;
        
        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        Double_t fit_range_left  = peak_pos - 1.5*stdev;
        Double_t fit_range_right = peak_pos + 2.0*stdev;
        Double_t conv_range_left  = (fit_range_left - 2.0*stdev > 0) ? fit_range_left - 2.0*stdev : 0.0;
        Double_t conv_range_right = fit_range_right + 2.0*stdev;

        // -- scan and search best sigma of conv gaussian -----
        Double_t best_sigma = stdev;
        Double_t max_p_value = 0.0;

        Double_t sigma_left  = 0.01;
        Double_t sigma_right = stdev;
        Double_t n_scan = 11.0;
        Double_t step = (sigma_right-sigma_left)/(n_scan-1.0);
        for (Int_t dummy = 0; dummy < 3; dummy++) {
            for (Int_t i = 0; i < n_scan; i++) {
                Double_t sigma = sigma_left + i*step;
                TF1Convolution* f_conv = fit_functions::landau_gauss_conv(sigma, conv_range_left, conv_range_right);
                TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), f_conv, fit_range_left, fit_range_right, f_conv->GetNpar());
                f_fit->SetParameter(0, h->GetBinCenter(h->GetMaximumBin())*stdev );
                f_fit->SetParameter(1, peak_pos- 0.5*stdev);
                f_fit->SetParameter(2, stdev - sigma);
                h->Fit(f_fit, "0Q", "", fit_range_left, fit_range_right);
                Double_t chi2 = f_fit->GetChisquare();
                Double_t ndf  = f_fit->GetNDF();
                Double_t p_val = TMath::Prob(chi2, ndf);
                if (max_p_value < p_val) {
                    best_sigma = sigma;
                    max_p_value = p_val;
                }
                delete f_conv; delete f_fit;
            }
            sigma_left  = (best_sigma - step > 0.01) ? best_sigma - step : 0.01;
            sigma_right = (best_sigma + step < stdev) ? best_sigma + step : stdev;
            step = (sigma_right-sigma_left)/(n_scan-1.0);
        }        
        // std::cout << stdev << ", " << best_sigma <<  std::endl;
        
        TF1Convolution* f_conv = fit_functions::landau_gauss_conv(best_sigma, fit_range_left-stdev, fit_range_right+stdev);
        TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), f_conv, fit_range_left, fit_range_right, f_conv->GetNpar());
        f_fit->SetParameter(0, h->GetBinCenter(h->GetMaximumBin())/stdev );
        f_fit->SetParameter(1, peak_pos- 0.5*stdev);
        f_fit->SetParameter(2, stdev/2.0);
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(2);
        h->Fit(f_fit, "0Q", "", fit_range_left, fit_range_right);

        FitResult result;
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
        h->GetXaxis()->SetRangeUser(result.par[1] - 10.0*result.par[2], result.par[1] + 10.0*result.par[2]);
        h->Draw();
        f_fit->Draw("same");

        TLine *line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        line->SetLineStyle(2);
        line->SetLineColor(kRed);
        line->Draw("same");

        c->Update();

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult pedestal_fit_with_gumbel(TH1D *h, TCanvas *c, Int_t n_c) {
        c->cd(n_c);
        std::vector<Double_t> par, err;

        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        Double_t fit_range_left  = peak_pos - 1.5*stdev;
        Double_t fit_range_right = peak_pos + 2.0*stdev;

        // -- fit -----
        TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), fit_functions::gumbel, fit_range_left, fit_range_right, 3);
        f_fit->SetParameter(0, h->GetBinCenter(h->GetMaximumBin())/stdev );
        f_fit->SetParameter(1, peak_pos- 0.5*stdev);
        f_fit->SetParameter(2, stdev/2.0);
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(2);
        h->Fit(f_fit, "0Q", "", fit_range_left, fit_range_right);

        FitResult result;
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
        h->GetXaxis()->SetRangeUser(result.par[1] - 10.0*result.par[2], result.par[1] + 10.0*result.par[2]);
        h->Draw();
        f_fit->Draw("same");

        TLine *line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        line->SetLineStyle(2);
        line->SetLineColor(kRed);
        line->Draw("same");

        c->Update();

        return result;
    }


    // ____________________________________________________________________________________________
    FitResult pedestal_fit_with_gauss(TH1D *h, TCanvas *c, Int_t n_c, Double_t n_sigma) {
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

        // -- draw -----
        h->GetXaxis()->SetRangeUser(result.par[1] - 10.0*result.par[2], result.par[1] + 10.0*result.par[2]);
        h->Draw();
        f_fit->Draw("same");

        TLine *line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        line->SetLineStyle(2);
        line->SetLineColor(kRed);
        line->Draw("same");

        c->Update();

        return result;
    }

}