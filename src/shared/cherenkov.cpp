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

}