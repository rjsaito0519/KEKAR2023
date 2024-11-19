#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    FitResult trig_counter_adc_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        c->cd(n_c);
        std::vector<Double_t> par, err;

        Double_t peak_pos = h->GetBinCenter( h->GetMaximumBin() );
        Double_t half_width = 50.0;

        // -- first fit -----
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-half_width, peak_pos+half_width);
        f_prefit->SetParameter(1, peak_pos);
        f_prefit->SetParameter(2, half_width/2);
        h->Fit(f_prefit, "0Q", "", peak_pos-half_width, peak_pos+half_width );
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
        delete f_prefit;

        // -- second fit -----
        Double_t fit_range_min = par[1] - 5.0*par[2];
        Double_t fit_range_max = par[1];
        TF1 *f_fit_erf = new TF1( Form("erf_fit_%s", h->GetName()), "[0]*TMath::Erf( (x-[1])/[2] ) + [3]", fit_range_min, fit_range_max);
        f_fit_erf->SetParameter(0, h->GetMaximum() / 2.0);
        f_fit_erf->SetParameter(1, par[1]-par[2]);
        f_fit_erf->SetParameter(2, 15.0);
        f_fit_erf->SetParameter(3, h->GetMaximum() / 2.0);
        f_fit_erf->SetLineColor(kOrange);
        f_fit_erf->SetLineWidth(2);
        h->Fit(f_fit_erf, "0Q", "", fit_range_min, fit_range_max);

        FitResult result;
        par.clear();
        Int_t n_par = f_fit_erf->GetNpar();
        for (Int_t i = 0; i < n_par; i++) {
            par.push_back(f_fit_erf->GetParameter(i));
            err.push_back(f_fit_erf->GetParError(i));
        }
        Double_t chi2 = f_fit_erf->GetChisquare();
        Double_t ndf  = f_fit_erf->GetNDF();
        result.par = par;
        result.err = err;
        result.reduced_chi2 = (Double_t) chi2/ndf;

        // 数値的に解を求めるためのRootFinderの設定, t0の値を調べる
        Config& conf = Config::getInstance();
        Double_t target_value = 2*par[0] * conf.trig_counter_adc_target_val_ratio;
        ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
        // 誤差関数の値は -par[0] ~ +par[0] の範囲になるので、par[0]分だけずらして最小値を0にする必要がある(したほうが解析しやすい)
        ROOT::Math::Functor1D erf_func([=](Double_t x) { return par[0]*TMath::Erf( (x-par[1])/par[2] ) + par[0] - target_value; });
        rootFinder.SetFunction(erf_func, fit_range_min, par[1]); // 探索範囲が第2, 3引数
        if (rootFinder.Solve()) {
            result.additional.push_back(rootFinder.Root());
            result.is_success = true;
        } else {
            std::cerr << "Error: Root finding failed.\n";
            result.additional.push_back(0.0);
            result.is_success = false;
        }

        
        h->GetXaxis()->SetRangeUser(result.additional[0] - 75.0, fit_range_max + 200.0);
        h->Draw();
        f_fit_erf->Draw("same");

        TLine *line = new TLine(result.additional[0], 0, result.additional[0], h->GetMaximum());
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kRed);
        line->Draw("same");

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult trig_counter_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        c->cd(n_c);

        std::vector<Double_t> par, err;
        Double_t peak_pos = h->GetBinCenter( h->GetMaximumBin() );
        Double_t half_width = 1000.0;
        Double_t n_sigma = 5.0;

        // -- first fit -----
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-half_width, peak_pos+half_width);
        f_prefit->SetParameter(1, peak_pos);
        f_prefit->SetParameter(2, half_width/2);
        h->Fit(f_prefit, "0Q", "", peak_pos-half_width, peak_pos+half_width );
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

        // -- draw -----
        Config& conf = Config::getInstance();
        h->GetXaxis()->SetRangeUser(
            result.par[1]-(conf.trig_counter_tdc_n_sigma +5.0)*result.par[2], 
            result.par[1]+(conf.trig_counter_tdc_n_sigma +5.0)*result.par[2]
        );
        h->Draw();
        f_fit->Draw("same");

        Double_t x1 = result.par[1] - conf.trig_counter_tdc_n_sigma  * result.par[2];
        Double_t x2 = result.par[1] + conf.trig_counter_tdc_n_sigma  * result.par[2];
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