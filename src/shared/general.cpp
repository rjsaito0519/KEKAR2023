#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    TCanvas* add_tab(TGTab *tab, const char* tabName) {
        // タブを作成し、キャンバスを埋め込む
        TGCompositeFrame *tf = tab->AddTab(tabName);
        TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas(tabName, tf, 1000, 800);
        tf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        return embeddedCanvas->GetCanvas();
    }

    // ____________________________________________________________________________________________
    std::vector<Int_t> get_should_hit_seg(Int_t run_number) {

        static const std::unordered_map<Int_t, std::vector<Int_t>> kev_seg_map = {
            // Condition 1
            {300, {0}},       {301, {0}},       {302, {1}},       {303, {1, 2}},    {304, {2}},
            {305, {2, 3}},    {306, {3}},       {307, {}},        {308, {}},        {310, {0}},
            {312, {0}},       {313, {1}},       {314, {1, 2}},    {315, {2}},       {316, {2, 3}},
            {317, {3}},       {318, {}},        {319, {}},        {322, {0}},       {323, {0}},
            {324, {1}},       {325, {1, 2}},    {327, {2, 3}},    {328, {3}},       {329, {2}},
            {330, {}},        {331, {}},        {333, {0}},       {334, {0}},       {335, {1}},
            {336, {1, 2}},    {338, {2}},       {339, {2, 3}},    {340, {3}},       {341, {}},
            {342, {}},        {344, {0}},       {345, {0}},       {346, {1}},       {347, {1, 2}},
            {348, {2}},       {349, {2, 3}},    {350, {3}},       {351, {}},        {352, {3}},
            {353, {2, 3}},    {354, {}},        {357, {}},        {360, {1, 2}},    {361, {1}},
            {362, {}},        {363, {0}},       {365, {0}},       {366, {1}},       {368, {1, 2}},
            {371, {2}},       {372, {2, 3}},    {373, {3}},       {374, {}},        {376, {}},
            {378, {0}},       {379, {0}},       {380, {1}},       {381, {1, 2}},    {382, {2}},
            {383, {2, 3}},    {384, {3}},       {385, {}},        {386, {}},        {388, {2}},
            {390, {0}},       {391, {}},        {392, {0, 1}},

            // Condition 2
            {447, {1}},       {448, {1}},       {449, {1}},       {450, {1}},       {451, {1}},
            {452, {1}},       {453, {1}},       {454, {1}},       {455, {1}},       {457, {1, 2}},
            {458, {1, 2}},    {459, {1, 2}},    {460, {1, 2}},    {461, {1, 2}},    {462, {1, 2}},
            {463, {1, 2}},    {464, {1, 2}},    {465, {1, 2}},    {467, {2, 3}},    {468, {2, 3}},
            {469, {2, 3}},    {470, {2, 3}},    {471, {2, 3}},    {472, {2, 3}},    {473, {2, 3}},
            {474, {2, 3}},    {475, {2, 3}},    {477, {3}},       {478, {3}},       {479, {3}},
            {480, {3}},       {481, {3}},       {482, {3}},       {483, {3}},       {484, {3}},
            {485, {3}},       {487, {0, 1}},    {488, {0, 1}},    {489, {0, 1}},    {490, {0, 1}},
            {491, {0, 1}},    {492, {0, 1}},    {493, {0, 1}},    {494, {0, 1}},    {495, {0, 1}},
            {497, {0}},       {498, {0}},       {499, {0}},       {500, {0}},       {501, {0}},
            {502, {0}},       {503, {0}},       {504, {0}},       {505, {0}},       {507, {}},
            {508, {}},        {509, {}},        {510, {}},        {511, {}},        {512, {}},
            {513, {}},        {514, {}},        {515, {}},        {517, {}},        {518, {}},
            {519, {}},        {521, {0}},       {522, {0}},       {523, {0}},       {525, {2, 3}},
            {526, {2, 3}},    {527, {2, 3}},    {529, {3}},       {530, {3}},       {531, {3}}
        };

        auto it = kev_seg_map.find(run_number);
        return (it != kev_seg_map.end()) ? it->second : std::vector<Int_t>{};
    }

    // ____________________________________________________________________________________________
    Double_t get_shower_adc_min(Int_t run_number, Int_t seg) {

        static const std::vector<Double_t> adc_min_condition1{  300.0, 300.0, 300.0, 300.0 };
        static const std::vector<Double_t> adc_min_condition2{ 1000.0, 600.0, 700.0, 600.0 };

        if (run_number <= 392) {
            return adc_min_condition1[seg];
        } else {
            return adc_min_condition2[seg];
        }

    }

    // ____________________________________________________________________________________________
    std::pair<Double_t, Double_t> cal_one_photon_gain(std::pair<Double_t, Double_t> mean, std::pair<Double_t, Double_t> pedestal, std::pair<Double_t, Double_t> n_pedestal, Double_t n_total) {
        Double_t mu = -TMath::Log( n_pedestal.first / n_total );
        Double_t gain = ( mean.first - pedestal.first ) / mu;
        
        // -- cal error propagation -----
        Double_t pdv_mean = 1.0/mu;
        Double_t pdv_ped = -1.0/mu;
        Double_t pdv_n_ped = gain / (n_pedestal.first * mu);
        Double_t error = TMath::Sqrt( TMath::Power(pdv_mean*mean.second, 2.0) + TMath::Power(pdv_ped*pedestal.second, 2.0) + TMath::Power(pdv_n_ped*n_pedestal.second, 2.0) );
        
        std::pair<Double_t, Double_t> one_photon_gain{gain, error};
        return one_photon_gain;
    }


    // ____________________________________________________________________________________________
    std::vector<Double_t> fit_pedestal(TH1D *h, TCanvas *c, Int_t n_c) {
        c->cd(n_c);

        std::vector<Double_t> result{};
        
        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        Double_t fit_range_left  = peak_pos - 1.5*stdev;
        Double_t fit_range_right = peak_pos + 2.0*stdev;

        Double_t conv_range_left  = (fit_range_left - 2.0*stdev > 0) ? fit_range_left - 2.0*stdev : 0.0;
        Double_t conv_range_right = fit_range_right + 2.0*stdev;

        // -- second fit -----
        Double_t best_sigma = stdev;
        Double_t max_p_value = 0.0;

        Double_t sigma_left  = 0.01;
        Double_t sigma_right = stdev;
        Double_t n_scan = 10.0;
        Double_t step = (sigma_right-sigma_left)/(n_scan-1.0);
        for (Int_t dummy = 0; dummy < 3; dummy++) {
            // for (Double_t sigma = sigma_left; sigma < sigma_right; sigma += step) {
            for (Int_t i = 0; i < n_scan; i++) {
                Double_t sigma = sigma_left + i*step;
                TF1Convolution* f_conv = fit_functions::landau_gauss_conv(sigma, conv_range_left, conv_range_right);
                TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), f_conv, fit_range_left, fit_range_right, f_conv->GetNpar());
                f_fit->SetParameter(0, h->GetBinCenter(h->GetMaximumBin())*stdev );
                f_fit->SetParameter(1, peak_pos- 0.5*stdev);
                f_fit->SetParameter(2, stdev - sigma);
                h->Fit(f_fit, "0", "", fit_range_left, fit_range_right);
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
        
        std::cout << stdev << ", " << best_sigma <<  std::endl;

        
        TF1Convolution* f_conv = fit_functions::landau_gauss_conv(best_sigma, fit_range_left-stdev, fit_range_right+stdev);
        TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), f_conv, fit_range_left, fit_range_right, f_conv->GetNpar());
        f_fit->SetParameter(0, h->GetBinCenter(h->GetMaximumBin())/stdev );
        f_fit->SetParameter(1, peak_pos- 0.5*stdev);
        f_fit->SetParameter(2, stdev/2.0);
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(2);
        h->Fit(f_fit, "0", "", fit_range_left, fit_range_right);

        result.clear();
        for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParameter(i));
        for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParError(i));
        Double_t chi2 = f_fit->GetChisquare();
        Double_t ndf  = f_fit->GetNDF();
        Double_t p_value = TMath::Prob(chi2, ndf);
        // result.push_back(p_value);

        Double_t n_sigma_left = 10.0;
        Double_t n_sigma_right = 15.0;
        
        // -- draw -----
        h->GetXaxis()->SetRangeUser(result[1] - (n_sigma_left+5)*result[2], result[1]+(n_sigma_right+5)*result[2]);
        h->Draw();
        f_fit->Draw("same");

        TLine *line = new TLine(result[1], 0, result[1], h->GetMaximum());
        line->SetLineStyle(2);
        line->SetLineColor(kRed);
        line->Draw("same");

        c->Update();

        return result;

    }


    // ____________________________________________________________________________________________
    std::vector<Double_t> fit_pedestal2(TH1D *h, TCanvas *c, Int_t n_c) {
        c->cd(n_c);

        std::vector<Double_t> result{};
        
        Double_t peak_pos = h->GetMean();
        Double_t stdev = h->GetStdDev();
        Double_t fit_range_left  = peak_pos - 1.5*stdev;
        Double_t fit_range_right = peak_pos + 2.0*stdev;

        // -- second fit -----
        TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), fit_functions::gumbel, fit_range_left, fit_range_right, 3);
        f_fit->SetParameter(0, h->GetBinCenter(h->GetMaximumBin())/stdev );
        f_fit->SetParameter(1, peak_pos- 0.5*stdev);
        f_fit->SetParameter(2, stdev/2.0);
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(2);
        h->Fit(f_fit, "0", "", fit_range_left, fit_range_right);

        result.clear();
        for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParameter(i));
        for (Int_t i = 0; i < 3; i++) result.push_back(f_fit->GetParError(i));
        Double_t chi2 = f_fit->GetChisquare();
        Double_t ndf  = f_fit->GetNDF();
        Double_t p_value = TMath::Prob(chi2, ndf);
        // result.push_back(p_value);

        Double_t n_sigma_left = 10.0;
        Double_t n_sigma_right = 15.0;
        
        // -- draw -----
        h->GetXaxis()->SetRangeUser(result[1] - (n_sigma_left+5)*result[2], result[1]+(n_sigma_right+5)*result[2]);
        h->Draw();
        f_fit->Draw("same");

        TLine *line = new TLine(result[1], 0, result[1], h->GetMaximum());
        line->SetLineStyle(2);
        line->SetLineColor(kRed);
        line->Draw("same");

        c->Update();

        return result;

    }


}