#ifndef ANA_HELPER_
#define ANA_HELPER_

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <tuple>

// ROOTヘッダー
#include <Rtypes.h>
#include <TCanvas.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
#include <TLine.h>
#include <TBox.h>
#include <TROOT.h>
#include <TLatex.h>

#include "config.h"
#include "fit_functions.h"


// -- fitting result container -----
struct FitResult {
    std::vector<Double_t> par;
    std::vector<Double_t> err;
    Double_t reduced_chi2;
    Bool_t is_success;
    std::vector<Double_t> additional; // 何か追加で値を返したいときのコンテナ
};

// raw と trig をペアで管理する構造体
struct HistPair {
    TH1D* raw;  // raw histogram
    TH1D* trig; // trigged histogram

    HistPair(const TString& object_name, const TString& title, Int_t bins, Double_t range_min, Double_t range_max) {
        raw  = new TH1D(object_name+"_row",  title, bins, range_min, range_max);
        trig = new TH1D(object_name+"_trig", title, bins, range_min, range_max);
    }

    // // デストラクタでメモリを解放 (これをするとグラフが消えてしまう。ちゃんと動作はするので、グラフ化しないんだったらOK)
    // ~HistPair() {
    //     delete raw;
    //     delete trig;
    // }
};

namespace ana_helper {
    // -- general -----
    TCanvas* add_tab(TGTab *tab, const char* tabName);
    std::vector<Int_t> get_should_hit_ch(Int_t run_number);
    Int_t get_pedestal_run_num(Int_t run_number);

    // -- one photon gain -----
    std::pair<Double_t, Double_t> cal_one_photon_gain(
        std::pair<Double_t, Double_t> mean,
        std::pair<Double_t, Double_t> pedestal,
        std::pair<Double_t, Double_t> n_pedestal,
        Double_t n_total
    );
    
    // -- pedestal -----
    FitResult pedestal_fit_with_landau_conv(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult pedestal_fit_with_gumbel(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult pedestal_fit_with_gauss(TH1D *h, TCanvas *c, Int_t n_c, Double_t n_sigma = 1.0);

    // -- trigger counter -----
    FitResult trig_counter_adc_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult trig_counter_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c);

    // -- cherenkov counter -----
    FitResult cherenkov_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult correlation_fit(TH2D *h, TCanvas *c, Int_t n_c);
    FitResult poisson_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult conv_poisson_fit(TH1D *h, TCanvas *c, Int_t n_c, Double_t pedestal_sigma);
}

#endif  // ANA_HELPER_