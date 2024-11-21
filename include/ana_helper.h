#ifndef ANA_HELPER_
#define ANA_HELPER_

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>

// ROOTヘッダー
#include <Rtypes.h>
#include <TCanvas.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
#include <TLine.h>
#include <TBox.h>

#include "config.h"
#include "fit_functions.h"


struct FitResult {
    std::vector<Double_t> par;
    std::vector<Double_t> err;
    Double_t reduced_chi2;
    Bool_t is_success;
    std::vector<Double_t> additional; // 何か追加で値を返したいときのコンテナ
};

namespace ana_helper {
    // -- general -----
    TCanvas* add_tab(TGTab *tab, const char* tabName);
    std::vector<Int_t> get_should_hit_seg(Int_t run_number);
    Double_t get_shower_adc_min(Int_t run_number, Int_t seg);

    // -- one photon gain -----
    std::pair<Double_t, Double_t> cal_one_photon_gain(
        std::pair<Double_t, Double_t> mean,
        std::pair<Double_t, Double_t> pedestal,
        std::pair<Double_t, Double_t> n_pedestal,
        Double_t n_total
    );
    
    std::vector<Double_t> fit_pedestal(TH1D *h, TCanvas *c, Int_t n_c);
    std::vector<Double_t> fit_pedestal2(TH1D *h, TCanvas *c, Int_t n_c) ;

    // -- trigger counter -----
    FitResult trig_counter_adc_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult trig_counter_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c);
}

#endif  // ANA_HELPER_