#ifndef ANA_HELPER_
#define ANA_HELPER_

TCanvas* add_tab(TGTab *tab, const char* tabName);
void write_result(TString file_path, Int_t run_num, const std::vector<std::vector<Double_t>>& result, Int_t header_mode);
Int_t pedestal_index(Int_t run_num);
Double_t peak_search(TH1D *h, Int_t n_peak = 1, Int_t nth_peak = 0);
std::vector<Double_t> search_range(TH1D *h, Int_t n_cluster = 5, Double_t ratio = 0.01);
std::vector<Double_t> fit_tdc(TH1D *h, TCanvas *c, Int_t n_c);
std::vector<Double_t> fit_trig_adc(TH1D *h, TCanvas *c, Int_t n_c, Int_t ch, TString key);
std::vector<Double_t> fit_pedestal(TH1D *h, TCanvas *c, Int_t n_c, Double_t n_sigma_left  = 2.0, Double_t n_sigma_right = 0.5);
std::vector<Double_t> fit_poisson(TH1D *h, TCanvas *c, Int_t n_c, Double_t fit_range_min = 2, Double_t fit_range_max = 10, Bool_t do_convolution = false);
std::vector<Double_t> fit_threshold(TH1D *h, TCanvas *c, Int_t n_c, Double_t fit_range_min = -5, Double_t fit_range_max = 15);
std::vector<Double_t> fit_adc(TH1D *h, TCanvas *c, Int_t n_c, Double_t peak_pos, Double_t half_width);

#endif  // ANA_HELPER_