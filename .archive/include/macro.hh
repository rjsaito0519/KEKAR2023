#ifndef MACRO_H_
#define MACRO_H_

void searchRange(TH1D *h, Int_t *bin_info, Int_t n = 5, Double_t n_threshold = 10);
void fit_TDC(TH1D *h, Double_t *par, TCanvas *c, Int_t n_c);
void fit_ADC(TH1D *h, Double_t *par, TCanvas *c, Int_t n_c, Int_t peak_index = 0);
void fit_pedestal(TH1D *h, Double_t *par, TCanvas *c, Int_t n_c, Double_t n_threshold = 1000, Int_t range_extend = 2, Bool_t isFullRange = false);
Int_t pedestal_index(Int_t run_num);
Double_t peak_search(TH1D *h, Int_t peak_n = 1, Int_t peak_index = 0);

#endif  // MACRO_H_