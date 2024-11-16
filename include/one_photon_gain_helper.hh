#ifndef ONE_PHOTON_GAIN_HELPER_
#define ONE_PHOTON_GAIN_HELPER_

std::pair<Double_t, Double_t> cal_one_photon_gain(std::pair<Double_t, Double_t> mean, std::pair<Double_t, Double_t> pedestal, std::pair<Double_t, Double_t> n_pedestal, Double_t n_total);
void write_opg_result(TString file_path, Int_t run_num, Int_t ch, const std::vector<Double_t>& result, Int_t header_mode);

#endif  // ONE_PHOTON_GAIN_HELPER_