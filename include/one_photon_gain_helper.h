#ifndef ONE_PHOTON_GAIN_HELPER_
#define ONE_PHOTON_GAIN_HELPER_

#include <utility>
#include <Rtypes.h>
#include <TMath.h>

std::pair<Double_t, Double_t> cal_one_photon_gain(std::pair<Double_t, Double_t> mean, 
                                                  std::pair<Double_t, Double_t> pedestal, 
                                                  std::pair<Double_t, Double_t> n_pedestal, 
                                                  Double_t n_total);

#endif  // ONE_PHOTON_GAIN_HELPER_
