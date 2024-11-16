#ifndef FIT_FUNCTIONS_H_
#define FIT_FUNCTIONS_H_

namespace fit_functions
{

    inline TF1Convolution f_poisson_conv(Double_t sigma, Double_t range_min, Double_t range_max, Int_t NofPointsFFT = 10000)
    {
        auto *gauss   = new TF1("gauss", Form("TMath::Gaus(x, 0, %f, true)", sigma));
        auto *poisson = new TF1("poisson", "[0]*TMath::Poisson(x, [1])", range_min, range_max);
        auto *poisson_conv = new TF1Convolution(poisson, gauss, range_min-100, range_max+100, true);
        poisson_conv->SetNofPointsFFT(NofPointsFFT);
        return *poisson_conv;
    }

}


#endif  // FIT_FUNCTIONS_H_