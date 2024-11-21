#ifndef FIT_FUNCTIONS_H_
#define FIT_FUNCTIONS_H_

#include "TF1.h"
#include "TF1Convolution.h"
#include "TMath.h"

namespace fit_functions
{

    inline TF1Convolution* landau_gauss_conv(Double_t gauss_sigma, Double_t conv_range_min, Double_t conv_range_max)
    {
        auto *gauss  = new TF1("gauss", Form("TMath::Gaus(x, 0, %f, true)", gauss_sigma));
        auto *landau = new TF1("landau", "landaun");

        auto *landau_conv = new TF1Convolution(landau, gauss, conv_range_min, conv_range_max, true);
        landau_conv->SetRange(conv_range_min, conv_range_max);
        landau_conv->SetNofPointsFFT(5000);
        return landau_conv;
    }

    inline Double_t gumbel(Double_t *x, Double_t *par) {
        Double_t scale = par[0];
        Double_t x0    = par[1];
        Double_t sigma = par[2];

        Double_t z = (x[0] - x0) / sigma;
        return scale / sigma * TMath::Exp(-z) * TMath::Exp(-TMath::Exp(-z));
    }

}

#endif  // FIT_FUNCTIONS_H_
