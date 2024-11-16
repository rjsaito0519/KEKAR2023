#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>

#include "TMath.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TColor.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include <TGClient.h>
#include <TGFrame.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TLine.h>
#include <TProfile.h>
#include <TBox.h>
#include <Math/RootFinder.h>
#include <TLatex.h>
#include <TF1Convolution.h>
#include <Math/RootFinder.h>

#include "variable.hh"
#include "param.hh"
#include "fit_functions.hh"

std::pair<Double_t, Double_t> cal_one_photon_gain(std::pair<Double_t, Double_t> mean, std::pair<Double_t, Double_t> pedestal, std::pair<Double_t, Double_t> n_pedestal, Double_t n_total) {
    Double_t mu = -TMath::Log( n_pedestal.first / n_total );
    Double_t gain = ( mean.first - pedestal.first ) / mu;
    
    // -- cal error propagation -----
    Double_t pdv_mean = 1/mu;
    Double_t pdv_ped = -1/mu;
    Double_t pdv_n_ped = gain / (n_pedestal.first * mu);
    Double_t error = TMath::Sqrt( TMath::Power(pdv_mean*mean.second, 2) + TMath::Power(pdv_ped*pedestal.second, 2) + TMath::Power(pdv_n_ped*n_pedestal.second, 2) );
    
    std::pair<Double_t, Double_t> one_photon_gain{gain, error};
    return one_photon_gain;
}

void write_opg_result(TString file_path, Int_t run_num, Int_t ch, const std::vector<Double_t>& result, Int_t header_mode) {
        // ファイルを削除
        if (std::ifstream(file_path.Data())) std::remove(file_path.Data());
        // 新しいファイルを作成し、追記モードで開く
        std::ofstream ofs(file_path.Data(), std::ios::app);
        if (header_mode == 0) {
            // sac one photon gain
            ofs << "run_num,ch,";
            for (Int_t i = 0, n = result.size(); i < n/2-1; i++) ofs << Form("p%d,", i);
            for (Int_t i = 0, n = result.size(); i < n/2-1; i++) ofs << Form("p%d_err,", i);
            ofs << "opg,opg_err\n";
        }

        ofs << run_num << "," << ch;
        for (Int_t i = 0, n = result.size(); i < n; i++) ofs << "," << result[i];
   
        ofs.close();
}
