#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
#include <sstream>
#include <string>
#include <unistd.h> // for getcwd

// ROOT libraries
#include <TFile.h>
#include <TTree.h>
#include <TEventList.h>
#include <TMath.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TColor.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TSpectrum.h>
#include <TPolyMarker.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TParticle.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TMarker.h>
#include <TSpline.h>
#include <RooBreitWigner.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMinuit.h>
#include <TF1Convolution.h>
#include <TComplex.h>
#include <TProfile.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TLine.h>
#include <Math/IntegratorOptions.h>

#include "config.h"
#include "variable.hh"
#include "ana_helper.hh"
#include "param.hh"
#include "one_photon_gain_helper.hh"

static const TString pdf_name = "../img/bac_one_photon_gain.pdf";

std::pair<Double_t, Double_t> analyze(Int_t run_num, Int_t ch, TVirtualPad *c, Int_t n_c)
{   
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x"); // x軸のタイトルサイズ
    gStyle->SetTitleSize(0.06, "y"); // y軸のタイトルサイズ
    gStyle->SetTitleFontSize(0.06);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // -- parameter -----
    Config& conf = Config::getInstance();

    // +----------------+
    // | load root file |
    // +----------------+
    TString root_file_path = Form( "../root/kekar_run%05d.root", run_num );
    auto *f = new TFile( root_file_path.Data() );
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path << std::endl;
        return std::make_pair(-1.0, -1.0);
    }
    TTreeReader reader("tree", f);
    TTreeReaderArray<Double_t> baca(reader, "baca");

    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    auto *h = new TH1D(Form("BACa_%d_%d", run_num, ch+1), Form("run%05d BAC(ADC) ch%d;ADC;", run_num, ch+1), conf.adc_bin_num, conf.adc_min, conf.adc_max);

    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    reader.Restart();
    while (reader.Next()){
        h->Fill(baca[ch]);
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare -----
    c->cd(n_c);
    
    std::vector<Double_t> result{};
    std::pair<Double_t, Double_t> mean{ h->GetMean(),  h->GetMeanError() };
    Double_t n_total = h->GetEntries();

    // -- load and set parameter -----
    TString key = Form("%05d-%d", run_num, ch);
    std::vector<Double_t> par = param::bac_opg.count(key.Data()) ? param::bac_opg.at(key.Data()) : param::bac_opg.at("default");
    Double_t n_gauss         = par[0];
    Double_t first_peak_pos  = par[1];
    Double_t fit_range_left  = par[2];
    Double_t fit_range_right = par[3];
    // std::cout << n_gauss << ", " << first_peak_pos << ", " << fit_range_left << ", " << fit_range_right << std::endl;
    h->GetXaxis()->SetRangeUser(fit_range_left-5.0, fit_range_right+5.0);
    Double_t half_width = 5.0;

    // -- first fit -----
    TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", first_peak_pos-half_width, first_peak_pos+half_width);
    f_prefit->SetParameter(1, first_peak_pos);
    f_prefit->SetParameter(2, half_width);
    h->Fit(f_prefit, "0Q", "", first_peak_pos-half_width, first_peak_pos+half_width );
    for (Int_t i = 0; i < 3; i++) result.push_back(f_prefit->GetParameter(i));
    delete f_prefit;

    // -- second fit -----
    TString func_str = "[0]*TMath::Gaus(x, [1], [2], true)";
    for (Int_t i = 1; i < n_gauss; i++) func_str += Form(" + [%d]*TMath::Gaus(x, [1]+%d*[3], TMath::Sqrt(%d)*[4], true)", i+4, i, i+1);
    TF1 *f_fit = new TF1( Form("gauss_%s", h->GetName()), func_str.Data(), fit_range_left, fit_range_right);
    f_fit->SetParameter(0, result[0]);
    f_fit->SetParameter(1, result[1]);
    f_fit->SetParameter(2, result[2]*0.9);
    Double_t peak_to_peak = (fit_range_right-fit_range_left-2*result[2])/n_gauss;
    f_fit->SetParameter(3, peak_to_peak);
    f_fit->SetParameter(4, result[2]);
    for (Int_t i = 1; i < n_gauss; i++) {
        h->GetBinCenter( h->GetMaximumBin() );
        Double_t height = h->GetBinContent(h->FindBin( result[1]+i*peak_to_peak ));
        f_fit->SetParameter(i+4, height*TMath::Sqrt(i+1)*result[2]);
        f_fit->SetParLimits(i+4, 0, 0x10000000);
    }
    f_fit->SetLineColor(kOrange);
    f_fit->SetLineWidth(2);
    h->Fit(f_fit, "0", "", fit_range_left, fit_range_right);
    result.clear();
    Int_t n_par = f_fit->GetNpar();
    for (Int_t i = 0; i < n_par; i++) result.push_back(i < f_fit->GetNpar() ? f_fit->GetParameter(i) : TMath::QuietNaN());
    for (Int_t i = 0; i < n_par; i++) result.push_back(i < f_fit->GetNpar() ? f_fit->GetParError(i)  : TMath::QuietNaN());

    std::cout << result[3]*n_gauss + 2*result[2] + fit_range_left << std::endl;

    // -- cal one photon gain -----
    std::pair<Double_t, Double_t> pedestal{ f_fit->GetParameter(1), f_fit->GetParError(1) };
    std::pair<Double_t, Double_t> n_pedestal{ f_fit->GetParameter(0), f_fit->GetParError(0) };
    std::pair<Double_t, Double_t> one_photon_gain = cal_one_photon_gain( mean, pedestal, n_pedestal, n_total);
    std::cout << one_photon_gain.first << ", " << one_photon_gain.second << std::endl;
    result.push_back(one_photon_gain.first);
    result.push_back(one_photon_gain.second);

    // -- draw -----
    h->Draw();
    f_fit->Draw("same");

    TF1 *f_first_gaus = new TF1(Form("gauss_%s_0", h->GetName()), "gausn", fit_range_left, fit_range_right);
    f_first_gaus->SetParameters(result[0], result[1], result[2]);
    f_first_gaus->SetLineColor(kBlue);
    f_first_gaus->SetLineStyle(2);
    f_first_gaus->SetFillColor(kBlue);
    f_first_gaus->SetFillStyle(3003);
    f_first_gaus->Draw("same");

    for (Int_t i = 1; i < n_gauss; i++) {
        TF1 *f_single_gaus = new TF1(Form("gauss_%s_%d", h->GetName(), i), "gausn", fit_range_left, fit_range_right);
        f_single_gaus->SetParameters(result[i+4], result[1]+i*result[3], TMath::Sqrt(i+1)*result[4]);
        f_single_gaus->SetLineColor(kBlue);
        f_single_gaus->SetLineStyle(2);
        f_single_gaus->Draw("same");    
    }    

    c->Update();

    // //  +--------------+
    // //  | write result |
    // //  +--------------+
    // Int_t unused;
    // write_opg_result(".tmp_data.csv", run_num, ch, result, 0);
    // unused = system("python3 write.py test");
    // (void)unused;  // 変数を使わないことを明示

    return one_photon_gain;
}

// void wrapper(TGTab* tab, std::vector<Int_t> run_num_container, Int_t ch, Double_t led_vol) {
//     TCanvas *c = add_tab(tab, Form("led%.1f_ch%d", led_vol, ch));
//     c->Divide(1, 2);
//     TVirtualPad *c_fit = c->cd(1);
//     c_fit->Divide(3,1);
//     if (ch == 0 && led_vol == 3.4) c->Print(pdf_name + "["); // start

//     // -- fit and draw -----
//     std::vector<Double_t> x{}, x_err{}, y{}, y_err{};
//     for (Int_t cond = 0; cond < 3; cond++) {
//         std::pair<Double_t, Double_t> opg = analyze(run_num_container[cond], ch, c_fit, cond+1);
//         x.push_back(bac_pmt_vol[ch][cond]); 
//         x_err.push_back(0);
//         y.push_back(opg.first);
//         y_err.push_back(opg.second);    
//     }

//     // -- draw result -----
//     c->cd(2);
//     gPad->SetBottomMargin(0.15);
//     gPad->SetLeftMargin(0.15);
//     auto g = new TGraphErrors(3, &x[0], &y[0], &x_err[0], &y_err[0]);
//     g->SetTitle( Form("LED %.1f V, ch%d", led_vol, ch+1) );
//     g->GetXaxis()->SetTitle("PMT HV [V]");
//     g->GetYaxis()->SetTitle("Gain [ch]");
//     g->SetMarkerStyle(24);
//     // g->GetYaxis()->SetRangeUser(0, 35);
//     g->Draw("AP");
    
//     c->Print(pdf_name);
//     if (ch == 5 && led_vol == 3.6) c->Print(pdf_name + "]"); // end
// }

Int_t main(int argc, char** argv) {
    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
        return 1;
    }
    Int_t run_num = std::atoi(argv[1]);

    TApplication *theApp = new TApplication("App", &argc, argv);

    // -- create window -----
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    TGTab *tab = new TGTab(main, 1000, 800);

    // -- test -----
    TCanvas *c = add_tab(tab, "bac");
    c->Divide(2, 2);
    for (Int_t ch = 0; ch < max_bac_ch; ch++) analyze(run_num, ch, c, ch+1);

    // // // -- fit -----
    // // std::vector<Int_t> led_34{ 70107, 104, 101 };
    // // std::vector<Int_t> led_35{   108, 105, 113 };
    // // std::vector<Int_t> led_36{   109, 103, 100 };

    // // for (Int_t ch = 0; ch < max_bac_ch; ch++) {
    // //     wrapper(tab, led_34, ch, 3.4);
    // //     wrapper(tab, led_35, ch, 3.5);
    // //     wrapper(tab, led_36, ch, 3.6);
    // // }

    // -- add tab and draw window
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();
    
    theApp->Run();

    return 0;
}