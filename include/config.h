// config.h
#ifndef CONFIG_H
#define CONFIG_H

class Config {
public:
    static Config& getInstance() {
        static Config instance;
        return instance;
    }

    // チャンネル数の設定
    Int_t max_sac_ch = 8;
    Int_t max_kvc_ch = 4;
    Int_t max_bac_ch = 4;
    Int_t max_nhit_tdc = 16;

    // ADCのbin設定
    Int_t adc_bin_num = 4096;
    Int_t adjust_adc_bin_num = 4096; // for adjustment
    Double_t adc_min = 0.0;
    Double_t adc_max = 4096.0;

    // TDCのbin設定
    Int_t tdc_bin_num = 0x10000;
    Double_t tdc_min = 0.0;
    Double_t tdc_max = 2097152.0;

    // NPEのbin設定
    Int_t npe_bin_num = 2100;
    Double_t npe_min = -5.;
    Double_t npe_max = 415.;

    // online sum(pedestal引いたやつ)のbin設定
    Int_t sumadc_bin_num = 8192;
    Double_t sumadc_min = 0.0;
    Double_t sumadc_max = 8192.0;


    // -- for checking shower -----
    Int_t separate_run_num = 392;
    std::vector<Double_t> adc_min_condition1{  300.0, 300.0, 300.0, 300.0 };
    std::vector<Double_t> adc_min_condition2{ 1000.0, 600.0, 700.0, 600.0 };
    std::vector<std::pair<Double_t, Double_t>> tdc_gate_condition1{
        {115000, 135000},
        {115000, 135000},
        {115000, 135000},
        {115000, 135000}
    };
    std::vector<std::pair<Double_t, Double_t>> tdc_gate_condition2{
        {115000, 135000},
        {115000, 135000},
        {115000, 135000},
        {115000, 135000}
    };

    // -- trigger counter -----
    Double_t trig_counter_adc_target_val_ratio = 0.01;
    Double_t trig_counter_tdc_n_sigma =5.0;

    // -- cherenkov counter -----
    Double_t cherenkov_tdc_n_sigma =5.0;

    // BAC one photon gain (KEKAR result)
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> bac_opg{
    //    HV    val,   err
        { 56, {{ 8.850, 0.099}, // run00239
               { 9.570, 0.049},
               { 9.777, 0.044},
               { 8.571, 0.084} }},
        { 57, {{11.329, 0.067}, // run00238
               {11.968, 0.036},
               {12.120, 0.037},
               {10.970, 0.050}} },
        { 58, {{13.669, 0.073}, // run00240
               {14.280, 0.037},
               {14.455, 0.041},
               {13.138, 0.056}} }
    };


    void bac_initialize() {
        adjust_adc_bin_num = 1024;
        sumadc_bin_num = 1024;
        npe_bin_num = 840;
    }

private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
