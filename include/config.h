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
    Int_t tdc_bin_num = 65536;
    Int_t adjust_tdc_bin_num = 65536; // for adjustment
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
    Double_t trig_counter_tdc_n_sigma = 5.0;
    Double_t trig_counter_adc_min_n_sigma = 3.0;
    Double_t trig_counter_adc_max_n_sigma = 4.0;
    

    // -- cherenkov counter -----
    Double_t default_cherenkov_tdc_gate_min = 105000;
    Double_t default_cherenkov_tdc_gate_max = 140000;
    
    Double_t cherenkov_tdc_n_sigma = 5.0;
    
    Bool_t log_flag = false;
    
    Double_t linear_fit_range_min = 0.0;
    Double_t linear_fit_range_max = 3500.0;

    // BAC one photon gain
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> bac_opg{
    //    HV    val,   err
        { 56, {{ 8.676, 0.015},
                {9.468, 0.016},
                {9.438, 0.016},
                {8.312, 0.018}} },
       
        { 57, {{11.106, 0.015},
               {11.841, 0.016},
               {11.699, 0.016},
               {10.639, 0.018}} },

        { 58, {{13.400, 0.015},
               {14.129, 0.016},
               {13.954, 0.016},
               {12.742, 0.018}} }       
    };

    // KVC one photon gain
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> kvc_thick_opg{
    //    HV     val,    err
        { 56, { {10.096, 0.087},    // board1 kvc seg1
                {10.088, 0.083},    // board1 kvc seg2
                {10.121, 0.092},    // board1 kvc seg3
                {10.059, 0.073},    // board1 kvc seg4

                { 9.259, 0.095},    // board2 kvc seg1
                { 9.628, 0.114},    // board2 kvc seg2
                { 9.613, 0.141},    // board2 kvc seg3
                { 9.090, 0.056},    // board2 kvc seg4

                {10.323, 0.091},    // board3 kvc seg1
                {10.846, 0.071},    // board3 kvc seg2
                {10.378, 0.099},    // board3 kvc seg3
                {10.515, 0.089},    // board3 kvc seg4

                { 9.511, 0.095},    // board4 kvc seg1
                { 9.657, 0.114},    // board4 kvc seg2
                { 9.616, 0.153},    // board4 kvc seg3
                { 9.429, 0.071}}},  // board4 kvc seg4

        { 57, { {13.028, 0.087},    // board1 kvc seg1
                {13.018, 0.083},    // board1 kvc seg2
                {13.060, 0.092},    // board1 kvc seg3
                {12.980, 0.073},    // board1 kvc seg4

                {11.983, 0.095},    // board2 kvc seg1
                {12.460, 0.114},    // board2 kvc seg2
                {12.440, 0.141},    // board2 kvc seg3
                {11.764, 0.056},    // board2 kvc seg4

                {13.314, 0.091},    // board3 kvc seg1
                {13.988, 0.071},    // board3 kvc seg2
                {13.385, 0.099},    // board3 kvc seg3
                {13.561, 0.089},    // board3 kvc seg4

                {12.777, 0.095},    // board4 kvc seg1
                {12.972, 0.114},    // board4 kvc seg2
                {12.917, 0.153},    // board4 kvc seg3
                {12.666, 0.071}}},  // board4 kvc seg4

        { 58, { {15.975, 0.087},    // board1 kvc seg1
                {15.962, 0.083},    // board1 kvc seg2
                {16.014, 0.092},    // board1 kvc seg3
                {15.916, 0.073},    // board1 kvc seg4

                {14.608, 0.095},    // board2 kvc seg1
                {15.189, 0.114},    // board2 kvc seg2
                {15.166, 0.141},    // board2 kvc seg3
                {14.341, 0.056},    // board2 kvc seg4

                {16.355, 0.091},    // board3 kvc seg1
                {17.183, 0.071},    // board3 kvc seg2
                {16.441, 0.099},    // board3 kvc seg3
                {16.658, 0.089},    // board3 kvc seg4

                {15.707, 0.095},    // board4 kvc seg1
                {15.947, 0.114},    // board4 kvc seg2
                {15.879, 0.153},    // board4 kvc seg3
                {15.571, 0.071}}}   // board4 kvc seg4

    };

    // KVC one photon gain
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> kvc_thin_opg{
    //    HV     val,    err
        { 56, { { 9.303, 0.050},    // up kvc seg1
                { 9.120, 0.040},    // up kvc seg2
                { 9.723, 0.047},    // up kvc seg3
                {10.490, 0.053},    // up kvc seg4

                { 9.336, 0.052},    // down kvc seg1
                { 9.345, 0.048},    // down kvc seg2
                { 9.325, 0.053},    // down kvc seg3
                { 9.457, 0.044}}},  // down kvc seg4

        { 57, { {11.962, 0.050},    // up kvc seg1
                {11.728, 0.040},    // up kvc seg2
                {12.502, 0.047},    // up kvc seg3
                {13.489, 0.053},    // up kvc seg4

                {12.004, 0.052},    // down kvc seg1
                {12.016, 0.048},    // down kvc seg2
                {11.990, 0.053},    // down kvc seg3
                {12.160, 0.044}}},  // down kvc seg4

        { 58, { {14.390, 0.050},    // up kvc seg1
                {14.108, 0.040},    // up kvc seg2
                {15.040, 0.047},    // up kvc seg3
                {16.227, 0.053},    // up kvc seg4

                {14.441, 0.052},    // down kvc seg1
                {14.455, 0.048},    // down kvc seg2
                {14.424, 0.053},    // down kvc seg3
                {14.629, 0.044}}}   // down kvc seg4
    };

    void bac_initialize() {
        adjust_adc_bin_num = 1024;
        sumadc_bin_num = 1024;
        npe_bin_num = 840;
    }

    void kvc_thin_initialize() {
        adjust_adc_bin_num = 1024;
        sumadc_bin_num = 512;
        npe_bin_num = 525;
        adjust_tdc_bin_num = 32768;
        log_flag = true;
        linear_fit_range_max = 1500.0;
    }
    
    void kvc_thick_initialize() {
        adjust_adc_bin_num = 1024;
        sumadc_bin_num = 512;
        npe_bin_num = 420;
        npe_max = 835.;
        adjust_tdc_bin_num = 32768;
        log_flag = true;
        linear_fit_range_max = 2000.0;
    }
    

private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
