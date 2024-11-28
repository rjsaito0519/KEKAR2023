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
    Double_t trig_counter_tdc_n_sigma = 5.0;
    Double_t trig_counter_adc_min_n_sigma = 3.0;
    Double_t trig_counter_adc_max_n_sigma = 4.0;

    

    // -- cherenkov counter -----
    Double_t cherenkov_tdc_n_sigma = 5.0;

    // BAC one photon gain
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> bac_opg{
    //    HV    val,   err
        { 56, {{ 8.598, 0.015},     // ch1
               { 9.429, 0.016},     // ch2
               { 9.314, 0.016},     // ch3
               { 8.198, 0.018} }},  // ch4
        { 57, {{11.077, 0.015},
               {11.827, 0.016},
               {11.656, 0.016},
               {10.597, 0.018}} },
        { 58, {{13.417, 0.015},
               {14.140, 0.016},
               {13.992, 0.016},
               {12.766, 0.018}} }       
    };

    // KVC one photon gain
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> kvc_thick_opg{
    //    HV     val,    err
        { 56, { {10.050, 0.087},    // board1 kvc seg1
                {10.039, 0.083},    // board1 kvc seg2
                {10.087, 0.092},    // board1 kvc seg3
                { 9.996, 0.073},    // board1 kvc seg4

                { 8.962, 0.095},    // board2 kvc seg1
                { 9.508, 0.114},    // board2 kvc seg2
                { 9.486, 0.141},    // board2 kvc seg3
                { 8.711, 0.056},    // board2 kvc seg4

                {10.107, 0.091},    // board3 kvc seg1
                {10.827, 0.071},    // board3 kvc seg2
                {10.182, 0.099},    // board3 kvc seg3
                {10.371, 0.089},    // board3 kvc seg4

                { 9.417, 0.095},    // board4 kvc seg1
                { 9.634, 0.114},    // board4 kvc seg2
                { 9.573, 0.153},    // board4 kvc seg3
                { 9.294, 0.071}}},  // board4 kvc seg4

        { 57, { {13.011, 0.087},    // board1 kvc seg1
                {12.999, 0.083},    // board1 kvc seg2
                {13.047, 0.092},    // board1 kvc seg3
                {12.956, 0.073},    // board1 kvc seg4

                {11.866, 0.095},    // board2 kvc seg1
                {12.413, 0.114},    // board2 kvc seg2
                {12.391, 0.141},    // board2 kvc seg3
                {11.615, 0.056},    // board2 kvc seg4

                {13.264, 0.091},    // board3 kvc seg1
                {13.984, 0.071},    // board3 kvc seg2
                {13.339, 0.099},    // board3 kvc seg3
                {13.528, 0.089},    // board3 kvc seg4

                {12.748, 0.095},    // board4 kvc seg1
                {12.965, 0.114},    // board4 kvc seg2
                {12.903, 0.153},    // board4 kvc seg3
                {12.625, 0.071}}},  // board4 kvc seg4

        { 58, { {15.986, 0.087},    // board1 kvc seg1
                {15.974, 0.083},    // board1 kvc seg2
                {16.022, 0.092},    // board1 kvc seg3
                {15.931, 0.073},    // board1 kvc seg4

                {14.666, 0.095},    // board2 kvc seg1
                {15.213, 0.114},    // board2 kvc seg2
                {15.191, 0.141},    // board2 kvc seg3
                {14.415, 0.056},    // board2 kvc seg4

                {16.473, 0.091},    // board3 kvc seg1
                {17.193, 0.071},    // board3 kvc seg2
                {16.548, 0.099},    // board3 kvc seg3
                {16.737, 0.089},    // board3 kvc seg4

                {15.736, 0.095},    // board4 kvc seg1
                {15.954, 0.114},    // board4 kvc seg2
                {15.892, 0.153},    // board4 kvc seg3
                {15.614, 0.071}}}   // board4 kvc seg4
    };

    // KVC one photon gain
    std::unordered_map<Int_t, std::vector<std::pair<Double_t, Double_t>>> kvc_thin_opg{
    //    HV     val,    err
        { 56, { { 9.417, 0.050},    // up kvc seg1
                { 9.174, 0.040},    // up kvc seg2
                { 9.976, 0.047},    // up kvc seg3
                {10.998, 0.053},    // up kvc seg4

                { 9.460, 0.052},    // down kvc seg1
                { 9.472, 0.048},    // down kvc seg2
                { 9.446, 0.053},    // down kvc seg3
                { 9.622, 0.044}}},  // down kvc seg4

        { 57, { {11.978, 0.050},    // up kvc seg1
                {11.735, 0.040},    // up kvc seg2
                {12.537, 0.047},    // up kvc seg3
                {13.559, 0.053},    // up kvc seg4

                {12.022, 0.052},    // down kvc seg1
                {12.033, 0.048},    // down kvc seg2
                {12.007, 0.053},    // down kvc seg3
                {12.183, 0.044}}},  // down kvc seg4

        { 58, { {14.316, 0.050},    // up kvc seg1
                {14.073, 0.040},    // up kvc seg2
                {14.876, 0.047},    // up kvc seg3
                {15.898, 0.053},    // up kvc seg4

                {14.360, 0.052},    // down kvc seg1
                {14.372, 0.048},    // down kvc seg2
                {14.346, 0.053},    // down kvc seg3
                {14.522, 0.044}}}   // down kvc seg4

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
