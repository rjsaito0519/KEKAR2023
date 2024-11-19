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
    Int_t max_sac_ch = 6;
    Int_t max_kvc_ch = 1;
    Int_t max_bac_ch = 4;
    Int_t max_nhit_tdc = 16;

    // ADCのbin設定
    Int_t adc_bin_num = 4096;
    Int_t adc_min = 0;
    Int_t adc_max = 4096;

    // ADCのbin設定
    Int_t tdc_bin_num = 0x10000;
    Int_t tdc_min = 0;
    Int_t tdc_max = 0x200000;

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


    // shower event cut condition (use KVC sum)
    Double_t shower_tdc_min = 115000.0;
    Double_t shower_tdc_max = 135000.0;
    Double_t shower_adc_min = 300.0;
    

    // // 必要に応じて初期化メソッドを追加
    // void initialize(double th, int run, int samples) {
    //     threshold = th;
    //     currentRunNumber = run;
    //     sampleCount = samples;
    // }

private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
