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
