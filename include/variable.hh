#ifndef VARIABLE_H_
#define VARIABLE_H_

// チャンネル数の設定
static const Int_t max_sac_ch = 6;
static const Int_t max_kvc_ch = 1;
static const Int_t max_bac_ch = 4;
static const Int_t max_nhit_tdc = 16;

// ADCのbin設定
static const Int_t adc_bin_num = 4096;
static const Int_t adc_min = 0;
static const Int_t adc_max = 4096;

// TDCのbin設定
static const Int_t tdc_bin_num = 1500;
static const Int_t tdc_min = 0;
static const Int_t tdc_max = 1500;

// NPEのbin設定
static const Int_t npe_bin_num = 2100;
static const Double_t npe_min = -5.;
static const Double_t npe_max = 415.;

// online sum(pedestal引いたやつ)のbin設定
static const Int_t onsum_bin_num = 8192;
static const Int_t onsum_min = 0;
static const Int_t onsum_max = 2048;

// triggerの範囲の設定
static const Int_t tdc_gate_n_sigma = 5;
static const Int_t adc_gate_n_sigma = 5;

// pedestalの値
// --- SAC ----------------------------------------------------
// static const Double_t pedestal_SACSUM_pos44  = 226.444;
// static const Double_t pedestal_SACSUM_pos57  = 247.242;
// static const Double_t pedestal_SACSUM_pos79  = 235.316;
// static const Double_t pedestal_SACSUM_pos285 = 212.811;
static const Double_t sacsum_pedestal_pos[4] = {
    226.444,
    247.242,
    235.316,
    212.811
};

// --- BAC ----------------------------------------------------
static const Double_t BACindiv_pedestal_pos[4][4] = {
    //    44,      57,      79,     285
    {355.215, 211.971, 170.746, 183.053},
    {355.476, 211.262, 170.455, 181.924},
    {361.448, 218.378, 175.134, 185.723},
    {347.513, 206.285, 171.014, 174.982}
};
// static const Double_t pedestal_BACSUM_pos44  = 268.988;
// static const Double_t pedestal_BACSUM_pos57  = 269.095;
// static const Double_t pedestal_BACSUM_pos79  = 277.837;
// static const Double_t pedestal_BACSUM_pos285 = 273.013;
static const Double_t BACsum_pedestal_pos[4] = {
    268.988,
    269.095,
    277.837,
    273.013
};

// --- KVC ----------------------------------------------------
// static const Double_t pedestal_KVCSUM_pos44  = 226.444;
// static const Double_t pedestal_KVCSUM_pos57  = 247.242;
// static const Double_t pedestal_KVCSUM_pos79  = 235.316;
// static const Double_t pedestal_KVCSUM_pos285 = 212.811;
static const Double_t KACsum_pedestal_pos[4] = {
    160.649,
    161.961,
    162.151,
    151.752
};


// one photon gainの設定
// --- SAC ----------------------------------------------------
static const Double_t sac_one_photon_gain[6][3] = {
    {11.2342,  15.4303,  20.7793},  // 3.6 V
    {18.9387,  24.5975,  28.8386},  // 3.6 V
    { 7.03288,  9.85967, 10.3469},  // 3.6 V
    { 8.29301, 10.6327,  12.1416},  // 3.5 V 
    { 7.46516, 13.2924,  18.2573},  // 3.4 V
    {10.1753,  13.8259,  19.5989}   // 3.5 V
};
// static const Double_t SACSUM_one_photon_gain[3] = { 12.5889, 16.2106, 19.7347};

// --- BAC ----------------------------------------------------
static const Double_t BAC_one_photon_gain[4][3] = {
    //   56V,     57V,     58V    
    {11.2307, 14.9373, 18.6440}, 
    {9.01049, 12.6678, 16.3251}, 
    {7.90286, 11.2044, 14.5060}, 
    {9.58750, 11.8408, 14.0941}
};

// -- SAC -----
static const Double_t sac_pmt_vol[max_sac_ch][3] = {
    {2000, 2100, 2200},
    {1861, 1927, 1993},
    {1783, 1835, 1886},
    {1866, 1932, 1999},
    {1817, 1958, 2098},
    {1937, 2026, 2115}
};

#endif  // VARIABLE_H_