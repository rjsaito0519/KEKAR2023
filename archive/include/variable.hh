#ifndef VARIABLE_H_
#define VARIABLE_H_

// チャンネル数の設定
static const Int_t maxSACch = 8;
static const Int_t maxKVCch = 4;
static const Int_t maxKVCSUMch = 4;
static const Int_t maxBACch = 4;
static const Int_t maxTDChit = 16;

// ADCのbin設定
static const Int_t adc_bin_num = 4096;
static const Int_t adc_min = 0;
static const Int_t adc_max = 4096;

// TDCのbin設定
static const Int_t tdc_bin_num = 0x10000;
static const Int_t tdc_min = 0;
static const Int_t tdc_max = 0x200000;
static const Double_t tdc_ticks_size = (tdc_max-tdc_min)/tdc_bin_num;

// TDCをfillするときのthreshold(SACでありえないところにピークがあったりする)
static const Double_t TdcCutThreshold = 100000;


// triggerのTDCのfill範囲とbin数
static const Double_t TdcRangeT1[3] = {150, 129000, 132000};
static const Double_t TdcRangeT2[3] = {450, 125000, 134000};
static const Double_t TdcRangeT3[3] = {450, 125000, 134000};
static const Double_t TdcRangeT4[3] = {450, 123000, 132000};
static const Double_t TdcRangeSACSUM[3] = {500, 100000, 150000};
static const Double_t TdcRangeBACSUM[3] = {500, 100000, 150000};
static const Double_t TdcRangeKVCSUM[3] = {500, 100000, 150000};

// TDCのgate設定
static const Double_t TdcGateT1[2] = {120000, 140000};
static const Double_t TdcGateT2[2] = {110000, 150000};
static const Double_t TdcGateT3[2] = {120000, 140000};
static const Double_t TdcGateT4[2] = {110000, 140000};
static const Double_t TdcGateSACSUM[2] = {110000, 130000};
static const Double_t TdcGateBACSUM[2] = {120000, 140000};
static const Double_t TdcGateKVCSUM[2] = {110000, 140000};

// NPEのbin設定
// static const Int_t npe_bin_num = 2100;
static const Int_t npe_bin_num = 840;
static const Double_t npe_min = -5.;
static const Double_t npe_max = 415.;

// triggerの範囲の設定
static const Int_t tdc_n_sigma = 20;
static const Int_t adc_n_sigma = 7;

// pedestalの値
// --- SAC ----------------------------------------------------
// static const Double_t pedestal_SACSUM_pos44  = 226.444;
// static const Double_t pedestal_SACSUM_pos57  = 247.242;
// static const Double_t pedestal_SACSUM_pos79  = 235.316;
// static const Double_t pedestal_SACSUM_pos285 = 212.811;
static const Double_t SACsum_pedestal_pos[4] = {
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
static const Double_t SAC_one_photon_gain[maxSACch][3] = {
    {22.3218, 19.5891, 17.3434},
    {20.0320, 18.2662, 17.1973},
    {17.2782, 15.6305, 14.0576},
    {16.1369, 14.0492, 12.7709},
    {18.8409, 17.0543, 14.9761},
    {16.4578, 14.2855, 12.5269},
    {13.2460, 12.7492, 12.0413},
    {14.3252, 12.2378, 10.8038}
};
// static const Double_t SACSUM_one_photon_gain[3] = { 12.5889, 16.2106, 19.7347};

// --- BAC ----------------------------------------------------
// static const Double_t BAC_one_photon_gain[4][3] = {
//     //   56V,     57V,     58V    
//     {11.2307, 14.9373, 18.6440}, 
//     {9.01049, 12.6678, 16.3251}, 
//     {7.90286, 11.2044, 14.5060}, 
//     {9.58750, 11.8408, 14.0941}
// };
static const Double_t BAC_one_photon_gain = 15.088540;

static const Double_t KVC_one_photon_gain = 14.181332;




#endif  // VARIABLE_H_