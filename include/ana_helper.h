#ifndef ANA_HELPER_
#define ANA_HELPER_

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>

// ROOTヘッダー
#include <TCanvas.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>

namespace ana_helper {
    TCanvas* add_tab(TGTab *tab, const char* tabName);
    std::vector<Int_t> get_should_hit_seg(Int_t run_number);
}

#endif  // ANA_HELPER_