void multiplicity(TString file_path){
    // -- parameter ------------------------------------------------------------------
    Int_t maxTDChit   = 16;
    Int_t maxKVCSUMch = 4;

    // -- load file ------------------------------------------------------------------
    auto *f = new TFile( file_path.Data() );
    TTreeReader reader("tree", f);

    // -- prepare reader and histogram -----------------------------------------------------------------------------------
    TTreeReaderArray<Double_t> kvcsumt(reader, "kvcsumt");
    auto *hKVCMulti = new TH1I("kvc_multi", "KVC multiplicity", 5, -0.5, 4.5);
    auto *hKVCSeg   = new TH1I("kvc_segment", "KVC segment", maxKVCSUMch, -0.5, maxKVCSUMch-0.5);
    
    // -- set TDC gate -----------------------------------------------------------------------------------
    Double_t MinTDCKVC = 0;
    Double_t MaxTDCKVC = 0x200000;

    // --- fill event -----------------------------------------------------------------------------------
    Bool_t hitFlag = false;
    Int_t  tmp_multi = 0;
    reader.Restart();
    while (reader.Next() ){
        tmp_multi = 0;
        for (Int_t ch = 0; ch < maxKVCSUMch; ch++) {
            hitFlag = false;
            for (Int_t n_hit = 0; n_hit < maxTDChit; n_hit++) {
                if ( MinTDCKVC < kvcsumt[ch*maxTDChit + n_hit] && kvcsumt[ch*maxTDChit + n_hit] < MaxTDCKVC ) hitFlag = true;
            }
            if (hitFlag) {
                tmp_multi++;
                hKVCSeg->Fill(ch);
            }
        }
        hKVCMulti->Fill(tmp_multi);
    }

    // --- draw figure -----------------------------------------------------------------------------------
    TCanvas *c = new TCanvas("", "", 600, 600);
    c->Divide(1, 2);
    c->cd(1);
    hKVCMulti->SetTitle("KVC multiplicity;multiplicity;");
    hKVCMulti->Draw();
    c->cd(2);
    hKVCSeg->SetTitle("KVC hit segment;ch;");
    hKVCSeg->Draw();

}