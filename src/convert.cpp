#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TKey.h>
#include <TString.h>
#include <vector>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <root_file_path>" << std::endl;
        return 1;
    }

    // 開く ROOT ファイル
    TString file_path = argv[1];
    TFile input_file(file_path, "UPDATE");
    if (input_file.IsZombie()) {
        std::cerr << "Failed to open file: " << file_path << std::endl;
        return 1;
    }

    // TTree の準備
    TTree* tree = new TTree("tree", "Tree containing TF1 evaluation results");

    // 各ブランチ用の変数
    std::vector<std::vector<Double_t>> x_values(10);
    std::vector<std::vector<Double_t>> y_values(10);

    // ファイル内のすべての TF1 を処理
    TIter next(input_file.GetListOfKeys());
    TKey* key;

    Int_t index = 0;
    while ((key = (TKey*)next())) {
        if (TString(key->GetClassName()) != "TF1") {
            continue; // TF1 以外はスキップ
        }

        TF1* func = (TF1*)key->ReadObj();
        if (!func) {
            std::cerr << "Failed to read object: " << key->GetName() << std::endl;
            continue;
        } 

        TString obj_name = func->GetName();
        Double_t x_min = func->GetXmin();
        Double_t x_max = func->GetXmax();

        // 等間隔に分割する点
        Double_t num_points = 10000;
        Double_t step = (x_max - x_min) / (num_points - 1.0);

        // x, y の値を計算して保存
        std::vector<Double_t> x_vec, y_vec;
        for (Double_t i = 0.0; i < num_points; ++i) {
            Double_t x = x_min + i * step;
            Double_t y = func->Eval(x);
            x_vec.push_back(x);
            y_vec.push_back(y);
        }

        // ブランチ用の vector に追加
        x_values[index] = x_vec;
        y_values[index] = y_vec;
        tree->Branch(Form("%s_x", obj_name.Data()), &x_values[index]);
        tree->Branch(Form("%s_y", obj_name.Data()), &y_values[index]);
        index++;
    }

    // 読み込んだ ROOT ファイルに保存
    tree->Fill();
    tree->Write();
    input_file.Close();

    std::cout << "Tree saved to " << file_path << std::endl;

    return 0;
}
