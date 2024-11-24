# KEKAR2023

## メモ
`ana_hepler`の名前空間の中に共通の関数は入れようと思うが、たくさんのコードを1つのソースファイルに書くのは見にくい。そこでヘッダーファイルは共通だけど、ソースファイルはgeneralっぽそうなものは`ana_helper.cpp`に、トリガーカウンター関連は`trigger_counter.cpp`に分けて書くみたいなことをしている

解析上rowヒストグラムとトリガーに引っかかったヒストグラムの用意をする必要があって、online sum/offline sumなどに至ってはヒストグラムが増えすぎて可読性が落ちる。そこでrawとtrigのヒストグラムをペアにした構造体を使うことに(シンプルさには欠けるが可読性は上がる)。ヒストグラムをnewで作るので、一応デストラクタも書いてみたが、Canvasに描画して最後に確認するみたいな場合にはちゃんと消されてしまうので、見えなくなる。まぁ大規模なコードではないのでメモリリークは起きないと思ってその辺には目をつぶる。



# matplotlibcpp セットアップ (Ubuntu)

このドキュメントでは、Ubuntu環境で`matplotlibcpp`をセットアップしてC++プロジェクトで利用する方法を説明します。

---

## 必要な環境

1. **Python3** と開発ライブラリ
2. Python用ライブラリ:
   - `matplotlib`
   - `numpy`
3. **pybind11** (Python-C++バインディング用)
4. **matplotlibcpp.h** ヘッダファイル

---

## セットアップ手順

### 1. Python3と必要なライブラリをインストール

以下のコマンドでPython3と関連ライブラリをインストールします。

```bash
sudo apt update
sudo apt install python3 python3-pip python3-dev
```

次に、Pythonライブラリをインストールします。

```bash
pip3 install matplotlib numpy
```

`pybind11`もインストールします。

```bash
sudo apt install pybind11-dev
```

---

### 2. `matplotlibcpp.h` を取得

`matplotlibcpp`リポジトリをクローンし、必要なヘッダファイルをプロジェクトにコピーします。

```bash
git clone https://github.com/lava/matplotlib-cpp.git
cp matplotlib-cpp/matplotlibcpp.h /your/project/include/
```

`/your/project/include/` はプロジェクトのインクルードディレクトリに置き換えてください。

---

### 3. コードのコンパイルとリンク

#### **CMakeを使用する場合**

以下を`CMakeLists.txt`に追加します。

```cmake
cmake_minimum_required(VERSION 3.10)
project(MatplotlibCppExample)

find_package(Python3 COMPONENTS Development NumPy REQUIRED)
find_package(pybind11 REQUIRED)

include_directories(${Python3_INCLUDE_DIRS})
include_directories(/path/to/your/include)

add_executable(example main.cpp)
target_link_libraries(example ${Python3_LIBRARIES})
```

#### **手動でコンパイルする場合**

以下のコマンドでコンパイルします。

```bash
g++ -std=c++17 main.cpp -o example -I/path/to/your/include -lpython3.8
```

`-lpython3.8`はシステムのPythonバージョンに応じて変更してください（例: `-lpython3.10`）。

---

## サンプルコード

以下を`main.cpp`として保存してください。

```cpp
#include <matplotlibcpp.h>
#include <vector>

namespace plt = matplotlibcpp;

int main() {
    // データ
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {1, 4, 9, 16, 25};

    // プロット
    plt::plot(x, y, "r-");
    plt::title("サンプルプロット");
    plt::show();

    return 0;
}
```

次に、以下のコマンドでコンパイル＆実行します。

```bash
g++ -std=c++17 main.cpp -o example -I/path/to/your/include -lpython3.8
./example
```

---

## トラブルシューティング

### **1. Pythonバージョンの問題**
Pythonのバージョンを確認し、適切なバージョンをリンクしているか確認してください。

```bash
python3 --version
```

### **2. 必要なPythonライブラリが見つからない**
以下のコマンドでライブラリを再インストールしてください。

```bash
pip3 install matplotlib numpy
```

### **3. pybind11が見つからない**
以下のコマンドでインストールしてください。

```bash
pip3 install pybind11
```

---

これでUbuntu環境で`matplotlibcpp`を利用する準備が整います！
