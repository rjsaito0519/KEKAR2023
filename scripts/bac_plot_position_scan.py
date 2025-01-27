import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import seaborn as sns
import glob
import cv2
import uproot
import os

import pos_scan_tool

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 28
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.axisbelow'] = True
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ


class BAC(pos_scan_tool.pos_scan):
    def __init__(self, root_file_path):
        super().__init__(root_file_path)

    def plot(self, key, cbar_range = (np.nan, np.nan), is_2layer = False, img_type = "pdf"):
        x_pos = np.array([ -48, -32, -16, 0, 16, 32, 48 ])
        y_pos = np.array([ -36, -18, 0, 18, 36 ])

        # -- data selection -----
        data = {
            "eff_val": [],
            "eff_err": [],
            "sum_npe_val": [],
            "sum_npe_err": [],
            "indiv_npe0_val": [],
            "indiv_npe0_err": [],
            "indiv_npe1_val": [],
            "indiv_npe1_err": [],
            "indiv_npe2_val": [],
            "indiv_npe2_err": [],
            "indiv_npe3_val": [],
            "indiv_npe3_err": [],
        }
        run_min, run_max = 0, 392
        if is_2layer:
            run_min, run_max = 392, np.inf 

        for y in y_pos:
            eff_container = []
            sum_npe_container = []
            indiv_npe_container = [[], [], [], []]
            for x in x_pos:
                indices = np.where((run_min < self.tree["run_num"]) * (self.tree["run_num"] <= run_max) * (self.tree["pos_y"] == y) * (self.tree["pos_x"] == x))[0]
                # 検出効率が一番よかったrunを採用することにする。個々の最大値を取るのではなく、基準を検出効率にしてrunをそろえる
                eff_max = pos_scan_tool.pair(0, 0)
                sum_npe = pos_scan_tool.pair(0, 0)
                indiv_npe = [pos_scan_tool.pair(0, 0), pos_scan_tool.pair(0, 0), pos_scan_tool.pair(0, 0), pos_scan_tool.pair(0, 0)]
                for index in indices:
                    eff = self.tree["n_hit"][index]/self.tree["n_trig"][index]
                    if eff_max.val < eff:
                        eff_max.val = eff
                        eff_max.err = np.sqrt(eff*(1-eff)/self.tree["n_trig"][index])
                        sum_npe.val = self.tree["onsum_npe_val"][index][0]
                        sum_npe.err = self.tree["onsum_npe_err"][index][0]
                        for ch in range(4):
                            indiv_npe[ch].val = self.tree["indiv_npe_val"][index][ch]
                            indiv_npe[ch].err = self.tree["indiv_npe_err"][index][ch]
                eff_container.append(eff_max)
                sum_npe_container.append(sum_npe)
                for ch in range(4):
                    indiv_npe_container[ch].append(indiv_npe[ch])
            data["eff_val"].append([it.val*100 for it in eff_container])
            data["eff_err"].append([it.err*100 for it in eff_container])
            data["sum_npe_val"].append([it.val for it in sum_npe_container])
            data["sum_npe_err"].append([it.err for it in sum_npe_container])
            for ch in range(4):
                data[f"indiv_npe{ch}_val"].append([it.val for it in indiv_npe_container[ch]])
                data[f"indiv_npe{ch}_err"].append([it.err for it in indiv_npe_container[ch]])

        # DataFrameの作成
        df_val = pd.DataFrame(data[f"{key}_val"], columns=x_pos, index=y_pos)
        df_err = pd.DataFrame(data[f"{key}_err"], columns=x_pos, index=y_pos)
        if np.isnan(cbar_range[0]):
            cbar_range = (np.min(df_val.values), np.max(df_val.values))

        # -- add margin left and right ----- 
        df_val[np.min(x_pos)-16] = [np.nan]*len(y_pos)
        df_val[np.max(x_pos)+16] = [np.nan]*len(y_pos)
        df_err[np.min(x_pos)-16] = [np.nan]*len(y_pos)
        df_err[np.max(x_pos)+16] = [np.nan]*len(y_pos)

        # -- add margin top and bottom ----- 
        df_val.loc[np.min(y_pos)-18] = np.nan
        df_val.loc[np.max(y_pos)+18] = np.nan
        df_err.loc[np.min(y_pos)-18] = np.nan
        df_err.loc[np.max(y_pos)+18] = np.nan

        # 両方の DataFrame をソート
        df_val = df_val.sort_index(axis=0, ascending=False)
        df_val = df_val.sort_index(axis=1)
        df_err = df_err.sort_index(axis=0, ascending=False)
        df_err = df_err.sort_index(axis=1)

        # annot を作成する
        annot = df_val.copy()  # 元の DataFrame をコピー

        # 各セルごとに値とエラーを組み合わせる
        for row_idx, col_idx in np.ndindex(df_val.shape):
            val = df_val.iloc[row_idx, col_idx]
            err = df_err.iloc[row_idx, col_idx]
            if pd.notna(val) and pd.notna(err):
                annot.iloc[row_idx, col_idx] = f" {val:.2f}\n" + r"$\pm$" + f"{err:.2f}"
            else:
                annot.iloc[row_idx, col_idx] = ""

        # 確認
        print(annot)

        if key == "eff":
            color_map = "viridis"
            title = "BAC Efficiency ({}-layer)".format(2 if is_2layer else 3)
            img_save_path = os.path.join(self.script_dir, "../results/img/bac/bac_{}layer_eff.{}".format(2 if is_2layer else 3, img_type))
        elif key == "sum_npe":
            color_map = "cividis"
            title = "BACSUM NPE ({}-layer)".format(2 if is_2layer else 3)
            img_save_path = os.path.join(self.script_dir, "../results/img/bac/bac_{}layer_sumnpe.{}".format(2 if is_2layer else 3, img_type))
        else:
            color_map = "cividis"
            title = "BAC ch{} NPE ({}-layer)".format(int(key[-1])+1, 2 if is_2layer else 3)
            img_save_path = os.path.join(self.script_dir, "../results/img/bac/bac_{}layer_ch{}npe.{}".format(2 if is_2layer else 3, int(key[-1])+1, img_type))
            
        os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
        
        plt.figure(figsize=(10, 8))
        plt.grid(which="major", alpha=0.3)
        sns.heatmap(df_val, cmap=color_map, annot=annot, vmin=cbar_range[0], vmax=cbar_range[1], cbar=True, fmt='', annot_kws={'fontsize': 18}, cbar_kws=dict(pad=-0.06, shrink=0.92))

        # -- set radiator size ------
        edge_left   = -115/2
        edge_right  =  115/2
        edge_top    =  115/2
        edge_bottom = -115/2
        # 左端
        plt.vlines(self.convert_x(edge_left, np.min(df_val.columns)), self.convert_y(edge_top, np.max(df_val.index)), self.convert_y(edge_bottom, np.max(df_val.index)), color = "red", ls = "dashed", lw = 2, zorder = 1)
        # 右端
        plt.vlines(self.convert_x(edge_right, np.min(df_val.columns)), self.convert_y(edge_top, np.max(df_val.index)), self.convert_y(edge_bottom, np.max(df_val.index)), color = "red", ls = "dashed", lw = 2, zorder = 1)
        # 上端
        plt.hlines(self.convert_y(edge_top, np.max(df_val.index)), self.convert_x(edge_left, np.min(df_val.columns)), self.convert_x(edge_right, np.min(df_val.columns)), color = "red", ls = "dashed", lw = 2, zorder = 1)
        # 下端
        plt.hlines(self.convert_y(edge_bottom, np.max(df_val.index)), self.convert_x(edge_left, np.min(df_val.columns)), self.convert_x(edge_right, np.min(df_val.columns)), color = "red", ls = "dashed", lw = 2, zorder = 1)

        plt.text(self.convert_x(0, np.min(df_val.columns)), self.convert_y(edge_top+3, np.max(df_val.index)), title, ha='center', va='bottom', zorder = 1)
        plt.xlabel("x position [mm]")
        plt.ylabel("y position [mm]")
        plt.subplots_adjust(left = 0.1, right = 1.0, top = 0.96, bottom = 0.1)
        plt.savefig(img_save_path,  format=img_type, bbox_inches='tight', dpi=600, transparent=True)
        plt.show()

if __name__ == '__main__':
    bac = BAC("../results/root/bac_pos_scan_analysis.root")
    bac.plot("eff", cbar_range=(97, 100), is_2layer=True, img_type="png")
    bac.plot("eff", cbar_range=(97, 100), is_2layer=False, img_type="png")
    bac.plot("sum_npe", cbar_range=(43, 93), is_2layer=True, img_type="png")
    bac.plot("sum_npe", cbar_range=(43, 93), is_2layer=False, img_type="png")
    
    for ch in range(4):
        bac.plot(f"indiv_npe{ch}", cbar_range=(3, 34), is_2layer=True)
        bac.plot(f"indiv_npe{ch}", cbar_range=(3, 34), is_2layer=False)
        
    
    