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


class KVC_thin(pos_scan_tool.pos_scan):
    def __init__(self, root_file_path):
        super().__init__(root_file_path)

    def plot(self, key, cbar_range = (np.nan, np.nan), img_type = "png"):
        x_pos = np.array([ -48, -32, -16, 0, 16, 32 ])
        # y_pos = np.array([ -54, -36, -18, 0, 18, 36, 54 ])
        y_pos = np.array([ -36, -18, 0, 18, 36 ])

        # -- data selection -----
        data = {
            "eff_val": [],
            "eff_err": [],
            "sum_npe0_val": [],
            "sum_npe0_err": [],
            "sum_npe1_val": [],
            "sum_npe1_err": [],
            "sum_npe2_val": [],
            "sum_npe2_err": [],
            "sum_npe3_val": [],
            "sum_npe3_err": [],
        }

        for y in y_pos:
            eff_container = []
            sum_npe_container = [[], [], [], []]
            for x in x_pos:
                indices = np.where((self.tree["pos_y"] == y) * (self.tree["pos_x"] == x))[0]
                # 検出効率が一番よかったrunを採用することにする。個々の最大値を取るのではなく、基準を検出効率にしてrunをそろえる
                eff_max = pos_scan_tool.pair(0, 0)
                sum_npe = [pos_scan_tool.pair(0, 0), pos_scan_tool.pair(0, 0), pos_scan_tool.pair(0, 0), pos_scan_tool.pair(0, 0)]
                for index in indices:
                    eff = self.tree["n_hit"][index]/self.tree["n_trig"][index]
                    if eff_max.val < eff:
                        eff_max.val = eff
                        eff_max.err = np.sqrt(eff*(1-eff)/self.tree["n_trig"][index])
                        for ch in range(4):
                            sum_npe[ch].val = self.tree["onsum_npe_val"][index][ch]
                            sum_npe[ch].err = self.tree["onsum_npe_err"][index][ch]
                eff_container.append(eff_max)
                for ch in range(4):
                    sum_npe_container[ch].append(sum_npe[ch])
            data["eff_val"].append([it.val*100 for it in eff_container])
            data["eff_err"].append([it.err*100 for it in eff_container])
            for ch in range(4):
                data[f"sum_npe{ch}_val"].append([it.val if it.val > 5 else np.nan for it in sum_npe_container[ch]])
                data[f"sum_npe{ch}_err"].append([it.err for it in sum_npe_container[ch]])

        # DataFrameの作成
        df_val = pd.DataFrame(data[f"{key}_val"], columns=x_pos, index=y_pos)
        df_err = pd.DataFrame(data[f"{key}_err"], columns=x_pos, index=y_pos)
        if np.isnan(cbar_range[0]):
            cbar_range = (np.min(df_val.values), np.max(df_val.values))

        # # -- 上下1列ずつ表示しない -----
        # df_val.loc[np.min(y_pos)] = np.nan
        # df_val.loc[np.max(y_pos)] = np.nan
        # df_err.loc[np.min(y_pos)] = np.nan
        # df_err.loc[np.max(y_pos)] = np.nan

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
            title = "KVC Efficiency (1 cm)"
            img_save_path = os.path.join(self.script_dir, f"../results/img/kvc/kvc_1cm_eff.{img_type}")
        else:
            color_map = "cividis"
            title = "KVC seg{} NPE (1 cm)".format(int(key[-1])+1)
            img_save_path = os.path.join(self.script_dir, "../results/img/kvc/kvc_1cm_seg{}_sumnpe.{}".format(int(key[-1])+1, img_type))
            
        os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
        
        plt.figure(figsize=(10, 8))
        plt.grid(which="major", alpha=0.3)
        # sns.heatmap(df_val, cmap=color_map, annot=annot, vmin=cbar_range[0], vmax=cbar_range[1], cbar=True, fmt='', annot_kws={'fontsize': 18}, cbar_kws=dict(pad=-0.08, shrink=0.78))
        sns.heatmap(df_val, cmap=color_map, annot=annot, vmin=cbar_range[0], vmax=cbar_range[1], cbar=True, fmt='', annot_kws={'fontsize': 18}, cbar_kws=dict(pad=-0.1, shrink=0.715))
        

        # -- set radiator size ------
        kvc_height = 120
        kvc_seg_width = 26
        edge_left   = -kvc_seg_width*2-kvc_seg_width/2
        edge_right  =  kvc_seg_width+kvc_seg_width/2
        edge_top    =  kvc_height/2
        edge_bottom = -kvc_height/2
        for i in range(5):
            plt.vlines(self.convert_x(edge_left+i*kvc_seg_width, np.min(df_val.columns)), self.convert_y(edge_top, np.max(df_val.index)), self.convert_y(edge_bottom, np.max(df_val.index)), color = "red", ls = "dashed", lw = 1.5, zorder = 1)
        # 上端
        plt.hlines(self.convert_y(edge_top, np.max(df_val.index)), self.convert_x(edge_left, np.min(df_val.columns)), self.convert_x(edge_right, np.min(df_val.columns)), color = "red", ls = "dashed", lw = 1.5, zorder = 1)
        # 下端
        plt.hlines(self.convert_y(edge_bottom, np.max(df_val.index)), self.convert_x(edge_left, np.min(df_val.columns)), self.convert_x(edge_right, np.min(df_val.columns)), color = "red", ls = "dashed", lw = 1.5, zorder = 1)

        # plt.text(self.convert_x((edge_left+edge_right)/2.0, np.min(df_val.columns)), self.convert_y(edge_top+6, np.max(df_val.index)), title, ha='center', va='bottom', zorder = 1)
        plt.text(self.convert_x((edge_left+edge_right)/2.0, np.min(df_val.columns)), self.convert_y(edge_top, np.max(df_val.index)), title, ha='center', va='bottom', zorder = 1)

        for ch in range(4):
            plt.text(
                self.convert_x(edge_left+(ch+0.5)*kvc_seg_width, np.min(df_val.columns)),
                self.convert_y(edge_bottom, np.max(df_val.index)), 
                f"seg. {ch+1}", ha='center', va='bottom', zorder = 1)


        plt.xlabel("x position [mm]")
        plt.ylabel("y position [mm]")
        plt.subplots_adjust(left = 0.1, right = 1.0, top = 0.96, bottom = 0.1)
        plt.savefig(img_save_path,  format=img_type, bbox_inches='tight', dpi=600, transparent=True)
        plt.show()
        
        return df_val, annot
    
    def plot_summary(self, df_val, annot, cbar_range, img_type = "png"):
        title = "KVCSUM NPE (1 cm)"
        img_save_path = os.path.join(self.script_dir, f"../results/img/kvc/kvc_1cm_sumnpe.{img_type}")

        plt.figure(figsize=(10, 8))
        plt.grid(which="major", alpha=0.3)
        # sns.heatmap(df_val, cmap="cividis", annot=annot, vmin=cbar_range[0], vmax=cbar_range[1], cbar=True, fmt='', annot_kws={'fontsize': 18}, cbar_kws=dict(pad=-0.08, shrink=0.78))
        sns.heatmap(df_val, cmap="cividis", annot=annot, vmin=cbar_range[0], vmax=cbar_range[1], cbar=True, fmt='', annot_kws={'fontsize': 18}, cbar_kws=dict(pad=-0.1, shrink=0.715))
        

        # -- set radiator size ------
        kvc_height = 120
        kvc_seg_width = 26
        edge_left   = -kvc_seg_width*2-kvc_seg_width/2
        edge_right  =  kvc_seg_width+kvc_seg_width/2
        edge_top    =  kvc_height/2
        edge_bottom = -kvc_height/2
        for i in range(5):
            plt.vlines(self.convert_x(edge_left+i*kvc_seg_width, np.min(df_val.columns)), self.convert_y(edge_top, np.max(df_val.index)), self.convert_y(edge_bottom, np.max(df_val.index)), color = "red", ls = "dashed", lw = 1.5, zorder = 1)
        # 上端
        plt.hlines(self.convert_y(edge_top, np.max(df_val.index)), self.convert_x(edge_left, np.min(df_val.columns)), self.convert_x(edge_right, np.min(df_val.columns)), color = "red", ls = "dashed", lw = 1.5, zorder = 1)
        # 下端
        plt.hlines(self.convert_y(edge_bottom, np.max(df_val.index)), self.convert_x(edge_left, np.min(df_val.columns)), self.convert_x(edge_right, np.min(df_val.columns)), color = "red", ls = "dashed", lw = 1.5, zorder = 1)

        # plt.text(self.convert_x((edge_left+edge_right)/2.0, np.min(df_val.columns)), self.convert_y(edge_top+6, np.max(df_val.index)), title, ha='center', va='bottom', zorder = 1)
        plt.text(self.convert_x((edge_left+edge_right)/2.0, np.min(df_val.columns)), self.convert_y(edge_top, np.max(df_val.index)), title, ha='center', va='bottom', zorder = 1)

        for ch in range(4):
            plt.text(
                self.convert_x(edge_left+(ch+0.5)*kvc_seg_width, np.min(df_val.columns)),
                self.convert_y(edge_bottom, np.max(df_val.index)), 
                f"seg. {ch+1}", ha='center', va='bottom', zorder = 1)


        plt.xlabel("x position [mm]")
        plt.ylabel("y position [mm]")
        plt.subplots_adjust(left = 0.1, right = 1.0, top = 0.96, bottom = 0.1)
        plt.savefig(img_save_path, format=img_type, bbox_inches='tight', dpi=600, transparent=True)
        plt.show()

if __name__ == '__main__':
    img_type = "pdf"

    kvc = KVC_thin("../results/root/kvc_thin_pos_scan_analysis.root")
    summary_df_val, summary_annot = kvc.plot("eff", cbar_range=(91, 100), img_type = img_type)
    
    for ch in range(4):
        # kvc.plot(f"sum_npe{ch}", cbar_range=(37, 81)) # use a, b
        tmp_df_val, tmp_annot = kvc.plot(f"sum_npe{ch}", cbar_range=(30, 59)) # use average one photon gain
        if ch == 0:
            summary_df_val[-48] = tmp_df_val[-48]
            summary_annot[-48]  = tmp_annot[-48]
        elif ch == 1:
            summary_df_val[-32] = tmp_df_val[-32]
            summary_annot[-32]  = tmp_annot[-32]
            summary_df_val[-16] = tmp_df_val[-16]
            summary_annot[-16]  = tmp_annot[-16]
        elif ch == 2:
            summary_df_val[0] = tmp_df_val[0]
            summary_annot[0]  = tmp_annot[0]
        elif ch == 3:
            summary_df_val[16] = tmp_df_val[16]
            summary_annot[16]  = tmp_annot[16]        
            summary_df_val[32] = tmp_df_val[32]
            summary_annot[32]  = tmp_annot[32]

    print(summary_annot)

    kvc.plot_summary(summary_df_val, summary_annot, (30, 59), img_type = img_type)