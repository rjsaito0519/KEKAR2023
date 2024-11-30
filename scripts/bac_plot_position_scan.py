import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import seaborn as sns
import glob
import cv2
import uproot
import os

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


script_dir = os.path.dirname(os.path.abspath(__file__))
root_file_path = os.path.join(script_dir, "../results/root/bac_pos_scan_analysis.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")


def convert_x(x, x_min, width = 16):
    return (x - x_min)/width + 1/2

def convert_y(y, y_max, width = 18):
    return (y_max - y)/width + 1/2

def concat_tile(im_list_2d):
    return cv2.vconcat([cv2.hconcat(im_list_h) for im_list_h in im_list_2d])

x_pos = np.array([ -48, -32, -16, 0, 16, 32, 48 ])
# y_pos = np.array([ -54, -36, -18, 0, 18, 36, 54 ])
y_pos = np.array([ -36, -18, 0, 18, 36 ])

# -- data selection -----
data = {
    "eff":[],
    "npe":[]
}
for y in y_pos:
    eff_container = []
    npe_container = []
    for x in x_pos:
        indices = np.where((tree["run_num"] < 393) * (tree["pos_y"] == y) * (tree["pos_x"] == x))[0]
        eff_max = 0
        npe_max = 0
        for index in indices:
            eff = tree["n_hit"][index]/tree["n_trig"][index] * 100
            if eff_max < eff:
                eff_max = eff
            if npe_max < tree["offsum_npe_val"][index][0]:
               npe_max = tree["offsum_npe_val"][index][0] 
        eff_container.append(eff_max)
        npe_container.append(npe_max)
    data["eff"].append(eff_container)
    data["npe"].append(npe_container)



plt.figure(figsize=(10, 8))
df = pd.DataFrame(data["npe"], columns=x_pos, index = y_pos)
# -- add margin left and right ----- 
df[np.min(x_pos)-16] = [np.nan]*len(y_pos)
df[np.max(x_pos)+16] = [np.nan]*len(y_pos)
# -- add margin top and bottom ----- 
df.loc[np.min(y_pos)-18] = np.nan
df.loc[np.max(y_pos)+18] = np.nan

df = df.sort_index(axis=0, ascending=False)
df = df.sort_index(axis=1)
print(df)
plt.grid(which="major", alpha=0.3)
# sns.heatmap(df, cmap="viridis", annot=True, cbar=True, fmt='.3f', annot_kws={'fontsize': 16}, cbar_kws=dict(pad=-0.06, shrink=0.92))
sns.heatmap(df, cmap="cividis", annot=True, cbar=True, fmt='.3f', annot_kws={'fontsize': 16}, cbar_kws=dict(pad=-0.06, shrink=0.92))


# -- set radiator size ------
edge_left   = -115/2
edge_right  =  115/2
edge_top    =  115/2
edge_bottom = -115/2
# 左端
plt.vlines(convert_x(edge_left, np.min(df.columns)), convert_y(edge_top, np.max(df.index)), convert_y(edge_bottom, np.max(df.index)), color = "red", ls = "dashed", lw = 2, zorder = 1)
# 右端
plt.vlines(convert_x(edge_right, np.min(df.columns)), convert_y(edge_top, np.max(df.index)), convert_y(edge_bottom, np.max(df.index)), color = "red", ls = "dashed", lw = 2, zorder = 1)
# 上端
plt.hlines(convert_y(edge_top, np.max(df.index)), convert_x(edge_left, np.min(df.columns)), convert_x(edge_right, np.min(df.columns)), color = "red", ls = "dashed", lw = 2, zorder = 1)
# 下端
plt.hlines(convert_y(edge_bottom, np.max(df.index)), convert_x(edge_left, np.min(df.columns)), convert_x(edge_right, np.min(df.columns)), color = "red", ls = "dashed", lw = 2, zorder = 1)

plt.text(convert_x(0, np.min(df.columns)), convert_y(edge_top+3, np.max(df.index)), "BAC Efficiency (3-layer)", ha='center', va='bottom', zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
# plt.savefig("./img/SAC_Eff_corr2.jpg", dpi=600)
plt.subplots_adjust(left = 0.12, right = 0.98, top = 0.93, bottom = 0.1)
plt.show()