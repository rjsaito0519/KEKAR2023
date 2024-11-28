import os
import numpy as np
import matplotlib.pyplot as plt
import uproot

import opg_tool

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 28
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = False
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
root_file_path = os.path.join(script_dir, "../results/root/bac_opg.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

mppc_map = {
    138: 1,
    142: 2,
    144: 3,
    147: 4,
    149: 5,
    150: 5,
    153: 6,
    154: 6,
    157: 7,
    158: 7,
    160: 8,
    162: 9,
    164: 10,
    167: 11,
    169: 12,
    174: 13,
    176: 14,
    178: 15,
    180: 16
}


HV_map = {
    57 : 238,
    56 : 239,
    58 : 240
}

ch1 = opg_tool.data_summarize(tree, 0, mppc_map)
ch2 = opg_tool.data_summarize(tree, 1, mppc_map)
ch3 = opg_tool.data_summarize(tree, 2, mppc_map)
ch4 = opg_tool.data_summarize(tree, 3, mppc_map)

indiv_data = [ch1, ch2, ch3, ch4]
use_ch = 0

for HV in [56, 57, 58]:
    offset = []
    for ch in range(4):
        for i in range(len(tree["run_num"])):
            if tree["run_num"][i] == HV_map[HV] and tree["ch"][i] == ch:
                offset.append(tree["result_val"][i][3] - indiv_data[ch][use_ch][0])

    fig = plt.figure(figsize=(10, 6))
    ax  = fig.add_subplot(111)
    for i, ch in enumerate(indiv_data):
        corrected_data = ch[:, 0] + np.full_like(ch[:, 0], offset[i])
        ax.errorbar(
            np.arange(1, 17), corrected_data, yerr = ch[:, 1], 
            fmt = "s", capsize = 0, markeredgecolor = "k", ms = 8, ecolor='black',  color=f'C{i}', markeredgewidth = 0.2, zorder = 3
        )

        mean, error = opg_tool.weighted_mean_and_error(corrected_data, 1/ch[:, 1]**2)
        print(i, format(mean, ".3f"), format(error, ".3f"))

        # 重み付き平均の横線を描画
        ax.hlines(mean, 0.5, 16.5, color=f'C{i}', linestyle='--', label = f"Amp board {i}\n(mean = {mean:.2f})")

        # エラー帯を描画
        ax.fill_between(
            x=[0.5, 16.5],  # x範囲をデータ範囲に拡張
            y1=mean - error,             # エラー範囲下限
            y2=mean + error,             # エラー範囲上限
            color=f'C{i}',
            alpha=0.2
        )
        
    ax.set_xticks(np.arange(1, 16, 2))
    ax.set_xlabel("MPPC ch")
    ax.set_ylabel("One Photon Gain [arb. unit]")
    ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(1.0, 1))
    plt.subplots_adjust(left = 0.13, right = 0.72, top = 0.98, bottom = 0.15)
    plt.savefig(os.path.join(script_dir, f"../results/img/bac_opg_{HV:.0f}.png"), dpi=600, transparent=True)
    plt.show()