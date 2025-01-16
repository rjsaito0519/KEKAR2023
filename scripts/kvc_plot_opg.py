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
root_file_path = os.path.join(script_dir, "../results/root/kvc_opg.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

mppc_map = {
    52: 1,
    56: 1,
    59: 1,
    62: 1,
    65: 2,
    66: 2,
    71: 3,
    72: 3,
    74: 3,
    77: 4,
    78: 4,
    79: 4,
    82: 5,
    83: 5,
    86: 6,
    87: 6,
    89: 7,
    90: 7,
    94: 8,
    95: 8,
    101: 9,
    102: 9,
    103: 9,
    105: 10,
    106: 10,
    107: 10,
    108: 10,
    110: 11,
    111: 11,
    114: 12,
    115: 12,
    116: 12,
    119: 13,
    120: 13,
    124: 14,
    125: 14,
    126: 14,
    128: 15,
    129: 15,
    131: 16,
    132: 16,
    133: 16
}

HV_map = {
    58 : 402,
    57 : 404,
    56 : 405
}

ch1 = opg_tool.data_summarize(tree, 0, mppc_map)
ch2 = opg_tool.data_summarize(tree, 1, mppc_map)
ch3 = opg_tool.data_summarize(tree, 2, mppc_map)
ch4 = opg_tool.data_summarize(tree, 3, mppc_map)

indiv_data = [ch1, ch2, ch3, ch4]
use_ch = 5

for HV in [56, 57, 58]:
    scale = []
    for ch in range(4):
        for i in range(len(tree["run_num"])):
            if tree["run_num"][i] == HV_map[HV] and tree["ch"][i] == ch:
                scale.append(tree["result_val"][i][3]/indiv_data[ch][use_ch][0])

    # 凡例ラベルを保存するリスト
    legend_labels = []

    fig = plt.figure(figsize=(12, 6))
    ax  = fig.add_subplot(111)
    for i, ch in enumerate(indiv_data):
        corrected_data = scale[i] * ch[:, 0]
        ax.errorbar(
            np.arange(1, 17), corrected_data, yerr = ch[:, 1], 
            fmt = "s", capsize = 0, markeredgecolor = "k", ms = 8, ecolor='black',  color=f'C{i}', markeredgewidth = 0.2, zorder = 3
        )

        # label_str = f"Amp board {i}"
        label_str = ""
        chunk_size = 4
        for j in range(0, len(corrected_data), chunk_size):
            tmp_val = corrected_data[j:j+chunk_size]
            tmp_err = ch[:, 1][j:j+chunk_size]
            
            mean, error = opg_tool.weighted_mean_and_error(tmp_val, 1/tmp_err**2)
            if j == 0:
                label_str += f"seg.{j//chunk_size+1}: mean = {mean:.2f}"
            else:
                label_str += f"\nseg.{j//chunk_size+1}: mean = {mean:.2f}"
            print(i, j, format(mean, ".3f"), format(error, ".3f"))

            ax.hlines(mean, j+0.5, j+4.5, color=f'C{i}', linestyle='--')

            # エラー帯を描画
            ax.fill_between(
                x=[j+0.5, j+4.5],  # x範囲をデータ範囲に拡張
                y1=mean - error,             # エラー範囲下限
                y2=mean + error,             # エラー範囲上限
                color=f'C{i}',
                alpha=0.2
            )

        legend_labels.append(label_str)
        
    ax.set_xticks(np.arange(1, 16, 2))
    ax.set_xlabel("MPPC number")
    ax.set_ylabel("One Photon Gain [arb. unit]")
    for i, label in enumerate(legend_labels):
        ax.plot([], [], color=f'C{i}', linestyle='--', label=label)  # ダミー線を作成して凡例に追加
    
    ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(1.0, 1.03))
    plt.subplots_adjust(left = 0.13, right = 0.72, top = 0.98, bottom = 0.15)
    plt.savefig(os.path.join(script_dir, f"../results/img/kvc_opg_{HV:.0f}.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
    plt.show()