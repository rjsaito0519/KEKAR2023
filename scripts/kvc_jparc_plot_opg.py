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
root_file_path1 = os.path.join(script_dir, "../results/root/kvc_opg.root")
root_file_path2 = os.path.join(script_dir, "../results/root/kvc_jparc_opg.root")

file1 = uproot.open(root_file_path1)
tree1 = file1["tree"].arrays(library="np")

file2 = uproot.open(root_file_path2)
tree2 = file2["tree"].arrays(library="np")

# 他のものと定義が違うので注意
mppc_map = {
    1: [23],
    2: [26],
    3: [30],
    4: [35],
    5: [24],
    6: [28],
    7: [31],
    8: [36],
    9: [45, 47],
    10: [51, 52],
    11: [55],
    12: [60],
    13: [48],
    14: [52],
    15: [57],
    16: [61],
}


HV_map = {
    58 : 251,
    57 : 249,
    56 : 250
}


ch1 = opg_tool.data_summarize_for_jparc(tree2, mppc_map, True)
ch2 = opg_tool.data_summarize_for_jparc(tree2, mppc_map, False)

indiv_data = [ch1, ch2]
use_ch = 4

for HV in [56, 57, 58]:
    scale = np.nan
    for i in range(len(tree1["run_num"])):
        if tree1["run_num"][i] == HV_map[HV]:
            scale = tree1["result_val"][i][3]/indiv_data[0][use_ch][0]

    # 凡例ラベルを保存するリスト
    legend_labels = []

    fig = plt.figure(figsize=(12, 6))
    ax  = fig.add_subplot(111)
    for i, ch in enumerate(indiv_data):
        corrected_data = scale * ch[:, 0]
        
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
    plt.savefig(os.path.join(script_dir, f"../results/img/kvc_jparc_opg_{HV:.0f}.png"), dpi=600, transparent=True)
    plt.show()