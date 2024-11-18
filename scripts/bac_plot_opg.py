import os
import numpy as np
import matplotlib.pyplot as plt
import uproot

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


def weighted_mean_and_error(values, weights):
    """
    重み付き平均とそのエラーを計算する関数

    Parameters:
        values (list or np.ndarray): データ点
        weights (list or np.ndarray): 重み

    Returns:
        tuple: (重み付き平均, エラー)
    """
    values = np.array(values)
    weights = np.array(weights)

    # 重み付き平均の計算
    weighted_mean = np.sum(weights * values) / np.sum(weights)

    # 標準誤差の計算
    weighted_error = np.sqrt(1 / np.sum(weights))

    return weighted_mean, weighted_error


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

# -- data summarize -----
def data_summarize(ch):
    data = []
    prev_mppc = -1
    val, err = 0, np.inf
    for i in range(len(tree["run_num"])):
        if tree["run_num"][i] in mppc_map.keys() and tree["ch"][i] == ch:
            if prev_mppc != mppc_map[tree["run_num"][i]]:
                if tree["result_err"][i][3] < err:
                    val = tree["result_val"][i][3]
                    err = tree["result_err"][i][3]
                data.append([ val, err ])
                prev_mppc = mppc_map[tree["run_num"][i]]
                val, err = 0, np.inf
            else:
                if tree["result_err"][i][3] < err:
                    val = tree["result_val"][i][3]
                    err = tree["result_err"][i][3]
    return np.array(data)

ch1 = data_summarize(0)
ch2 = data_summarize(1)
ch3 = data_summarize(2)
ch4 = data_summarize(3)

result = []

fig = plt.figure(figsize=(10, 6))
ax  = fig.add_subplot(111)
for i, ch in enumerate([ch1, ch2, ch3, ch4]):
    ax.errorbar(
        np.arange(1, 17), ch[:, 0], yerr = ch[:, 1], 
        fmt = "s", capsize = 0, markeredgecolor = "k", ms = 8, ecolor='black',  color=f'C{i}', markeredgewidth = 0.2, zorder = 3
    )

    mean, error = weighted_mean_and_error(ch[:, 0], 1/ch[:, 1]**2)

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
    
    result.append([ mean, error ])

ax.set_xticks(np.arange(1, 16, 2))
ax.set_xlabel("MPPC ch")
ax.set_ylabel("One Photon Gain [arb. unit]")
ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(1.0, 1))
plt.subplots_adjust(left = 0.13, right = 0.72, top = 0.98, bottom = 0.15)
plt.savefig(os.path.join(script_dir, "../results/img/bac_opg.png"), dpi=600, transparent=True)
plt.show()

result = np.array(result)


HV_map = {
    238: 57.0,
    239: 56.0,
    240: 58.0
}

a = [ch1, ch2, ch3, ch4]

for ch in range(4):
    for i in range(len(tree["run_num"])):
        if tree["run_num"][i] in HV_map.keys() and tree["ch"][i] == ch:
            print( HV_map[tree["run_num"][i]], format(tree["result_val"][i][3], ".3f"), format(tree["result_err"][i][3], ".3f") )
            plt.errorbar(
                HV_map[tree["run_num"][i]], tree["result_val"][i][3], yerr = tree["result_err"][i][3],
                fmt = "s", capsize = 0, markeredgecolor = "k", ms = 8, ecolor='black',  color=f'C{tree["ch"][i]}', markeredgewidth = 0.2, zorder = 3
            )
    # 重み付き平均の横線を描画
    mean, error = a[ch][0]
    plt.axhline(mean, color=f'C{ch}', linestyle='--', label = f"Amp board {ch}")

    # エラー帯を描画
    plt.fill_between(
        x=[56, 58],  # x範囲をデータ範囲に拡張
        y1=mean - error,             # エラー範囲下限
        y2=mean + error,             # エラー範囲上限
        color=f'C{ch}',
        alpha=0.2
    )

plt.show()