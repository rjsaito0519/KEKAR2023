import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import uproot
import os
import lmfit as lf
import lmfit.models as lfm
import sys
from scipy.special import erf

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

def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values

def gauss(x, amp, mu, sigma):
    return amp/(np.sqrt(2*np.pi)*sigma) * np.exp(-(x - mu)**2 / (2 * sigma**2))

def multi_gauss(x, center, sigma1, sigma2, opg, *amps):
    result = gauss(x, amps[0], center, sigma1)
    for i in range(1, len(amps)):
        mu = center + i * opg
        sigma = np.sqrt(i) * sigma2
        result += gauss(x, amps[i], mu, sigma)
    return result

run_num = 162
ch = 0
n = 7
file = uproot.open(os.path.join(script_dir, "../../root/kekar_run{:0=5}.root".format(run_num)))
tree = file["tree"].arrays(library="np")

fit_file = uproot.open(os.path.join(script_dir, "../results/root/bac_opg.root"))
fit_tree = fit_file["tree"].arrays(library="np")

index = np.where((fit_tree["run_num"] == run_num) * (fit_tree["ch"] == ch))[0]


fig = plt.figure(figsize=(10, 8))
ax  = fig.add_subplot(111)

hist_data = ax.hist( tree["baca"][:, 0], bins = np.linspace(0, 4096, 4097), histtype='step', color = "k", lw = 1.5, zorder = 3)

par = fit_tree["result_val"][index][0]
x = np.linspace(175.0, 273.2, 5000)
y = gauss(x, par[0], par[1], par[2])
ax.plot(x, y, color = "C0", ls = "dashed", zorder = 1)
y_tot = y
for i in range(1, n):
    y = gauss(x, par[i+4], par[1] + i*par[3], np.sqrt(i+1)*par[4])
    ax.plot(x, y, color = "C0", ls = "dashed", zorder = 1)
    y_tot += y
ax.plot(x, y_tot, color = "C1", lw = 2, zorder = 5)

y_min, y_max = ax.get_ylim()
ax.set_ylim(y_min, y_max*1.1)

ax.axvline(par[1], color = "red", ls = "dashed", zorder = 7)
ax.axvline(par[1] + par[3], color = "red", ls = "dashed", zorder = 7)
ax.annotate("",
    xy=(par[1], y_max),                  # 矢印の先端 (終点)
    xytext=(par[1]+par[3], y_max),            # 矢印の始点
    arrowprops=dict(
        arrowstyle="<->",     # 矢印のスタイル
        lw=1.5,             # 矢印の線の太さ
        color="k"        # 矢印の色
    ),
)
ax.text(par[1]+par[3]*1.03, y_max, "one photon gain", va = "center", ha = "left", fontsize = 24)



ax.set_xlim(160, 290)
ax.yaxis.set_major_formatter(ptick.EngFormatter())
ax.set_xlabel("ADC [arb. unit]")

ax.set_title("BAC 3-layer (Ch.1, MPPC No.9)")

plt.subplots_adjust(left = 0.13, right = 0.98, top = 0.93, bottom = 0.12)
plt.savefig(os.path.join(script_dir, f"../results/img/explain/opg_fit_example.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.show()