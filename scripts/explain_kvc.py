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

# Constant                  =      7419.88   +/-   105.645
# Mean                      =      49.2767   +/-   0.160137
# Sigma                     =      10.4676   +/-   0.20308         (limited)
# 33.7597, 64.912

file = uproot.open(os.path.join(script_dir, "../results/root/kvc_thin_pos_scan_analysis_test.root"))
tree = file["tree"].arrays(library="np")

fig = plt.figure(figsize=(10, 8))
ax  = fig.add_subplot(111)

center, edge, value = get_hist_data(file, "KVConsumnpe_304_3_trig")
ax.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="k", zorder = 1)

# x = np.linspace(0, 33.7597, 5000)
# y = gauss(x, 7419.88, 49.2767, 10.4676)
# ax.plot(x, y, color = "C0", ls = "dashed", lw = 2, zorder = 2)

x = np.linspace(33.7597, 64.912, 5000)
y = gauss(x, 7419.88, 49.2767, 10.4676)
ax.plot(x, y, color = "C1", ls = "solid", lw = 2, zorder = 3)

ax.set_xlim(0, 99)
ax.axvline(20, ls = "dashed", color = "red")
ax.yaxis.set_major_formatter(ptick.EngFormatter())
ax.set_xlabel(r"$N_{\rm p.e.}$")

ax.set_title("KVC SUM seg. 3 (1 cm)")

plt.subplots_adjust(left = 0.13, right = 0.98, top = 0.93, bottom = 0.12)
plt.savefig(os.path.join(script_dir, f"../results/img/explain/kvc_fit_example.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.show()