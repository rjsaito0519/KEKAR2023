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

def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values

def gauss(x, amp, mu, sigma):
    return amp/(np.sqrt(2*np.pi)*sigma) * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Constant                  =      8215.06   +/-   110.825
# Mean                      =      64.1356   +/-   0.196569
# Sigma                     =      13.6826   +/-   0.251407        (limited)
# 43.9051, 84.5152

file = uproot.open(os.path.join(script_dir, "../results/root/bac_pos_scan_analysis_test.root"))
tree = file["tree"].arrays(library="np")

fig = plt.figure(figsize=(9, 7))
ax  = fig.add_subplot(111)

center, edge, value = get_hist_data(file, "BAConsumnpe_330_trig")
ax.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="k", zorder = 1, label = r"SUM $N_{\rm p. e.}$")

x = np.linspace(43.9051, 84.5152, 5000)
y = gauss(x, 8215.06, 64.1356, 13.6826)
ax.plot(x, y, color = "C1", ls = "solid", lw = 3, zorder = 3, label = r"$N_{\rm p. e.}^{\rm mean}$" + " = 64.14\n" + r"$\sigma$" + " = 13.7")
ax.axvline(15, ls ="dashed", color = "red")

x = np.linspace(0, 100, 5000)
y = gauss(x, 8215.06, 64.1356*0.82, 13.6826*0.82)
ax.plot(x, y, color = "C0", ls = "dashed", lw = 2, zorder = 2, label = r"scaled gaussian")

ax.legend(fontsize = 24)
ax.set_xlim(0, 149)
ax.yaxis.set_major_formatter(ptick.EngFormatter())
ax.set_xlabel(r"$N_{\rm p.e.}$")

ax.set_title("BAC SUM (3-layer)")

plt.subplots_adjust(left = 0.1, right = 0.98, top = 0.9, bottom = 0.12, hspace=0.01)
plt.savefig(os.path.join(script_dir, f"../results/img/explain/bac_fit_example.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()