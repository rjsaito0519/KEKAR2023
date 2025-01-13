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

def erf_func(x, amplitude, mean, sigma):
    return amplitude * erf((x - mean) / sigma) + 0.5


file = uproot.open(os.path.join(script_dir, "../results/root/bac_hv_threshold_scan.root"))
tree = file["tree"].arrays(library="np")

threshold_list = [50, 60, 75, 100, 125, 150]

marker = {
    56: "o",
    57: "s",
    58: "^"
}
color = {
    56: "C0",
    57: "C1",
    58: "C2"
}


# +-------------+
# | for explain |
# +-------------+
explain_file = uproot.open(os.path.join(script_dir, "../results/root/explain_bac_hv_threshold_scan.root"))
explain_tree = explain_file["tree"].arrays(library="np")

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# -- hist -----
center, edge, raw_value = get_hist_data(explain_file, f"BAConsumnpe_278_raw")
_, _, trig_value = get_hist_data(explain_file, f"BAConsumnpe_278_trig")

ax1.hist(center, bins=edge, weights=raw_value, lw = 1., histtype='step', color="C0", zorder = 3, label = "raw")
ax1.hist(center, bins=edge, weights=trig_value, lw = 1., histtype='stepfilled', color="C1", alpha = 0.8, zorder = 3, label = "trig.")
ax1.set_xlim(0, 160)
ax1.yaxis.set_major_formatter(ptick.EngFormatter())
ax1.set_xticklabels([])
ax1.axvline(67.6401, color = "k", ls = "dashed", zorder = 5)
ax1.legend()

# -- ratio -----
center, edge, ratio_value = get_hist_data(explain_file, f"onsum_ratio")
ax2.hist(center, bins=edge, weights=ratio_value, lw = 2., histtype='step', color="C0", zorder = 0)
ax2.set_xlim(0, 160)
ax2.set_ylim(-0.05, 1.08)

x = np.linspace(0, 220, 1000)
y = erf_func(x, 0.5, 67.6401, 5.65509)
ax2.plot(x, y, color= "red", lw = 1.5, ls = "dashed")
ax2.axvline(67.6401, color = "k", ls = "dashed", zorder = 5)

ax2.set_xlabel(r"$N_{\rm p.e.}$", fontsize = 30)
ax2.set_ylabel(r"$N_{\rm trig.}/N_{\rm raw}$", fontsize = 30)

fig.suptitle("BAC 3-layer (SUM NPE)", x = 0.15 + (0.98-0.15)/2)

plt.subplots_adjust(left = 0.15, right = 0.98, top = 0.9, bottom = 0.12, hspace=0.02)
plt.savefig(os.path.join(script_dir, f"../results/img/bac/explain_bac_3layer_threshold.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.show()


sys.exit()

# +---------+
# | 3 layer |
# +---------+
info = np.array([
    [58, 50, False],
    [58, 75, False],
    [58, 100, False],
    [58, 125, True],
    [58, 150, True],
    [57, 50, False],
    [57, 75, False],
    [57, 100, True],
    [57, 125, True],
    [57, 150, True],
    [56, 50, True],
    [56, 75, True],
    [56, 100, True],
    [56, 125, True],
    [56, 150, True]
])

data = {56:[], 57:[], 58:[]}
err  = {56:[], 57:[], 58:[]}
for hv, threshold, do_use in info:
    if (do_use):
        indices = np.where((tree["run_num"] <= 400) * (tree["hv"] == hv) * (tree["threshold"] == threshold))[0]
        
        data[hv].append([threshold, tree["onsum_thre_val"][indices[0]][0]])
        err[hv].append(tree["onsum_thre_err"][indices[0]][0])


fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
for i, key in enumerate(data.keys()):
    data[key] = np.array(data[key])
    err[key]  = np.array(err[key])
    ax.errorbar(data[key][:, 0], data[key][:, 1], yerr = err[key], marker = marker[key], ls = "none", color = f"C{i}", markeredgecolor = "k", ecolor= "k", ms = 10, label = f"HV = {key} V", zorder = 5)

    model = lfm.LinearModel()
    params = model.guess(x = data[key][:, 0], data = data[key][:, 1],)
    result = model.fit(x = data[key][:, 0], data = data[key][:, 1], weights = 1/err[key], params=params, method='leastsq')
    # print(result.fit_report())
    fit_x = np.linspace(50, 150, 10)
    fit_y = result.eval_components(x=fit_x)["linear"]
    ax.plot(fit_x, fit_y, color = color[key])
    print(result.eval_components(x=75.0))

ax.set_xlabel("Threshold [mV]")
ax.set_ylabel(r"Threshold [$N_{\rm p. e.}$]")
ax.set_title("BAC 3-layer")
ax.axvline(75, ls = "dashed", color = "red")
ax.legend()

plt.subplots_adjust(left = 0.15, right = 0.98, top = 0.93, bottom = 0.12)
plt.savefig(os.path.join(script_dir, f"../results/img/bac/bac_3layer_threshold.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.show()



# +---------+
# | 2 layer |
# +---------+
info = np.array([
    [58, 50, False],
    [58, 60, False],
    [58, 75, False],
    [58, 100, True],
    [58, 125, True],
    [58, 150, True],
    [57, 50, False],
    [57, 75, True],
    [57, 100, True],
    [57, 125, True],
    [57, 150, True],
    [56, 50, True],
    [56, 75, True],
    [56, 100, True],
    [56, 125, True],
    [56, 150, True]
])

data = {56:[], 57:[], 58:[]}
err  = {56:[], 57:[], 58:[]}
for hv, threshold, do_use in info:
    if (do_use):
        indices = np.where((400 <= tree["run_num"]) * (tree["hv"] == hv) * (tree["threshold"] == threshold))[0]
        
        data[hv].append([threshold, tree["onsum_thre_val"][indices[0]][0]])
        err[hv].append(tree["onsum_thre_err"][indices[0]][0])


fig = plt.figure(figsize=(8, 8))
ax  = fig.add_subplot(111)
for i, key in enumerate(data.keys()):
    data[key] = np.array(data[key])
    err[key]  = np.array(err[key])
    ax.errorbar(data[key][:, 0], data[key][:, 1], yerr = err[key], marker = marker[key], ls = "none", color = f"C{i}", markeredgecolor = "k", ecolor= "k", ms = 10, label = f"HV = {key} V", zorder = 5)

    model = lfm.LinearModel()
    params = model.guess(x = data[key][:, 0], data = data[key][:, 1],)
    result = model.fit(x = data[key][:, 0], data = data[key][:, 1], weights = 1/err[key], params=params, method='leastsq')
    # print(result.fit_report())
    fit_x = np.linspace(50, 150, 10)
    fit_y = result.eval_components(x=fit_x)["linear"]
    ax.plot(fit_x, fit_y, color = color[key])
    print(result.eval_components(x=60.0))

ax.set_xlabel("Threshold [mV]")
ax.set_ylabel(r"Threshold [$N_{\rm p. e.}$]")
ax.set_title("BAC 2-layer")
ax.axvline(60, ls = "dashed", color = "red")
ax.legend()

plt.subplots_adjust(left = 0.15, right = 0.98, top = 0.93, bottom = 0.12)
plt.savefig(os.path.join(script_dir, f"../results/img/bac/bac_2layer_threshold.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.show()