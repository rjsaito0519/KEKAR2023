import numpy as np
import matplotlib.pyplot as plt
import uproot
import os
import lmfit as lf
import lmfit.models as lfm

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
for key in data.keys():
    data[key] = np.array(data[key])
    err[key]  = np.array(err[key])
    ax.errorbar(data[key][:, 0], data[key][:, 1], yerr = err[key], marker = marker[key], ls = "--", label = f"HV = {key} V")

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
for key in data.keys():
    data[key] = np.array(data[key])
    err[key]  = np.array(err[key])
    ax.errorbar(data[key][:, 0], data[key][:, 1], yerr = err[key], marker = marker[key], ls = "--", label = f"HV = {key} V")

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