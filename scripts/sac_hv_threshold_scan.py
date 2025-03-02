import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import uproot

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
root_file_path = os.path.join(script_dir, "../results/root/sac_hv_threshold_scan.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

marker = {
    0: "o",
    1: "s",
    2: "^"
}

label = {
    0: "Low",
    1: "Middle",
    2: "High"
}

# +---------+
# | 3 layer |
# +---------+
fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for HV in [0, 1, 2]:
    data = []
    for i in range(len(tree["hv"])):
        if tree["run_num"][i] < 400 and HV == tree["hv"][i]:
            threshold = tree["threshold"][i]
            n_entry_beam = tree["n_entry_beam"][i]
            n_entry_ped  = tree["n_entry_ped"][i]
            n_hit_beam = tree["n_hit_beam"][i]
            n_hit_ped  = tree["n_hit_ped"][i]
            
            data.append([ threshold, n_hit_beam[0]/n_entry_beam*100, n_hit_ped[0]/n_entry_ped*100 ])
    data = np.array(data)
    ax1.plot(data[:, 0], data[:, 1], marker[HV], ms = 10)    
    ax2.plot(data[:, 0], data[:, 2], marker[HV], ms = 10, label = label[HV])

ax1.set_ylim(0, 105)
ax1.set_ylabel("Efficiency [%]")
ax1.yaxis.set_label_coords(-0.11, 0.5)
ax1.set_xticks([30, 40, 50, 60, 70])
ax1.set_xticklabels(["", "", "", "", ""])

ax2.set_ylim(-0.23, 4.9)
ax2.set_ylabel("Ped. Efficiency [%]")
# ax2.set_ylabel("ノイズ混入率[%]", family='BIZ UDGothic')
ax2.yaxis.set_label_coords(-0.11, 0.5)

ax2.set_xlabel(r"$V_{\rm th}$ [mV]")
ax2.legend(fontsize = 24, handletextpad = 0.5, handlelength=0.5)

ax1.set_title("SAC 3-layer")
plt.subplots_adjust(left = 0.17, right = 0.98, top = 0.93, bottom = 0.12, hspace=0.02)
plt.savefig(os.path.join(script_dir, f"../results/img/sac/sac_3layer_hv_threshold_scan.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.savefig(os.path.join(script_dir, f"../results/img/sac/sac_3layer_hv_threshold_scan.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()



# +---------+
# | 2 layer |
# +---------+
fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for HV in [0, 1, 2]:
    data = []
    for i in range(len(tree["hv"])):
        if tree["run_num"][i] > 400 and HV == tree["hv"][i]:
            threshold = tree["threshold"][i]
            n_entry_beam = tree["n_entry_beam"][i]
            n_entry_ped  = tree["n_entry_ped"][i]
            n_hit_beam = tree["n_hit_beam"][i]
            n_hit_ped  = tree["n_hit_ped"][i]
            
            data.append([ threshold, n_hit_beam[0]/n_entry_beam*100, n_hit_ped[0]/n_entry_ped*100 ])
    data = np.array(data)
    ax1.plot(data[:, 0], data[:, 1], marker[HV], ms = 10)    
    ax2.plot(data[:, 0], data[:, 2], marker[HV], ms = 10, label = label[HV])

ax1.set_ylim(0, 105)
ax1.set_ylabel("Efficiency [%]")
ax1.yaxis.set_label_coords(-0.11, 0.5)
ax1.set_xticks([30, 40, 50, 60, 70])
ax1.set_xticklabels(["", "", "", "", ""])

# ax2.set_ylim(-0.04, 0.84)
ax2.set_ylim(-0.23, 4.9)
ax2.set_ylabel("Ped. Efficiency [%]")
# ax2.set_ylabel("ノイズ混入率[%]", family='BIZ UDGothic')
ax2.yaxis.set_label_coords(-0.11, 0.5)

ax2.set_xlabel(r"$V_{\rm th}$ [mV]")
ax2.legend(fontsize = 24, handletextpad = 0.5, handlelength=0.5)

ax1.set_title("SAC 2-layer")
plt.subplots_adjust(left = 0.17, right = 0.98, top = 0.93, bottom = 0.12, hspace=0.02)
plt.savefig(os.path.join(script_dir, f"../results/img/sac/sac_2layer_hv_threshold_scan.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.savefig(os.path.join(script_dir, f"../results/img/sac/sac_2layer_hv_threshold_scan.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()