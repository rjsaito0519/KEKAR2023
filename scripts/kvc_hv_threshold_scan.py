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

marker = {
    56: "o",
    57: "s",
    58: "^"
}

# +----------+
# | KVC 1 cm |
# +----------+
root_file_path = os.path.join(script_dir, "../results/root/kvc_1cm_hv_threshold_scan.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for HV in [56, 57, 58]:
    data = []
    for i in range(len(tree["hv"])):
        if HV == tree["hv"][i]:
            threshold = tree["threshold"][i]
            n_entry_beam = tree["n_entry_beam"][i]
            n_entry_ped  = tree["n_entry_ped"][i]
            n_hit_beam = tree["n_hit_beam"][i]
            n_hit_ped  = tree["n_hit_ped"][i]
            
            data.append([ threshold, n_hit_beam[2]/n_entry_beam*100, n_hit_ped[2]/n_entry_ped*100 ])
    data = np.array(data)
    ax1.plot(data[:, 0], data[:, 1], marker[HV], ms = 10)    
    ax2.plot(data[:, 0], data[:, 2], marker[HV], ms = 10, label = f"HV = {HV} V")

ax1.set_ylim(0, 105)
ax1.set_ylabel("R [%] (w/ beam)")
ax1.fill_betweenx([0, 105], 72.5, 77.5, color='C3', alpha=0.1, zorder = 0)
ax1.set_xticks([50, 75, 100, 125, 150])
ax1.set_xticklabels(["", "", "", "", ""])

ax2.set_ylim(-0.02, 0.23)
ax2.set_ylabel("R [%] (w/o beam)")
ax2.fill_betweenx([-1, 1], 72.5, 77.5, color='C3', alpha=0.1, zorder = 0)

ax2.set_xlabel(r"$V_{\rm th}$ [mV]")
ax2.legend(fontsize = 24)

ax1.set_title("KVC 1 cm")
plt.subplots_adjust(left = 0.15, right = 0.98, top = 0.93, bottom = 0.12, hspace=0.02)
plt.savefig(os.path.join(script_dir, f"../results/img/kvc/kvc_1cm_hv_threshold_scan.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.show()



# +----------+
# | KVC 1 cm |
# +----------+
root_file_path = os.path.join(script_dir, "../results/root/kvc_2cm_hv_threshold_scan.root")
file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for HV in [56, 57, 58]:
    data = []
    for i in range(len(tree["hv"])):
        if HV == tree["hv"][i]:
            threshold = tree["threshold"][i]
            n_entry_beam = tree["n_entry_beam"][i]
            n_entry_ped  = tree["n_entry_ped"][i]
            n_hit_beam = tree["n_hit_beam"][i]
            n_hit_ped  = tree["n_hit_ped"][i]
            
            data.append([ threshold, n_hit_beam[1]/n_entry_beam*100, n_hit_ped[1]/n_entry_ped*100 ])
    data = np.array(data)
    ax1.plot(data[:, 0], data[:, 1], marker[HV], ms = 10)    
    ax2.plot(data[:, 0], data[:, 2], marker[HV], ms = 10, label = f"HV = {HV} V")

ax1.set_ylim(0, 105)
ax1.set_ylabel("R [%] (w/ beam)")
ax1.fill_betweenx([0, 105], 100-2.5, 100+2.5, color='C3', alpha=0.1, zorder = 0)
ax1.set_xticks([50, 75, 100, 125, 150])
ax1.set_xticklabels(["", "", "", "", ""])

ax2.set_ylim(-0.02, 0.23)
ax2.set_ylabel("R [%] (w/o beam)")
ax2.fill_betweenx([-1, 1], 100-2.5, 100+2.5, color='C3', alpha=0.1, zorder = 0)

ax2.set_xlabel(r"$V_{\rm th}$ [mV]")
ax2.legend(fontsize = 24)

ax1.set_title("KVC 2 cm")
plt.subplots_adjust(left = 0.15, right = 0.98, top = 0.93, bottom = 0.12, hspace=0.02)
plt.savefig(os.path.join(script_dir, f"../results/img/kvc/kvc_2cm_hv_threshold_scan.pdf"), format='pdf', bbox_inches='tight', dpi=600, transparent=True)
plt.show()

