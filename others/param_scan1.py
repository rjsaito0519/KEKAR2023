import numpy as np
import matplotlib.pyplot as plt
import lmfit as lf
import lmfit.models as lfm
import os
import csv
from matplotlib.ticker import ScalarFormatter
from numpy.lib.ufunclike import isposinf
import datetime
import uncertainties

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = True
plt.rcParams['figure.subplot.left'] = 0.10
plt.rcParams['figure.subplot.right'] = 0.98
plt.rcParams['figure.subplot.top'] = 0.98
plt.rcParams['figure.subplot.bottom'] = 0.10
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ

def plot(data, ax, isEff = True):

    ax.plot( data[0], data[1]*100, "o", ms = 7, color = "C1", label = "high" )
    ax.plot( data[0], data[2]*100, "^", ms = 7, color = "C0", label = "middle" )
    ax.plot( data[0], data[3]*100, "s", ms = 7, color = "C2", label = "low" )

    ax.set_xlabel("threshold [mV]")
    if isEff:
        ax.set_ylabel("Efficiency")
        ax.legend(borderaxespad=0.1, handletextpad = 0.5, handlelength=1., fontsize = 18, loc = "lower left")
    else:
        # ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        # ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
        ax.set_ylabel("Ped. Efficiency")
    


SAC_eff = np.array([
    [   30,    40,    50,    60,    70],
    [0.997, 0.985, 0.960, 0.909, 0.818],    # top
    [0.996, 0.983, 0.945, 0.893, 0.817],    # middle
    [0.994, 0.973, 0.911, 0.846, 0.751]  # low
])
SAC_ineff = np.array([
    [   30,    40,    50,    60,    70],
    [0.02120, 0.00081, 0.00004, 0.00004, 0.00008],
    [0.01835, 0.00038, 0.00008, 0.00000, 0.00004],
    [0.01900, 0.00026, 0.00050, 0.00035, 0.00000]
])

BAC_eff = np.array([
    [   50,    75,   100,   125,   150],
    [0.995, 1.000, 0.999, 0.991, 0.894],
    [1.000, 0.999, 0.982, 0.850, 0.537],
    [1.000, 0.977, 0.730, 0.317, 0.159]
])
BAC_ineff = np.array([
    [   50,    75,   100,   125,   150],
    [0.00291, 0.00004, 0.00000, 0.00000, 0.00008],
    [0.00026, 0.00000, 0.00000, 0.00000, 0.00000],
    [0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
])

KVC_eff = np.array([
    [   50,    75,   100,   125,   150],
    [1.000, 1.000, 0.998, 0.986, 0.911],
    [1.000, 0.998, 0.977, 0.877, 0.681],
    [1.000, 0.952, 0.719, 0.422, 0.261]
])
KVC_ineff = np.array([
    [   50,    75,   100,   125,   150],
    [0.13156, 0.00035, 0.00046, 0.00046, 0.00027],
    [0.01143, 0.00038, 0.00035, 0.00015, 0.00038],
    [0.00035, 0.00023, 0.00023, 0.00054, 0.00008]
])


fig = plt.figure( figsize=(12, 8) )
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)


plot(SAC_eff, ax1)
plot(BAC_eff, ax2)
plot(KVC_eff, ax3)
plot(SAC_ineff, ax4, False)
plot(BAC_ineff, ax5, False)
plot(KVC_ineff, ax6, False)


# plt.subplots_adjust(hspace=.5)
# plt.legend()
# plt.savefig("./img/gas/oxygen.svg", dpi=150, transparent=True)
# plt.savefig("../img/Eff.jpg", dpi=600)
plt.show()