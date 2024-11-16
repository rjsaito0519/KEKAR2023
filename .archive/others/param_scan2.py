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
    [30.     , 40.     , 50.     , 60.     , 70.     ],
    [ 0.98851,  0.96559,  0.96559,  0.83765,  0.73805],
    [ 0.98551,  0.9604 ,  0.88404,  0.79037,  0.68116],
    [ 0.98426,  0.93764,  0.85829,  0.72737,  0.60412]
])
SAC_ineff = np.array([
    [30, 40, 50, 60, 70],
    [0.01625, 0.00043, 0.00004, 0, 0],
    [0.01645, 0.00046, 0, 0.00000, 0],
    [0.01646, 0.00046, 0.00000, 0, 0.00000]
])

BAC_eff = np.array([
    [50, 75, 100, 125, 150],
    [0.99513, 0.995, 0.99457, 0.68273, 0.38706],
    [0.998, 0.946, 0.62855, 0.268, 0.13228],
    [0.974, 0.61257, 0.193, 0.09882, 0.06638]
])
BAC_ineff = np.array([
    [50, 75, 100, 125, 150],
    [0.00350, 0, 0.00000, 0.00000, 0],
    [0.00023, 0.00000, 0.00000, 0.00000, 0.00000],
    [0.00000, 0.00000, 0.00000, 0.00000, 0.00000]
])


KVC_eff = np.array([
    [30, 75, 100, 125, 150],
    [1., 1, 1, 1, 0.99995],
    [1., 1, 0.99994, 0.99992, 1],
    [1., 1, 0.99981, 0.999850, 0.99995]
])
KVC_ineff = np.array([
    [30, 75, 100, 125, 150],
    [0.20762, 0.00004, 0, 0, 0],
    [0.03020, 0, 0, 0.00004, 0],
    [0.00146, 0.00004, 0, 0, 0]
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
ax2.plot( 60, 99.844, "o", ms = 7, color = "C1" )
ax5.plot( 60, 0.00041*100, "o", ms = 7, color = "C1" )
ax3.set_ylim(60, 102)
# plt.savefig("../img/Eff.jpg", dpi=600)
plt.show()