import numpy as np
import matplotlib.pyplot as plt
import lmfit as lf
import lmfit.models as lfm

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 22
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = True
plt.rcParams['figure.subplot.left'] = 0.13
plt.rcParams['figure.subplot.right'] = 0.98
plt.rcParams['figure.subplot.top'] = 0.95
plt.rcParams['figure.subplot.bottom'] = 0.15
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ

class CheckGain:
    def __init__(self, data):
        self.data = data
        self.fit_x = np.array([])
        self.fit_y = np.array([])

    def fit(self, isLinear = False):
        if isLinear:
            model = lfm.LinearModel()
        else:
            model = lfm.ExponentialModel()
        
        params = model.guess( x = self.data[:, 0], data = self.data[:, 1] )
        result = model.fit(x = self.data[:, 0], data = self.data[:, 1], params=params, method='leastsq')
        print(result.fit_report())
        self.fit_x = np.linspace(np.min(self.data[:, 0]), np.max(self.data[:, 0]), 10000)
        if isLinear:
            self.fit_y = result.eval_components(x=self.fit_x)["linear"]
        else:
            self.fit_y = result.eval_components(x=self.fit_x)["exponential"]

    def cal_voltage(self, gain):
        diff = np.abs(self.fit_y - np.full_like(self.fit_y, gain))
        index = np.argmin(diff)
        return self.fit_x[index]


T1 = CheckGain(np.array([
    [1501,  676],
    [1551,  720],
    [1600,  800],
    [1650,  880],
    [1700,  992],
    [1751, 1060],
]))

T2 = CheckGain(np.array([
    [1500,  904],
    [1550, 1020],
    [1600, 1140],
    [1650, 1220],
    [1700, 1380],
    [1750, 1480],
]))

T3 = CheckGain(np.array([
    [1500,  920],
    [1550, 1020],
    [1600, 1140],
    [1650, 1280],
    [1700, 1380],
    [1750, 1480],
]))

T4 = CheckGain(np.array([
    [1500,  776],
    [1550,  840],
    [1600,  920],
    [1650,  984],
    [1700, 1060],
    [1750, 1120],
]))

fig, ax = plt.subplots(figsize=(8,8))

for i, T in enumerate([T1, T2, T3, T4]):
    T.fit()
    x = T.fit_x
    y = T.fit_y
    vol = T.cal_voltage(1050)

    ax.plot(T.data[:, 0], T.data[:, 1], "o", color = "C{}".format(i), ms = 8, label = "T{} HV = {:.0f}".format(i+1, vol))
    ax.plot( x, y, "--", color = "C{}".format(i) )

ax.hlines(1050, 1500, 1750, ls = "solid", color = "gray")
ax.set_xlabel("HV [V]")
ax.set_ylabel("Pulse Height [mV]")
ax.legend(fontsize = 18, handletextpad = -0.0)
plt.savefig("../img/gain_trigger_counter.jpg", dpi=600)
plt.show()

