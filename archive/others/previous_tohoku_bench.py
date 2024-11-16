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

def cal_voltage(gain, a, b):
    x = np.linspace(1500, 2400, 500000)
    y = a*x + b
    diff = np.abs(y - np.full_like(y, gain))
    index = np.argmin(diff)
    return x[index]

data = np.genfromtxt(
    "./csv_data/previous_tohoku_bench.csv",
    skip_header=1,
    delimiter=","
)

# ---  use old ch1 as a reference ----------------------
ref1 = data[0][1] * 2200 + data[0][3]
ref2 = data[0][1] * 2150 + data[0][3]
ref3 = data[0][1] * 2100 + data[0][3]

param1 = []
param2 = []
param3 = []

for tmp_data in data:
    tmp_a = tmp_data[1]
    tmp_b = tmp_data[3]
    for ref, param in zip([ref1, ref2, ref3], [param1, param2, param3]):
        vol = cal_voltage(ref, tmp_a, tmp_b)
        param.append([ int(tmp_data[0]), int(np.round(vol)) ])

param1 = np.array(param1)
param2 = np.array(param2)
param3 = np.array(param3)


print("-"*20)
print(ref1)
print(param1)
print("-"*20)
print(ref2)
print(param2)
print("-"*20)
print(ref3)
print(param3)
print("-"*20)
