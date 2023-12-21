import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
import matplotlib as mpl
import matplotlib
import seaborn as sns

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


x_pos = np.array([ -64, -48, -32, -16, 0, 16, 32, 48, 64 ])
y_pos = np.array([ 54, 36, 18, 0, -18, -36, -54 ])

SAC = np.array([
    [0.94185, 0.9428, 0.98273, 0.9692, 0.96363, 0.97636, 0.98441, 0.98608, 0.97011],
    [0.98722, 0.99244, 0.99552, 0.99623, 0.99732, 0.99577, 0.9946, 0.99264, 0.98171],
    [0.99465, 0.99015, 0.99031, 0.98865, 0.99069, 0.98863, 0.99082, 0.99279, 0.99309],
    [0.99937, 0.99232, 0.98634, 0.98619, 0.98599, 0.98575, 0.98809, 0.99488, 0.99606],
    [0.99354, 0.98575, 0.98797, 0.9878, 0.98862, 0.98999, 0.99063, 0.99378, 0.99249],
    [0.9839, 0.99017, 0.99408, 0.9951, 0.9974, 0.99682, 0.99727, 0.99699, 0.99017],
    [0.97734, 0.99079, 0.98915, 0.97116, 0.97987, 0.97755, 0.99446, 0.99423, 0.99087]
])

BAC = np.array([
    [ 0.08217, 0.99872, 0.99976, 0.99976, 0.99985, 0.99969, 0.99964, 0.99933, 0.05572 ],
])

KVC = np.array([
    [ 0.70527, 0.99918, 0.99945, 0.97355, 0.99985, 0.95611, 0.99907, 0.0769, 0.09787 ],
])

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, figsize=(12,8))

for i in range(len(SAC)):
    ax1.plot(x_pos, SAC[i], "--o", color = "C0")
ax1.text( -64, 0.992, "SAC", color = "k", fontsize = 24 )
ax_info = ax1.axis()
ax1.set_ylim( ax_info[2], ax_info[3] )
# SAC = patches.Rectangle(xy=(-145/2, 0), width=145, height=1, fc ="C3", fill=True, alpha = 0.2)
# ax1.add_patch(SAC)

ax2.plot(x_pos, BAC[0], "--^", color = "C1")
ax2.text( -64, 0.81, "BAC", color = "k", fontsize = 24 )
ax_info = ax2.axis()
ax2.set_ylim( ax_info[2], ax_info[3] )
BAC = patches.Rectangle(xy=(-115/2, 0), width=115, height=1, fc ="C3", fill=True, alpha = 0.2)
ax2.add_patch(BAC)

ax3.plot(x_pos, KVC[0], "--s", color = "C2")
ax3.text( -64, 0.91, "KVC", color = "k", fontsize = 24 )
ax_info = ax3.axis()
ax3.set_ylim( ax_info[2], ax_info[3] )
# KVC = patches.Rectangle(xy=(-26*2-13, 0), width=26*4, height=1, fc ="C3", fill=True, alpha = 0.2)
# ax3.add_patch(KVC)
for i in range(4):
    KVC = patches.Rectangle(xy=(-26*2-13 + 26*i, 0), width=26, height=1, ec = "k", fc ="C3", fill=True, alpha = 0.2)
    ax3.add_patch(KVC)

plt.subplots_adjust(hspace=.0)
# plt.subplots_adjust(hspace=.5)
# plt.legend()
# plt.savefig("./img/gas/oxygen.svg", dpi=150, transparent=True)
# plt.savefig("../img/Eff.jpg", dpi=600)
plt.show()

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

for i in range(len(SAC)):
    ax.plot(x_pos, np.full_like(x_pos, y_pos[i]), SAC[i], color = "C0")
plt.show()

plt.rcParams['axes.grid'] = False
# im = plt.imshow(SAC, interpolation='none', origin='lower', extent=(-64-8, 64+8, -54-9, 54+9))
# plt.show()

plt.figure(figsize=(10, 8))

df = pd.DataFrame(SAC, columns=x_pos, index = y_pos)
sns.heatmap(df, cmap="viridis", annot=True, cbar=False, fmt='.3f')
sns.lineplot(x=[0, 1], y=[-10, 50], zorder = 1)
plt.show()