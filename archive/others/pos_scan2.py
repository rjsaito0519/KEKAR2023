import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import seaborn as sns
import glob
import cv2

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = True
plt.rcParams['figure.subplot.left'] = 0.12
plt.rcParams['figure.subplot.right'] = 0.98
plt.rcParams['figure.subplot.top'] = 0.93
plt.rcParams['figure.subplot.bottom'] = 0.10
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ

def convert_x(x):
    return (x + 88)/16

def convert_y(y):
    return (81 - y)/18

def concat_tile(im_list_2d):
    return cv2.vconcat([cv2.hconcat(im_list_h) for im_list_h in im_list_2d])

x_pos = np.array([ -64, -48, -32, -16, 0, 16, 32, 48, 64 ])
y_pos = np.array([ -54, -36, -18, 0, 18, 36, 54 ])

SAC = np.array([
    [0.97119, 0.9366, 0.96524, 0.95589, 0.93285, 0.96555, 0.9635, 0.96477, 0.93791],
    [0.96948, 0.98779, 0.98542, 0.98124, 0.98317, 0.97501, 0.97367, 0.96645, 0.93452],
    [0.97692, 0.97837, 0.9734, 0.96761, 0.96803, 0.96697, 0.96529, 0.96569, 0.9499],
    [0.99345, 0.97787, 0.96842, 0.96501, 0.96269, 0.96394, 0.96456, 0.97677, 0.99795],
    [0.97550, 0.96916, 0.96602, 0.9674, 0.96421, 0.96441, 0.9618, 0.96419, 0.96058],
    [0.96034, 0.97277, 0.97292, 0.977, 0.98179, 0.97647, 0.97173, 0.96453, 0.94071],
    [0.9539, 0.99383, 0.98573, 0.96589, 0.98908, 0.95845, 0.98633, 0.9731, 0.90818]
])

BAC = np.array([
    [0.05017, 0.54832, 0.63763, 0.66872, 0.68306, 0.69379, 0.696389, 0.67837, 0.04595],
    [0.05235, 0.96434, 0.97828, 0.98755, 0.99183, 0.9887, 0.98535, 0.96851, 0.03949],
    [0.05394, 0.98136, 0.98739, 0.99439, 0.99705, 0.99294, 0.9228, 0.98463, 0.04532],
    [0.054431, 0.99604, 0.99747, 0.99797, 0.99743, 0.99779, 0.99806, 0.99517, 0.05084],
    [0.05546, 0.99764, 0.99854, 0.99826, 0.99853, 0.99876, 0.99903, 0.99811, 0.04719],
    [0.05008, 0.99741, 0.99931, 0.9988, 0.99877, 0.99902, 0.99878, 0.99688, 0.04771],
    [0.03986, 0.72262, 0.74776, 0.74342, 0.67257, 0.69598, 0.71501, 0.66393, 0.04235]
])

KVC = np.array([
    [0.4007, 0.99781, 0.999, 0.99932, 0.99922, 0.99948, 0.99896, 0.99812, 0.45598],
    [0.31595, 0.9582, 0.95995, 0.95957, 0.95708, 0.96198, 0.96314, 0.96365, 0.43839],
    [0.29412, 0.99602, 0.99426, 0.99082, 0.98691, 0.98528, 0.97994, 0.97386, 0.45563],
    [0.28165, 0.99916, 0.99989, 1, 1, 0.99991, 1, 0.99928, 0.44732],
    [0.28975, 0.94717, 0.94889, 0.95219, 0.95156, 0.95145, 0.953, 0.95603, 0.43122],
    [0.22784, 0.93762, 0.92919, 0.92274, 0.91134, 0.90217, 0.89921, 0.88805, 0.36945],
    [0.05987, 0.04348, 0.04572, 0.04754, 0.04949, 0.05165, 0.04868, 0.04485, 0.06541]
])


for y_index in range(len(y_pos)):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, figsize=(10,9))

    text_x = -75
    text_y = np.max(SAC[y_index]) - ( np.max(SAC[y_index]) - np.min(SAC[y_index]) )/7
    # ax1.errorbar( x_pos, SAC[y_index], xerr= 9/2, fmt = "o", capsize=5, ms = 12, color = "w", markeredgecolor = "k", ecolor = "k", markeredgewidth = 0.5, zorder = 2)
    ax1.plot( x_pos, SAC[y_index], ":o", color = "C0", markersize = 8, zorder = 2)
    ax1.set_ylabel("SAC")
    ax_info = ax1.axis()
    ax1.set_ylim( np.min(np.ravel(SAC))*0.9, np.max(np.ravel(SAC))*1.05 )
    ax1.vlines(-145/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
    ax1.vlines(+145/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
    # ax1.fill_betweenx(np.linspace(0, 1), -145/2, +145/2, fc="r", ec = "none", alpha = 0.1, zorder = 1)

    text_x = -75
    text_y = np.max(BAC[y_index]) - ( np.max(BAC[y_index]) - np.min(BAC[y_index]) )/7
    # ax2.errorbar( x_pos, BAC[y_index], xerr= 9/2, fmt = "o", capsize=5, ms = 12, color = "w", markeredgecolor = "k", ecolor = "k", markeredgewidth = 0.5, zorder = 2)
    ax2.plot( x_pos, BAC[y_index], ":^", color = "C1", markersize = 8, zorder = 2)
    ax2.set_ylabel("BAC")
    ax_info = ax2.axis()
    ax2.set_ylim( np.min(np.ravel(BAC))*0.9, np.max(np.ravel(BAC))*1.05 )
    # ax2.set_ylim( ax_info[2], ax_info[3] )
    ax2.vlines(-115/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
    ax2.vlines(+115/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
    # ax2.fill_betweenx(np.linspace(0, 1), -115/2, +115/2, fc="r", ec = "none", alpha = 0.1, zorder = 1)

    text_x = -75
    text_y = np.max(KVC[y_index]) - ( np.max(KVC[y_index]) - np.min(KVC[y_index]) )/7
    # ax3.errorbar( x_pos, KVC[y_index], xerr= 9/2, fmt = "o", capsize=5, ms = 12, color = "w", markeredgecolor = "k", ecolor = "k", markeredgewidth = 0.5, zorder = 2)
    ax3.plot( x_pos, KVC[y_index], ":s", color = "C2", markersize = 8, zorder = 2)
    ax3.set_ylabel("KVC")
    ax_info = ax3.axis()
    # ax3.set_ylim( ax_info[2], ax_info[3] )
    ax3.set_ylim( np.min(np.ravel(KVC))*0.9, np.max(np.ravel(KVC))*1.05 )
    for i in range(5):
        ax3.vlines(-26*2-13 + 26*i, 0, 1, ls = "dashed", color = "r", zorder = 1)
    # ax3.fill_betweenx(np.linspace(0, 1), -26*2-13, 26+13, fc="r", ec = "none", alpha = 0.1, zorder = 1)

    plt.xticks(x_pos)
    fig.suptitle("y position = {} [mm]".format(y_pos[y_index]))
    plt.subplots_adjust(hspace=.0)
    fig.supxlabel("x position [mm]")
    fig.supylabel("Efficiency")
    plt.savefig("./img/Eff_{}.jpg".format(y_index), dpi=600)
    plt.clf()
    plt.close()
    # plt.show()

# files = glob.glob("./img/Eff_*.jpg")
# files.sort()
# im_tile = concat_tile([
#     [ cv2.imread(files[0]), cv2.imread(files[1]), cv2.imread(files[2]), cv2.imread(files[3]) ],
#     [ cv2.imread(files[4]), cv2.imread(files[5]), cv2.imread(files[6]), np.full_like(cv2.imread(files[0]), 255) ]
# ])
# cv2.imwrite("./slice_eff.jpg", im_tile)


plt.figure(figsize=(10, 8))
df = pd.DataFrame(SAC, columns=x_pos, index = y_pos)
df[-80] = [np.nan]*7
df[80] = [np.nan]*7
df.loc[72] = np.nan
df.loc[-72] = np.nan
df = df.sort_index(axis=0, ascending=False)
df = df.sort_index(axis=1)
print(df)
plt.grid(which="major", alpha=0.3)
sns.heatmap(df, cmap="viridis", annot=True, cbar=False, fmt='.3f', annot_kws={'fontsize': 15})
plt.vlines(convert_x( -145/2 ), convert_y( -113.8/2 ), convert_y( 113.8/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( 145/2 ), convert_y( -113.8/2 ), convert_y( 113.8/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( -113.8/2 ), convert_x( -145/2 ), convert_x( 145/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( 113.8/2 ), convert_x( -145/2 ), convert_x( 145/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.text(convert_x( -72 ), convert_y( 70 ), "SAC Efficiency", zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
plt.savefig("./img/SAC_Eff2.jpg", dpi=600)
# plt.clf()
# plt.close()
plt.show()

plt.figure(figsize=(10, 8))
df = pd.DataFrame(BAC, columns=x_pos, index = y_pos)
df[-80] = [np.nan]*7
df[80] = [np.nan]*7
df.loc[72] = np.nan
df.loc[-72] = np.nan
df = df.sort_index(axis=0, ascending=False)
df = df.sort_index(axis=1)
print(df)
plt.grid(which="major", alpha=0.3)
sns.heatmap(df, cmap="viridis", annot=True, cbar=False, fmt='.3f', annot_kws={'fontsize': 15})
plt.vlines(convert_x( -115/2 ), convert_y( -115/2 ), convert_y( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( 115/2 ), convert_y( -115/2 ), convert_y( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( -115/2 ), convert_x( -115/2 ), convert_x( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( 115/2 ), convert_x( -115/2 ), convert_x( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.text(convert_x( -72 ), convert_y( 70 ), "BAC Efficiency", zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
plt.savefig("./img/BAC_Eff2.jpg", dpi=600)
# plt.clf()
# plt.close()
plt.show()

plt.figure(figsize=(10, 8))
df = pd.DataFrame(KVC, columns=x_pos, index = y_pos)
df[-80] = [np.nan]*7
df[80] = [np.nan]*7
up = pd.DataFrame([[np.nan for _ in range(len(x_pos+2))]], columns=x_pos, index = [72])
df.loc[72] = np.nan
df.loc[-72] = np.nan
df = df.sort_index(axis=0, ascending=False)
df = df.sort_index(axis=1)
print(df)
plt.grid(which="major", alpha=0.3)
sns.heatmap(df, cmap="viridis", annot=True, cbar=False, fmt='.3f', annot_kws={'fontsize': 15})
for i in range(5):
    plt.hlines(convert_y( -26*2-13 + 26*i), convert_x( -120/2 ), convert_x( 120/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( -120/2 ), convert_y( -26*2-13 ), convert_y( 26+13 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( 120/2 ), convert_y( -26*2-13 ), convert_y( 26+13 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.text(convert_x( -72 ), convert_y( 72 ), "KVC Efficiency", zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
plt.savefig("./img/KVC_Eff2.jpg", dpi=600)
# plt.clf()
# plt.close()
plt.show()

# im_tile = concat_tile([
#     [ cv2.imread( "./img/SAC_Eff.jpg" ), cv2.imread( "./img/BAC_Eff.jpg" ) ],
#     [ cv2.imread( "./img/KVC_Eff.jpg" ), np.full_like(cv2.imread( "./img/KVC_Eff.jpg" ), 255) ],
# ])
# cv2.imwrite("./eff_map.jpg", im_tile)