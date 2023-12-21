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
    [0.97734, 0.99079, 0.98915, 0.97116, 0.97987, 0.97755, 0.99446, 0.99423, 0.99087],
    [0.9839, 0.99017, 0.99408, 0.9951, 0.9974, 0.99682, 0.99727, 0.99699, 0.99017],
    [0.99354, 0.98575, 0.98797, 0.9878, 0.98862, 0.98999, 0.99063, 0.99378, 0.99249],
    [0.99937, 0.99232, 0.98634, 0.98619, 0.98599, 0.98575, 0.98809, 0.99488, 0.99606],
    [0.99465, 0.99015, 0.99031, 0.98865, 0.99069, 0.98863, 0.99082, 0.99279, 0.99309],
    [0.98722, 0.99244, 0.99552, 0.99623, 0.99732, 0.99577, 0.9946, 0.99264, 0.98171],
    [0.94185, 0.9428, 0.98273, 0.9692, 0.96363, 0.97636, 0.98441, 0.98608, 0.97011]
])

BAC = np.array([
     [0.0707, 0.76871, 0.82636, 0.83856, 0.86, 0.87061, 0.86419, 0.86763, 0.04629],
    [0.08119, 0.99644, 0.99879, 0.99928, 0.99906, 0.99951, 0.99823, 0.99507, 0.04039],
    [0.07295, 0.99832, 0.9987, 0.99945, 0.99986, 0.99923, 0.99941, 0.99833, 0.04895],
    [0.08217, 0.99872, 0.99976, 0.99976, 0.99985, 0.99969, 0.99964, 0.99933, 0.05572],
    [0.0751, 0.99906, 0.99982, 0.99991, 0.99969, 0.99974, 0.99977, 0.99923, 0.04862],
    [0.06424, 0.99921, 0.99972, 0.99989, 0.99984, 0.99974, 0.9996, 0.99264, 0.04725],
    [0.05193, 0.61017, 0.671, 0.62791, 0.567, 0.57487, 0.59977, 0.54355, 0.04101]
])

KVC = np.array([
    [0.72492, 0.9627, 0.9778, 0.96306, 1, 0.96514, 0.98855, 0.08845, 0.08586],
    [0.75316, 0.99901, 0.99951, 0.96914, 0.99953, 0.96184, 0.99823, 0.07353, 0.08205],
    [0.71960, 0.99916, 0.99951, 0.97394, 0.99954, 0.95751, 0.99896, 0.07784, 0.09146],
    [0.70527, 0.99918, 0.99945, 0.97355, 0.99985, 0.95611, 0.99907, 0.0769, 0.09787],
    [0.64304, 0.99888, 0.99862, 0.98655, 0.99964, 0.95435, 0.99924, 0.080035, 0.09113],
    [0.60476, 0.99937, 0.999, 0.99597, 1.00000, 0.96583, 0.9991, 0.08403, 0.08468],
    [0.55907, 0.97478, 0.96839, 0.96812, 0.96065, 0.92119, 0.94758, 0.0891, 0.07548]
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
plt.savefig("./img/SAC_Eff1.jpg", dpi=600)
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
plt.savefig("./img/BAC_Eff1.jpg", dpi=600)
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
    plt.vlines(convert_x( -26*2-13 + 26*i), convert_y( -120/2 ), convert_y( 120/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( -120/2 ), convert_x( -26*2-13 ), convert_x( 26+13 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( 120/2 ), convert_x( -26*2-13 ), convert_x( 26+13 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.text(convert_x( -72 ), convert_y( 72 ), "KVC Efficiency", zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
plt.savefig("./img/KVC_Eff1.jpg", dpi=600)
# plt.clf()
# plt.close()
plt.show()

# im_tile = concat_tile([
#     [ cv2.imread( "./img/SAC_Eff.jpg" ), cv2.imread( "./img/BAC_Eff.jpg" ) ],
#     [ cv2.imread( "./img/KVC_Eff.jpg" ), np.full_like(cv2.imread( "./img/KVC_Eff.jpg" ), 255) ],
# ])
# cv2.imwrite("./eff_map.jpg", im_tile)