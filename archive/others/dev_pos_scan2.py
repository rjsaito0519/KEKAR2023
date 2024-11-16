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

data = np.genfromtxt(
        '../csv_data/efficiency2.csv',
        delimiter=',',
        skip_header=1
    )

SAC = []
BAC = []
KVC = []
for y in y_pos:
    SAC_buf = []
    BAC_buf = []
    KVC_buf = []
    for x in x_pos:
        tmp_data = data[ (data[:, 10] == x) * (data[:, 11] == y) * ( data[:, 3] != 0)  ]
        i_max = np.argmax(tmp_data[:][:, 7])
        SAC_buf.append( tmp_data[i_max][7] )
        BAC_buf.append( tmp_data[i_max][8] )
        KVC_buf.append( tmp_data[i_max][9] )
    SAC.append(SAC_buf)
    BAC.append(BAC_buf)
    KVC.append(KVC_buf)

SAC = np.array(SAC)
BAC = np.array(BAC)
KVC = np.array(KVC)


# for y_index in range(len(y_pos)):
#     fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, figsize=(10,9))

#     text_x = -75
#     text_y = np.max(SAC[y_index]) - ( np.max(SAC[y_index]) - np.min(SAC[y_index]) )/7
#     # ax1.errorbar( x_pos, SAC[y_index], xerr= 9/2, fmt = "o", capsize=5, ms = 12, color = "w", markeredgecolor = "k", ecolor = "k", markeredgewidth = 0.5, zorder = 2)
#     ax1.plot( x_pos, SAC[y_index], ":o", color = "C0", markersize = 8, zorder = 2)
#     ax1.set_ylabel("SAC")
#     ax_info = ax1.axis()
#     ax1.set_ylim( np.min(np.ravel(SAC))*0.9, np.max(np.ravel(SAC))*1.05 )
#     ax1.vlines(-145/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
#     ax1.vlines(+145/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
#     # ax1.fill_betweenx(np.linspace(0, 1), -145/2, +145/2, fc="r", ec = "none", alpha = 0.1, zorder = 1)

#     text_x = -75
#     text_y = np.max(BAC[y_index]) - ( np.max(BAC[y_index]) - np.min(BAC[y_index]) )/7
#     # ax2.errorbar( x_pos, BAC[y_index], xerr= 9/2, fmt = "o", capsize=5, ms = 12, color = "w", markeredgecolor = "k", ecolor = "k", markeredgewidth = 0.5, zorder = 2)
#     ax2.plot( x_pos, BAC[y_index], ":^", color = "C1", markersize = 8, zorder = 2)
#     ax2.set_ylabel("BAC")
#     ax_info = ax2.axis()
#     ax2.set_ylim( np.min(np.ravel(BAC))*0.9, np.max(np.ravel(BAC))*1.05 )
#     # ax2.set_ylim( ax_info[2], ax_info[3] )
#     ax2.vlines(-115/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
#     ax2.vlines(+115/2, 0, 1, ls = "dashed", color = "r", zorder = 1)
#     # ax2.fill_betweenx(np.linspace(0, 1), -115/2, +115/2, fc="r", ec = "none", alpha = 0.1, zorder = 1)

#     text_x = -75
#     text_y = np.max(KVC[y_index]) - ( np.max(KVC[y_index]) - np.min(KVC[y_index]) )/7
#     # ax3.errorbar( x_pos, KVC[y_index], xerr= 9/2, fmt = "o", capsize=5, ms = 12, color = "w", markeredgecolor = "k", ecolor = "k", markeredgewidth = 0.5, zorder = 2)
#     ax3.plot( x_pos, KVC[y_index], ":s", color = "C2", markersize = 8, zorder = 2)
#     ax3.set_ylabel("KVC")
#     ax_info = ax3.axis()
#     # ax3.set_ylim( ax_info[2], ax_info[3] )
#     ax3.set_ylim( np.min(np.ravel(KVC))*0.9, np.max(np.ravel(KVC))*1.05 )
#     for i in range(5):
#         ax3.vlines(-26*2-13 + 26*i, 0, 1, ls = "dashed", color = "r", zorder = 1)
#     # ax3.fill_betweenx(np.linspace(0, 1), -26*2-13, 26+13, fc="r", ec = "none", alpha = 0.1, zorder = 1)

#     plt.xticks(x_pos)
#     fig.suptitle("y position = {} [mm]".format(y_pos[y_index]))
#     plt.subplots_adjust(hspace=.0)
#     fig.supxlabel("x position [mm]")
#     fig.supylabel("Efficiency")
#     plt.savefig("./img/Eff_{}.jpg".format(y_index), dpi=600)
#     plt.clf()
#     plt.close()
#     # plt.show()

# # files = glob.glob("./img/Eff_*.jpg")
# # files.sort()
# # im_tile = concat_tile([
# #     [ cv2.imread(files[0]), cv2.imread(files[1]), cv2.imread(files[2]), cv2.imread(files[3]) ],
# #     [ cv2.imread(files[4]), cv2.imread(files[5]), cv2.imread(files[6]), np.full_like(cv2.imread(files[0]), 255) ]
# # ])
# # cv2.imwrite("./slice_eff.jpg", im_tile)


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
sns.heatmap(df, cmap="viridis", annot=True, cbar=True, fmt='.3f', annot_kws={'fontsize': 15}, cbar_kws=dict(pad=-0.05, shrink=0.79))
plt.vlines(convert_x( -145/2 ), convert_y( -113.8/2 ), convert_y( 113.8/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( 145/2 ), convert_y( -113.8/2 ), convert_y( 113.8/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( -113.8/2 ), convert_x( -145/2 ), convert_x( 145/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( 113.8/2 ), convert_x( -145/2 ), convert_x( 145/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.text(convert_x( -72 ), convert_y( 70 ), "SAC Efficiency at KEK AR test", zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
plt.savefig("./img/SAC_Eff_corr2.jpg", dpi=600)
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
sns.heatmap(df, cmap="viridis", annot=True, cbar=True, fmt='.3f', annot_kws={'fontsize': 15}, cbar_kws=dict(pad=-0.05, shrink=0.79))
plt.vlines(convert_x( -115/2 ), convert_y( -115/2 ), convert_y( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( 115/2 ), convert_y( -115/2 ), convert_y( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( -115/2 ), convert_x( -115/2 ), convert_x( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.hlines(convert_y( 115/2 ), convert_x( -115/2 ), convert_x( 115/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.text(convert_x( -72 ), convert_y( 70 ), "BAC Efficiency at KEK AR test", zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
# plt.savefig("./img/BAC_Eff_corr2.jpg", dpi=600)
plt.savefig("./img/BAC_Eff_corr2.jpg")
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
sns.heatmap(df, cmap="viridis", annot=True, cbar=True, fmt='.3f', annot_kws={'fontsize': 15}, cbar_kws=dict(pad=-0.05, shrink=0.79))
for i in range(5):
    plt.hlines(convert_y( -26*2-13 + 26*i), convert_x( -120/2 ), convert_x( 120/2 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( -120/2 ), convert_y( -26*2-13 ), convert_y( 26+13 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.vlines(convert_x( 120/2 ), convert_y( -26*2-13 ), convert_y( 26+13 ), color = "red", ls = "dashed", lw = 2, zorder = 1)
plt.text(convert_x( -72 ), convert_y( 72 ), "KVC Efficiency at KEK AR test", zorder = 1)
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
plt.savefig("./img/KVC_Eff_corr2.jpg", dpi=600)
# plt.clf()
# plt.close()
plt.show()

# im_tile = concat_tile([
#     [ cv2.imread( "./img/SAC_Eff.jpg" ), cv2.imread( "./img/BAC_Eff.jpg" ) ],
#     [ cv2.imread( "./img/KVC_Eff.jpg" ), np.full_like(cv2.imread( "./img/KVC_Eff.jpg" ), 255) ],
# ])
# cv2.imwrite("./eff_map.jpg", im_tile)