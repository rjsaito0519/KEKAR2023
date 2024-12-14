import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import uproot
import lmfit as lf
import lmfit.models as lfm
from scipy.stats import norm
from scipy.integrate import quad

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
root_file_path = os.path.join(script_dir, "../results/root/explain.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

run_num = 348
info = {
    "T1": [[130260, 130687], [439.604, 608.493]],
    "T2": [[129139, 131213], [371.753, 516.398]],
    "T3": [[127561, 130330], [370.283, 590.894]],
    "T4": [[126430, 129430], [402.628, 630.286]],
    "BACSUM": [[124856, 135130]]
}


def get_hist_data(file, key):
    hist_data = file[key].to_numpy()
    bin_values = hist_data[0]
    bin_edges = hist_data[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, bin_edges, bin_values

def trigger_counter():

    for key in ["T1", "T2", "T3", "T4"]:
        fig = plt.figure(figsize=(12, 8))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        # +-----+
        # | TDC |
        # +-----+
        # -- histogram -----
        center, edge, value = get_hist_data(file, f"{key}t_{run_num}")
        ax1.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="k", zorder = 2)
        ax1.xaxis.set_major_formatter(ptick.EngFormatter())
        ax1.yaxis.set_major_formatter(ptick.EngFormatter())

        # -- gauss fit result -----
        x = tree[f"gauss_{key}t_{run_num}_x"][0]
        y = tree[f"gauss_{key}t_{run_num}_y"][0]
        ax1.plot(x, y, lw = 2, color = "C1", zorder = 3)
        y_min, y_max = ax1.get_ylim()
        ax1.set_ylim(y_min, y_max*1.1)

        # -- draw gate -----
        mu_t = (info[key][0][0]+info[key][0][1])/2
        sigma_t = (info[key][0][1]-info[key][0][0])/10
        ax1.set_xlim(mu_t-6*sigma_t,mu_t+6*sigma_t)
        ax1.fill_betweenx([y_min, y_max*1.1], mu_t-5*sigma_t,mu_t+5*sigma_t, color='C0', alpha=0.1, zorder = 1)
        ax1.annotate("",
            xy=(mu_t-5*sigma_t, y_max),                  # 矢印の先端 (終点)
            xytext=(mu_t+5*sigma_t, y_max),            # 矢印の始点
            arrowprops=dict(
                arrowstyle="<->",     # 矢印のスタイル
                lw=1.5,             # 矢印の線の太さ
                color="k"        # 矢印の色
            ),
        )
        ax1.text(mu_t, y_max*1.01, r"$\pm 5\sigma$", va = "bottom", ha = "center", fontsize = 24)

        # +-----+
        # | ADC |
        # +-----+
        # -- histogram -----
        center, edge, value = get_hist_data(file, f"{key}a_{run_num}")
        ax2.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="k", zorder = 2)
        ax2.xaxis.set_major_formatter(ptick.EngFormatter())
        ax2.yaxis.set_major_formatter(ptick.EngFormatter())

        # -- gauss fit result -----
        x = tree[f"gauss_{key}a_{run_num}_x"][0]
        y = tree[f"gauss_{key}a_{run_num}_y"][0]
        ax2.plot(x, y, lw = 2, color = "C1", zorder = 3)
        y_min, y_max = ax2.get_ylim()
        ax2.set_ylim(y_min, y_max*1.1)

        # -- draw gate -----
        mu_a = (4*info[key][1][0]+3*info[key][1][1])/7
        sigma_a = (info[key][1][1]-info[key][1][0])/7
        ax2.set_xlim(mu_a-4*sigma_a, mu_a+13*sigma_a)
        ax2.fill_betweenx([y_min, y_max*1.1], mu_a-3*sigma_a, mu_a+4*sigma_a, color='C0', alpha=0.1, zorder = 1)
        ax2.annotate("",
            xy=(mu_a-3*sigma_a, y_max),                  # 矢印の先端 (終点)
            xytext=(mu_a, y_max),            # 矢印の始点
            arrowprops=dict(
                arrowstyle="->",     # 矢印のスタイル
                lw=1.5,             # 矢印の線の太さ
                color="k"        # 矢印の色
            ),
        )
        ax2.text(mu_a-1.5*sigma_a, y_max*1.01, r"$-3\sigma$", va = "bottom", ha = "center", fontsize = 24)

        ax2.annotate("",
            xy=(mu_a, y_max),                  # 矢印の先端 (終点)
            xytext=(mu_a+4*sigma_a, y_max),            # 矢印の始点
            arrowprops=dict(
                arrowstyle="<-",     # 矢印のスタイル
                lw=1.5,             # 矢印の線の太さ
                color="k"        # 矢印の色
            ),
        )
        ax2.text(mu_a+2*sigma_a, y_max*1.01, r"$+4\sigma$", va = "bottom", ha = "center", fontsize = 24)

        if key == "T3":
            # -- for multi peak -----
            ax2.annotate("",
                xy=(650, 45),                  # 矢印の先端 (終点)
                xytext=(650, 70),            # 矢印の始点
                arrowprops=dict(
                    arrowstyle="->",     # 矢印のスタイル
                    lw=2,             # 矢印の線の太さ
                    color="C3"        # 矢印の色
                ),
            )
            ax2.text(650, 70, r"$2e^-$hit", va = "bottom", ha = "center", fontsize = 24, color="C3")

            ax2.annotate("",
                xy=(790, 25),                  # 矢印の先端 (終点)
                xytext=(790, 50),            # 矢印の始点
                arrowprops=dict(
                    arrowstyle="->",     # 矢印のスタイル
                    lw=2,             # 矢印の線の太さ
                    color="C3"        # 矢印の色
                ),
            )
            ax2.text(790, 50, r"$3e^-$hit", va = "bottom", ha = "center", fontsize = 24, color="C3")
        elif key == "T4":
            ax2.annotate("",
                xy=(735, 35),                  # 矢印の先端 (終点)
                xytext=(735, 60),            # 矢印の始点
                arrowprops=dict(
                    arrowstyle="->",     # 矢印のスタイル
                    lw=2,             # 矢印の線の太さ
                    color="C3"        # 矢印の色
                ),
            )
            ax2.text(735, 60, r"$2e^-$hit", va = "bottom", ha = "center", fontsize = 24, color="C3")

        ax1.set_xlabel("TDC [arb. unit]")
        ax1.set_title(f"{key} TDC")
        ax2.set_xlabel("ADC [arb. unit]")
        ax2.set_title(f"{key} ADC")
        plt.subplots_adjust(left = 0.1, right = 0.95, top = 0.9, bottom = 0.12)
        img_save_path = os.path.join(script_dir, f"../results/img/explain/{key}.pdf")
        os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
        plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)
        plt.show()

def BACSUMt():
    key = "BACSUM"
    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111)

    # +-----+
    # | TDC |
    # +-----+
    # -- histogram -----
    center, edge, value = get_hist_data(file, f"{key}t_{run_num}")
    ax.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="k", zorder = 2)
    ax.xaxis.set_major_formatter(ptick.EngFormatter())
    ax.yaxis.set_major_formatter(ptick.EngFormatter())

    # -- gauss fit result -----
    x = tree[f"gauss_{key}t_{run_num}_x"][0]
    y = tree[f"gauss_{key}t_{run_num}_y"][0]
    ax.plot(x, y, lw = 2, color = "C1", zorder = 3)
    y_min, y_max = ax.get_ylim()
    ax.set_ylim(y_min, y_max*1.1)

    # -- draw gate -----
    mu_t = (info[key][0][0]+info[key][0][1])/2
    sigma_t = (info[key][0][1]-info[key][0][0])/10
    ax.set_xlim(mu_t-6*sigma_t,mu_t+6*sigma_t)
    ax.fill_betweenx([y_min, y_max*1.1], mu_t-5*sigma_t,mu_t+5*sigma_t, color='C0', alpha=0.1, zorder = 1)
    ax.annotate("",
        xy=(mu_t-5*sigma_t, y_max),                  # 矢印の先端 (終点)
        xytext=(mu_t+5*sigma_t, y_max),            # 矢印の始点
        arrowprops=dict(
            arrowstyle="<->",     # 矢印のスタイル
            lw=1.5,             # 矢印の線の太さ
            color="k"        # 矢印の色
        ),
    )
    ax.text(mu_t, y_max*1.01, r"$\pm 5\sigma$", va = "bottom", ha = "center", fontsize = 24)
    ax.set_xlabel("TDC [arb. unit]")
    ax.set_title(f"{key} TDC")
    plt.subplots_adjust(left = 0.1, right = 0.95, top = 0.9, bottom = 0.12)
    img_save_path = os.path.join(script_dir, f"../results/img/explain/{key}.pdf")
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)
    plt.show()

def BAC_onsumnpe():
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    # +-----+
    # | TDC |
    # +-----+
    # -- histogram -----
    center, edge, value = get_hist_data(file, f"BAConsumnpe_{run_num}_trig")
    ax1.hist(center, bins=edge, weights=value, lw = 1.5, histtype='step', color="k", zorder = 2, label = "Raw")

    # -- histogram -----
    _, _, value_shower = get_hist_data(file, f"BAConsumnpe_shower_{run_num}")
    ax1.hist(center, bins=edge, weights=value_shower, lw = 1.5, histtype='step', color="C2", zorder = 2, label = r"Trigger ADC > +4$\sigma$")
    ax1.xaxis.set_major_formatter(ptick.EngFormatter())
    ax1.yaxis.set_major_formatter(ptick.EngFormatter())
    ax1.set_xlim(0, 280)
    y_min, y_max = ax1.get_ylim()
    ax1.set_xticks([0, 50, 100, 150, 200, 250])
    ax1.set_xticklabels(["", "", "", "", "", ""])
    ax1.legend(fontsize=24)

    # -- explain saturate -----
    ax1.annotate("",
        xy=(260, 35),                  # 矢印の先端 (終点)
        xytext=(260, 100),            # 矢印の始点
        arrowprops=dict(
            arrowstyle="->",     # 矢印のスタイル
            lw=1.5,             # 矢印の線の太さ
            color="C3"        # 矢印の色
        ),
    )
    ax1.text(255, 100, "saturated\nevents", va = "bottom", ha = "center", fontsize = 24, color="k")

    # -- subtracted ----
    ax2.hist(center, bins=edge, weights=value-value_shower, lw = 1.5, histtype='step', color="C0", zorder = 2, label="Subtracted dist.")
    ax2.xaxis.set_major_formatter(ptick.EngFormatter())
    ax2.yaxis.set_major_formatter(ptick.EngFormatter())

    # -- for fitting -----
    mask = (60 < center) * (center <120)
    for _ in range(3):
        x = center[mask]
        y = (value-value_shower)[mask]
        model = lfm.GaussianModel()
        params = model.guess(x = x, data = y)
        result = model.fit(x = x, data = y, params=params, method='leastsq')

        mask = (result.result.params["center"].value - 1.5*result.result.params["sigma"].value < center) * (center < result.result.params["center"].value + 1.5*result.result.params["sigma"].value)

    x = center[mask]
    y = (value-value_shower)[mask]
    model = lfm.GaussianModel()
    params = model.guess(x = x, data = y)
    result = model.fit(x = x, data = y, params=params, method='leastsq')
    print(result.fit_report())
    fit_x = np.linspace(np.min(x), np.max(x), 100)
    fit_y = result.eval_components(x=fit_x)["gaussian"]
    ax2.plot(fit_x, fit_y, color  = "C1", lw = 2.5, label = "Gaussian fit\n" + r"$N_{\rm p. e.}^{\rm mean}$" + " = {:.1f}".format(result.result.params["center"].value))
    ax2.axvline(result.result.params["center"].value, ls = "dashed", color = "gray", zorder = 0)


    # def gaussian(x, mu, sigma):
    #     return (1 / (np.sqrt(2 * np.pi * sigma**2))) * np.exp(-0.5 * ((x - mu) / sigma)**2)

    # probability, error = quad(gaussian, 20, np.inf, args=(result.result.params["center"].value, result.result.params["sigma"].value))
    # print(probability)

    ax2.legend(fontsize = 24)
    ax2.set_xlim(0, 280)
    ax2.set_ylim(y_min, y_max)

    ax1.set_title(r"BAC online-sum $N_{\rm p. e.}$")
    ax2.set_xlabel(r"$N_{\rm p. e.}$")
    plt.subplots_adjust(left = 0.1, right = 0.98, top = 0.9, bottom = 0.12, hspace=0.01)
    img_save_path = os.path.join(script_dir, f"../results/img/explain/NPE_fit.pdf")
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='pdf', bbox_inches='tight', dpi=600, transparent=True)
    plt.show()


def linear_corr():

    values, edges_x, edges_y = file["on_off_correlation"].to_numpy()  # ヒストグラムのデータを取得

    # 値が0またはNaNの部分をマスクする
    masked_values = np.ma.masked_where(values == 0, values)

    plt.figure(figsize=(10, 8))
    mesh = plt.pcolormesh(
        edges_x, edges_y, masked_values.T,  # 転置が必要
        shading="flat", cmap='viridis'
    )
    plt.colorbar(mesh)

    # -- linear func -----
    x = tree["linear_x"][0]
    y = tree["linear_y"][0]
    plt.plot(x, y, color = "C3", lw = 1.5, zorder = 3)

    plt.xlim(0, 3000)
    plt.ylim(0, 220)
    plt.xlabel("Online-Sum ADC [arb. unit]")
    plt.ylabel(r"Offline-Sum $N_{\rm p. e.}$")

    plt.subplots_adjust(left = 0.13, right = 1.04, top = 0.98, bottom = 0.12)
    img_save_path = os.path.join(script_dir, f"../results/img/explain/linear_fit.png")
    os.makedirs(os.path.dirname(img_save_path), exist_ok=True)
    plt.savefig(img_save_path, format='png', bbox_inches='tight', dpi=600, transparent=True)
    plt.show()


if __name__ == '__main__':
    # BACSUMt()
    BAC_onsumnpe()
    # linear_corr()