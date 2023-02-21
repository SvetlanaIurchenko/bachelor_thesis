import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from loguru import logger
from src.tools.build_common_statistic_function import add_median_labels

def plot_electric_params(r_dict, signal, signal_features, **kwargs):
    fig, ax = plt.subplots(1, figsize=(19, 7))
    for k, v in r_dict.items():
        x_k = [k * signal.sweep_length + 1 / signal.sampling_rate]
        col = "grey"
        for c_k, c_v in signal_features.colors.items():
            if k == c_v[0]:
                col = c_v[2]
                plt.plot(x_k, v, color=col, linewidth=3, label=c_k)
            if c_v[0] < k < c_v[1]:
                col = c_v[2]
        plt.scatter(x_k, v, color=col, linewidth=3)
    plt.xlabel(kwargs.get("xlabel"), fontdict={"fontsize": 20})
    plt.ylabel(kwargs.get("ylabel"), fontdict={"fontsize": 20})
    plt.ylim(kwargs.get("ylim"))
    plt.legend(fontsize=20)
    if kwargs.get("plotname"):
        plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} time_dynamics.png")
        plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} time_dynamics.pdf")
    plt.close()


def plt_I_hold(signal, signal_features, is_noise=False, **kwargs):
    fig, ax = plt.subplots(1, figsize=(19, 7))
    for k, v in signal_features.distribution_param_sweep.items():
        x_k = [
            k * signal.sweep_length + i / signal.sampling_rate for i in range(len(v))
        ]
        col = "grey"
        for c_k, c_v in signal_features.colors.items():
            if k == c_v[0]:
                col = c_v[2]
                plt.plot(x_k[0], v[int(is_noise)], color=col, linewidth=3, label=c_k)
            if c_v[0] < k < c_v[1]:
                col = c_v[2]
        plt.scatter(x_k[0], v[int(is_noise)], color=col, linewidth=3)
    plt.xlabel("Time, s", fontdict={"fontsize": 20})
    plt.ylabel("I_hold, pA", fontdict={"fontsize": 20})
    plt.legend(fontsize=20)
    plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} time_dynamics.png")
    plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} time_dynamics.png")
    plt.close()


def plt_event_param(feature, signal, signal_features, **kwargs):
    fig, ax = plt.subplots(1, figsize=(19, 7))
    for k, v in feature.items():
        if v is None or len(v) == 0:
            for c_k, c_v in signal_features.colors.items():
                if k == c_v[0]:
                    col = c_v[2]
                    plt.scatter(k * signal.sweep_length, 0, color=col, label=c_k)
                if c_v[0] < k < c_v[1]:
                    col = c_v[2]
                plt.scatter(k * signal.sweep_length, 0, color=col, linewidth=3)
        else:
            x_k = signal_features.SpikeTrainTimes[k]
            col = "grey"
            for c_k, c_v in signal_features.colors.items():
                if k == c_v[0]:
                    col = c_v[2]
                    plt.scatter(x_k, v, color=col, label=c_k)
                if c_v[0] < k < c_v[1]:
                    col = c_v[2]
            plt.scatter(x_k, v, color=col, linewidth=3)

    plt.xlabel(kwargs.get("xlabel"), fontdict={"fontsize": 20})
    plt.ylabel(kwargs.get("ylabel"), fontdict={"fontsize": 20})
    plt.ylim(kwargs.get("ylim"))
    plt.legend(fontsize=20)
    if kwargs.get("plotname"):
        plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} time_dynamics.png")
        plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} time_dynamics.pdf")
    plt.close()


def box_plot(signal, data, y, **kwargs):
    logger.info("Box plotting.")
    plt.figure(figsize=(10, 10))
    ax = sns.boxplot(
        x="series",
        y=y,
        data=data,
        palette=["rosybrown", "lightcoral", "bisque", "paleturquoise"],
        linewidth=1.3,
        showfliers=False,
    )
    ax.set_title(f"{kwargs.get('title')}", fontsize=18)
    ax.set_xlabel("series", fontsize=16)
    ax.set_ylabel(f"{kwargs.get('ylabel')}", fontsize=16)
    add_median_labels(ax)
    plt.xticks(fontsize=15, rotation=0)
    plt.yticks(fontsize=15, rotation=0)
    plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} box plot.png")
    plt.savefig(signal.report_folder / f"{kwargs.get('plotname')} box plot.png")
    plt.close()


def plt_box_plot(parametres, lparametres, res, signal):
    for p, l in zip(parametres, lparametres):
        box_plot(
            signal, res, p, **{"title": p.capitalize(), "ylabel": l, "plotname": p}
        )


def plt_cumulative_distribution(
    column, llabel, spont_freq, res, signal, signal_features
):
    logger.info("Plotting cumulative distributions.")
    for c, l in zip(column, llabel):
        if c == "spont_freq":
            fig, ax = plt.subplots(1, figsize=(10, 7))
            for k in signal_features.colors.keys():
                count, bins_count = np.histogram(
                    np.array(spont_freq[spont_freq["series"] == k][c]), bins=50
                )
                pdf = count / sum(count)
                cdf = np.cumsum(pdf)
                plt.plot(
                    np.log(bins_count[1:]),
                    cdf,
                    color=signal_features.colors[k][2],
                    label=k,
                )
        else:
            fig, ax = plt.subplots(1, figsize=(10, 7))
            for k in signal_features.colors.keys():
                count, bins_count = np.histogram(
                    np.array(res[res["series"] == k][c]), bins=50
                )
                pdf = count / sum(count)
                cdf = np.cumsum(pdf)
                plt.plot(
                    np.log(bins_count[1:]),
                    cdf,
                    color=signal_features.colors[k][2],
                    label=k,
                )

        plt.ylabel("p", fontdict={"fontsize": 20})
        plt.xlabel(f"{l}", fontdict={"fontsize": 20})
        plt.legend(fontsize=20)
        plt.savefig(signal.report_folder / f"{c} cumulative distribution.png")
        plt.savefig(signal.report_folder / f"{c} cumulative distribution.pdf")
        plt.close()
