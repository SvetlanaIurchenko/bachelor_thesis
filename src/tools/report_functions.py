import pandas as pd
import numpy as np
import math
from functools import partial


def sweep_seria_correspondence(x, colors):
    for k in colors.keys():
        if colors[k][0] <= x <= colors[k][1]:
            return k


def electrical_param_table(signal, signal_features):
    r_s, r_in, c = [], [], []  # composite all electric cell parameters
    for k in range(signal.sweep_number):
        r_s.append(signal_features.r_s_dict[k])
        r_in.append(signal_features.r_in_dict[k])
    #     c.append(capacity[k])
    electric_cell_param = pd.DataFrame(r_s, columns=["R_s, Mohm"])
    electric_cell_param["R_in, Mohm"] = r_in
    # electric_cell_param['C, pF'] = c
    electric_cell_param["sweep_indx"] = [i for i in range(signal.sweep_number)]

    sweep_seria_correspondence_fun = partial(
        sweep_seria_correspondence, colors=signal_features.colors
    )

    electric_cell_param["series"] = electric_cell_param.sweep_indx.apply(
        sweep_seria_correspondence_fun
    )
    electric_cell_param.to_csv(
        signal.report_folder / "electric_cell_param.csv", index=True
    )


def baseline_report(signal, signal_features):
    Ihold = []
    Inoise = []
    for k in range(signal.sweep_number):
        Ihold.append(signal_features.distribution_param_sweep[k][0])
        Inoise.append(signal_features.distribution_param_sweep[k][1])
    baseline_param = pd.DataFrame(Ihold, columns=["Ihold, pA"])
    baseline_param["Inoise, pA"] = Inoise
    baseline_param["sweep_indx"] = [i for i in range(signal.sweep_number)]

    sweep_seria_correspondence_fun = partial(
        sweep_seria_correspondence, colors=signal_features.colors
    )

    baseline_param["series"] = baseline_param.sweep_indx.apply(
        sweep_seria_correspondence_fun
    )
    baseline_param.to_csv(signal.report_folder / "baseline_param.csv", index=True)


def tonic_cur_report(signal, signal_features):
    ton_cur = pd.DataFrame(signal_features.tonic_cur, columns=["Iton, pA"])
    ton_cur.to_csv(signal.report_folder / "tonic_cur.csv", index=True)


def make_result_feature_table(signal, signal_features):
    event_param = pd.DataFrame(
        np.concatenate(list(signal_features.new_spike_train_ind.values())),
        columns=["event_idx"],
    )
    event_param["spike_time"] = np.round(
        np.concatenate(list(signal_features.SpikeTrainTimes.values())), 3
    )
    event_param["amplitudes"] = np.round(
        np.concatenate(list(signal_features.ampls_spks_pos.values())), 3
    )
    event_param["tau_decay"] = np.round(
        np.concatenate(list(signal_features.tau_decay.values())), 3
    )
    event_param["rise_rate_10_90"] = np.round(
        np.concatenate(list(signal_features.rise_rate.values())), 3
    )

    start_sweep = [
        k * signal.sweep_length * signal.sampling_rate
        for k in range(signal.sweep_number)
    ]
    stop_sweep = [
        (k + 1) * signal.sweep_length * signal.sampling_rate
        for k in range(signal.sweep_number)
    ]

    def sweep_numeration(x):
        results_keys = [k for k in range(signal.sweep_number)]
        for i in range(len(start_sweep)):
            if start_sweep[i] <= x <= stop_sweep[i]:
                return results_keys[i]

    event_param["sweep"] = event_param.event_idx.apply(sweep_numeration)

    sweep_seria_correspondence_fun = partial(
        sweep_seria_correspondence, colors=signal_features.colors
    )

    event_param["series"] = event_param.sweep.apply(sweep_seria_correspondence_fun)
    event_param.to_csv(signal.report_folder / "event_param.csv", index=True)

    event_param_ctrl = event_param[event_param["series"] == "CTRL"]
    event_param_gaba1 = event_param[event_param["series"] == "GABA1"]
    event_param_gaba5 = event_param[event_param["series"] == "GABA5"]

    frames = [event_param_ctrl, event_param_gaba1, event_param_gaba5]
    res = pd.concat(frames).reset_index(drop=True)  # without picro

    res.to_csv(signal.report_folder / "event_param_3_series.csv", index=True)
    return res


def make_result_table_spont_freq(signal, signal_features):
    spont_frequencies = pd.DataFrame.from_dict(
        signal_features.new_spont_freq, orient="index"
    )
    spont_frequencies["sweep"] = signal_features.spont_freq.keys()

    sweep_seria_correspondence_fun = partial(
        sweep_seria_correspondence, colors=signal_features.colors
    )

    spont_frequencies["series"] = spont_frequencies.sweep.apply(sweep_seria_correspondence_fun)

    new_spont_freq = [[] for l in range(len(signal_features.colors.keys()))]

    for i, k in zip(
        range(len(signal_features.colors.keys())), signal_features.colors.keys()
    ):
        new_spont_freq[i] = spont_frequencies.loc[
            spont_frequencies["series"] == k
        ].to_numpy()

    new_spont_frequencies = [[] for l in range(len(signal_features.colors.keys()))]
    for l in range(len(new_spont_frequencies)):
        for i in range(len(new_spont_freq[l])):
            for j in range(len(new_spont_freq[l][i])):
                if type(new_spont_freq[l][i][j]) == float:
                    new_spont_frequencies[l].append(new_spont_freq[l][i][j])

    new_spont_frequencies_no_nan = [
        [] for l in range(len(signal_features.colors.keys()))
    ]
    spont_freqs = [[] for l in range(len(signal_features.colors.keys()))]
    for l, k in zip(
        range(len(new_spont_frequencies_no_nan)), signal_features.colors.keys()
    ):
        new_spont_frequencies_no_nan[l] = [
            item for item in new_spont_frequencies[l] if not (math.isnan(item)) == True
        ]
        spont_freqs[l] = pd.DataFrame(
            new_spont_frequencies_no_nan[l], columns=["spont_freq"]
        )
        spont_freqs[l]["series"] = [k] * len(spont_freqs[l])

    spont_freq = pd.concat([spont_freqs[l] for l in range(len(spont_freqs))]).reset_index(drop=True)
    spont_freq_ = spont_freq[spont_freq['series'] != 'PTX']
    spont_freq_.to_csv(signal.report_folder / "spont_freqs.csv", index=True)
    return spont_freq_
