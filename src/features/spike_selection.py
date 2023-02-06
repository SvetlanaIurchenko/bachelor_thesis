import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths
import warnings

warnings.filterwarnings("ignore")


def half_width(x, final_spike_train, block):
    SamplingRate = block.segments[0].analogsignals[0].sampling_rate.item()
    prominences, left_bases, right_bases = peak_prominences(x, final_spike_train)
    offset = np.ones_like(prominences)
    widths_ind, h_eval, left_ips, right_ips = peak_widths(
        x, final_spike_train, rel_height=0.5
    )
    right_ips_times = (
        right_ips / block.segments[0].analogsignals[0].sampling_rate.item()
    )
    widths = widths_ind / SamplingRate
    return widths, right_ips_times


def spike_selection(Sweep_number, SpikeTrainInd, corrected_signal, SamplingRate, block):
    Sp_for_ampl = {k: [] for k in range(Sweep_number)}
    for k in range(Sweep_number):
        for i in range(len(SpikeTrainInd[k])):
            Sp_for_ampl[k].append(int(SpikeTrainInd[k][i]))
    ampls_spks = {k: corrected_signal[Sp_for_ampl[k]] for k in range(Sweep_number)}
    widths = {}  # 5ms <= tau_decay <= 40 mc
    right_ips_times = {}
    for k in range(Sweep_number):
        widths[k], right_ips_times[k] = half_width(
            corrected_signal * (-1), SpikeTrainInd[k], block
        )
    tau_decay = {k: np.zeros(len(SpikeTrainInd[k])) for k in range(Sweep_number)}
    for k in range(Sweep_number):
        if len(SpikeTrainInd[k]) != 0:
            for l in range(len(tau_decay[k])):
                x = []
                y = []
                Amax = ampls_spks[k][l]
                for i in range(
                    SpikeTrainInd[k][l], int(right_ips_times[k][l] * SamplingRate) + 2
                ):
                    if (
                        corrected_signal[i] / Amax < 1
                        and corrected_signal[i] / Amax > 0
                    ):
                        y.append(corrected_signal[i] / Amax)
                        x.append(i / SamplingRate)
                    else:
                        continue

                x, y = np.array(x), np.array(y)
                if len(x) > 0:
                    coef, b = np.polyfit(
                        x, np.log(np.array(y)), 1, w=np.sqrt(np.array(y))
                    )
                    tau_dec = (-1 / coef + x[0] / b) / 2
                    leng = len(x)
                    while tau_dec < 0 and leng > 1:
                        leng = int(leng / 1.5)
                        x, y = x[:leng], y[:leng]
                        coef, b = np.polyfit(
                            x, np.log(np.array(y)), 1, w=np.sqrt(np.array(y))
                        )
                        tau_dec = (-1 / coef + x[0] / b) / 2
                    tau_decay[k][l] = tau_dec
                else:
                    tau_decay[k][l] = None
        else:
            tau_decay[k] = []
        tau_decay[k] = tau_decay[k] * 1000
    NewSpikeTrainInd = {}
    for k in tau_decay.keys():
        NewSpikeTrainInd[k] = []
        for i in range(len(tau_decay[k])):
            if 5 <= tau_decay[k][i] <= 40:
                NewSpikeTrainInd[k].append(SpikeTrainInd[k][i])
            else:
                continue
    return NewSpikeTrainInd
