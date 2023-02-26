from pybaselines import misc
from scipy import signal as sg
import numpy as np
import pandas as pd
import quantities as pq
import elephant
from collections import defaultdict
from loguru import logger
import typing
from scipy.signal import find_peaks, peak_prominences, peak_widths

from src.data.signal_model import Signal

BASELINE_FREQ_CUTOFF = 0.005
BASELINE_LAMBDAS = [0.7, 0.05, 0.2]
BASELINE_ASYMMETRY = 0.6
BRAZHE_CONSTANT = 1.4826


def bound_counting(sweep_number, signal, offset=0):
    return int((sweep_number * signal.sweep_length + offset) * signal.sampling_rate)


def mad_std(v):
    mad = np.median(np.abs(v - np.median(v)))
    return mad * BRAZHE_CONSTANT


class SignalFeatures:
    def __init__(self, signal, colors):
        self.signal = signal
        self.colors = colors
        logger.info("Baseline calculation started.")
        self.baseline = self.baseline_correction(
            array=np.array(self.signal.raw_signal)[0]
        )

        self.save_baseline_to_csv(signal=signal, baseline=self.baseline)
        self.corrected_signal = self.make_corrected_signal(
            signal=self.signal, baseline=self.baseline
        )

    def baseline_correction(self, array):
        signal_arr = array
        return misc.beads(
            signal_arr,
            freq_cutoff=BASELINE_FREQ_CUTOFF,
            lam_0=BASELINE_LAMBDAS[0],
            lam_1=BASELINE_LAMBDAS[1],
            lam_2=BASELINE_LAMBDAS[2],
            asymmetry=BASELINE_ASYMMETRY,
        )

    def save_baseline_to_csv(self, signal, baseline):
        baseline_export = pd.DataFrame(baseline[0], columns=["I, pA"])
        baseline_export.to_csv(signal.report_folder / "baseline.csv", index=True)

    def make_corrected_signal(self, signal, baseline):
        logger.info("Corrected signal calculation started.")
        corrected_signal = np.array(signal.raw_signal)[0] - baseline[0]
        return corrected_signal

    def sweeps_slices_dict_with_offcet(self, array, signal, offset=(0, 0)):
        return {
            k: array[
                bound_counting(k, signal, offset=offset[0]) : bound_counting(
                    k + 1, signal, offset=offset[1]
                )
            ]
            for k in range(signal.sweep_number)
        }

    def _make_inopeak_steps_mean(self, signal, baseline):
        return {
            k: baseline[0][
                bound_counting(k, signal, offset=0) : bound_counting(
                    k, signal, offset=0.9 * signal.step[1]
                )
            ].mean()
            for k in range(signal.sweep_number)
        }

    def series_resistance(self, signal, baseline):
        logger.info("Series resistance calculation started.")
        peak_ampl = {
            k: np.array(signal.raw_signal)[0][
                bound_counting(
                    k, signal, offset=signal.step[1] - 0.002
                ) : bound_counting(k, signal, offset=signal.step[1] + 0.02)
            ].min()
            for k in range(signal.sweep_number)
        }
        inopeak_steps_mean = self._make_inopeak_steps_mean(
            signal=signal, baseline=baseline
        )
        ipeak = {
            k: inopeak_steps_mean[k] - peak_ampl[k] for k in range(signal.sweep_number)
        }

        r_s_dict = {
            k: signal.step[0] * 10**3 / ipeak[k] for k in range(signal.sweep_number)
        }
        return r_s_dict

    def input_resistance(self, signal, baseline, r_s_dict):
        logger.info("Input resistance calculation started.")
        inopeak_steps_mean = self._make_inopeak_steps_mean(
            signal=signal, baseline=baseline
        )
        iss_steps = self.sweeps_slices_dict_with_offcet(
            array=baseline[0],
            signal=signal,
            offset=(signal.step[1] + 0.2, signal.step[1] + 0.45),
        )
        iss = {
            k: self.baseline_correction(iss_steps[k])
            for k in range(signal.sweep_number)
        }
        iss_mean = {
            k: inopeak_steps_mean[k] - iss[k][0].mean()
            for k in range(signal.sweep_number)
        }
        return {
            k: signal.step[0] * 10**3 / abs(iss_mean[k]) - r_s_dict[k]
            for k in range(signal.sweep_number)
        }

    def tonic_current_calc(self, signal, colors, distribution_param_sweep):
        logger.info("Tonic current calculation.")
        temp_tonic_cur = {k: [] for k in colors.keys()}
        for k in temp_tonic_cur.keys():
            for j in range(signal.sweep_number):
                if colors[k][0] <= j <= colors[k][1]:
                    temp_tonic_cur[k].append(distribution_param_sweep[j][0])

        mean_hold_cur = {
            k: np.array(temp_tonic_cur[k]).mean() for k in temp_tonic_cur.keys()
        }
        # ton_cur = pd.DataFrame(tonic_cur, columns=["Iton, pA"])
        # ton_cur.to_csv(report_folder / "_tonic_cur.csv", index=True)

        return [
            mean_hold_cur[k] - mean_hold_cur["CTRL"]
            for k in mean_hold_cur.keys()
            if k != "CTRL"
        ]

    def low_pass_filter(
        self, signal, corrected_signal, corrected_signal_sweeps_no_res_part
    ):
        logger.info("Signal filtering.")
        sig = corrected_signal
        sos = sg.butter(1, 50, "lp", fs=signal.sampling_rate, output="sos")
        filtered_signal = sg.sosfilt(sos, sig)  # filter whole signal

        filtered_signal_dict = {}  # filter sweeps
        for k in corrected_signal_sweeps_no_res_part.keys():
            sig = corrected_signal_sweeps_no_res_part[k]
            sos = sg.butter(1, 50, "lp", fs=signal.sampling_rate, output="sos")
            filtered_signal_dict[k] = sg.sosfilt(sos, sig)

        return filtered_signal, filtered_signal_dict

    def event_detection(self, signal, corrected_signal, filtered_signal, baseline):
        logger.info("Event detection started.")
        sigma = mad_std(corrected_signal)
        # Event detection in filtered signal
        block1 = signal.create_block()
        block1.segments[0].analogsignals[0].base[:] = signal.create_block().segments[
            0
        ].analogsignals[0].units * filtered_signal.reshape(-1, 1)
        events_times_reset = {}
        for i in range(signal.sweep_number):
            if i < signal.sweep_number - 1:
                a = block1.segments[0].time_slice(
                    (i * signal.sweep_length + 0.7) * pq.s,
                    (i * signal.sweep_length + signal.sweep_length) * pq.s,
                    reset_time=False,
                )
            else:
                a = block1.segments[0].time_slice(
                    (i * signal.sweep_length + 0.7) * pq.s, reset_time=False
                )
            Train_times = elephant.spike_train_generation.peak_detection(
                a.analogsignals[0],
                threshold=-sigma * 5 * pq.pA,
                sign="below",
                as_array=True,
            )
            events_times_reset[i] = Train_times
        events_inds_reset = {
            i: np.around(events_times_reset[i] * signal.sampling_rate).astype(int)
            for i in range(signal.sweep_number)
        }

        # Event detection in not filtered signal
        block2 = signal.create_block()
        block2.segments[0].analogsignals[0].base[:] = signal.create_block().segments[
            0
        ].analogsignals[0].units * baseline[0].reshape(-1, 1)
        block3 = signal.create_block()
        corrected_signal_neo = (
            block3.segments[0].analogsignals[0] - block2.segments[0].analogsignals[0]
        )
        events_times_cor = {}
        for i in range(signal.sweep_number):
            if i < signal.sweep_number - 1:
                a = corrected_signal_neo.time_slice(
                    (i * signal.sweep_length + signal.step[2]) * pq.s,
                    (i * signal.sweep_length + signal.sweep_length) * pq.s,
                )
            else:
                a = corrected_signal_neo.time_slice(
                    (i * signal.sweep_length + 0.7) * pq.s, None
                )
            Train_times = elephant.spike_train_generation.peak_detection(
                a, threshold=-sigma * 5 * pq.pA, sign="below", as_array=True
            )
            events_times_cor[i] = Train_times[:-1]

        events_inds_cor = {
            i: np.around(events_times_cor[i] * signal.sampling_rate).astype(int)
            for i in range(signal.sweep_number)
        }

        spike_train_ind = defaultdict(
            list
        )  # correspondence between filtered and not filtered signals
        for k in range(signal.sweep_number):
            prevbel = 0
            prev_ind = 0
            for rel in events_inds_reset[k]:
                for i in range(prev_ind, len(events_inds_cor[k])):
                    bel = events_inds_cor[k][i]
                    if bel > rel:
                        spike_train_ind[k].append(prevbel)
                        prev_ind = i - 1
                        break
                    prevbel = bel

        return spike_train_ind

    # noinspection PyTupleAssignmentBalance
    def half_width(self, signal, corrected_signal, spike_train_ind):
        prominences, left_bases, right_bases = peak_prominences(
         corrected_signal, spike_train_ind
        )
        offset = np.ones_like(prominences)
        widths_ind, h_eval, left_ips, right_ips = peak_widths(
            -1 * corrected_signal, spike_train_ind, rel_height=0.5
        )
        right_ips_times = right_ips / signal.sampling_rate
        widths = widths_ind / signal.sampling_rate
        return widths, right_ips_times

    def spike_selection(self, corrected_signal, signal, spike_train_ind):
        logger.info("Event selection started.")
        Sp_for_ampl = {k: [] for k in range(signal.sweep_number)}
        for k in range(signal.sweep_number):
            for i in range(len(spike_train_ind[k])):
                Sp_for_ampl[k].append(int(spike_train_ind[k][i]))
        ampls_spks = {
            k: corrected_signal[Sp_for_ampl[k]] for k in range(signal.sweep_number)
        }
        widths = {}  # 5ms <= tau_decay <= 40 mc
        right_ips_times = {}
        for k in range(signal.sweep_number):
            widths[k], right_ips_times[k] = self.half_width(
                signal=signal,
                corrected_signal=corrected_signal,
                spike_train_ind=spike_train_ind[k],
            )

        tau_decay = {
            k: np.zeros(len(spike_train_ind[k])) for k in range(signal.sweep_number)
        }
        for k in range(signal.sweep_number):
            if len(spike_train_ind[k]) != 0:
                for l in range(len(tau_decay[k])):
                    x = []
                    y = []
                    Amax = ampls_spks[k][l]
                    for i in range(
                        spike_train_ind[k][l],
                        int(right_ips_times[k][l] * signal.sampling_rate) + 2,
                    ):
                        if (
                            corrected_signal[i] / Amax < 1
                            and corrected_signal[i] / Amax > 0
                        ):
                            y.append(corrected_signal[i] / Amax)
                            x.append(i / signal.sampling_rate)
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
        new_spike_train_ind = {}
        for k in tau_decay.keys():
            new_spike_train_ind[k] = []
            for i in range(len(tau_decay[k])):
                if 5 <= tau_decay[k][i] <= 40:
                    new_spike_train_ind[k].append(spike_train_ind[k][i])
                else:
                    continue
        return new_spike_train_ind

    def frequency(
        self, signal, new_spike_train_ind, corrected_signal_sweeps_no_res_part
    ):
        logger.info("Frequencies calculation.")
        time_len_sweeps = {
            k: len(corrected_signal_sweeps_no_res_part) * signal.sampling_rate
            for k in range(signal.sweep_number)
        }
        freqs = {
            k: len(new_spike_train_ind[k]) * 1000 / time_len_sweeps[k]
            for k in range(signal.sweep_number)
        }

        spont_freqs = {
            k: 1 / elephant.statistics.isi(new_spike_train_ind[k]) * signal.sampling_rate
            for k in range(signal.sweep_number)
        }
        return freqs, spont_freqs

    def amplitude(self, signal, corrected_signal, new_spike_train_ind):
        logger.info("Amplitudes calculation...")
        Sp_for_ampl = {k: [] for k in range(signal.sweep_number)}
        for k in range(signal.sweep_number):
            for i in range(len(new_spike_train_ind[k])):
                Sp_for_ampl[k].append(int(new_spike_train_ind[k][i]))

        ampls_spks = {
            k: corrected_signal[Sp_for_ampl[k]] for k in range(signal.sweep_number)
        }

        ampls_spks_pos = {}
        for k in range(signal.sweep_number):
            ampls_spks_pos[k] = ampls_spks[k] * (-1)

        SpikeTrainTimes = {k: [] for k in range(signal.sweep_number)}
        for k in range(signal.sweep_number):
            SpikeTrainTimes[k] = np.array(new_spike_train_ind[k]) / signal.sampling_rate
        return ampls_spks_pos, SpikeTrainTimes, ampls_spks

    def calculate_tau_decay(self, signal, corrected_signal, new_spike_train_ind, ampls_spks):
        logger.info("Tau decays calculation...")
        widths = {}
        right_ips_times = {}
        for k in range(signal.sweep_number):
            widths[k], right_ips_times[k] = self.half_width(signal=self.signal,
                                                            corrected_signal=corrected_signal,
                                                            spike_train_ind=new_spike_train_ind[k])

        tau_decay = {
            k: np.zeros(len(new_spike_train_ind[k])) for k in range(signal.sweep_number)
        }
        for k in range(signal.sweep_number):
            if len(new_spike_train_ind[k]) != 0:
                for l in range(len(tau_decay[k])):
                    x = []
                    y = []
                    Amax = ampls_spks[k][l]
                    for i in range(
                        new_spike_train_ind[k][l],
                        int(right_ips_times[k][l] * signal.sampling_rate) + 2,
                    ):
                        if (
                            corrected_signal[i] / Amax < 1
                            and corrected_signal[i] / Amax > 0
                        ):
                            y.append(corrected_signal[i] / Amax)
                            x.append(i / signal.sampling_rate)
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
        return tau_decay

    def rise_rate10_90(self, signal, final_spike_train):
        rise_time_10_90 = []
        widths_ind_10, h_eval_10, left_ips_10, right_ips_10 = peak_widths(
            self.corrected_signal, final_spike_train, rel_height=0.9
        )
        widths_ind_90, h_eval_90, left_ips_90, right_ips_90 = peak_widths(
            self.corrected_signal, final_spike_train, rel_height=0.1
        )
        left_ips_times_10 = left_ips_10 / signal.sampling_rate
        left_ips_times_90 = left_ips_90 / signal.sampling_rate
        rise_time_10_90.append(
            (h_eval_90 - h_eval_10)
            * 10 ** (-3)
            / (left_ips_times_90 - left_ips_times_10)
        )

        return rise_time_10_90

    def rise_rate(self, signal, final_spike_train):
        rise_time_10_90 = []
        widths_ind_10, h_eval_10, left_ips_10, right_ips_10 = peak_widths(
            self.corrected_signal * (-1), final_spike_train, rel_height=0.9)
        widths_ind_90, h_eval_90, left_ips_90, right_ips_90 = peak_widths(
            self.corrected_signal * (-1), final_spike_train,
            rel_height=0.1)
        left_ips_times_10 = left_ips_10 / signal.sampling_rate
        left_ips_times_90 = left_ips_90 / signal.sampling_rate
        rise_time_10_90.append((h_eval_90 - h_eval_10) * 10 ** (-3) / (left_ips_times_90 - left_ips_times_10))

        return rise_time_10_90

    def calculate_rise_rate(self, signal, new_spike_train_ind):
        logger.info("Rise rate calculation.")
        rise_rate = {
            k: self.rise_rate(signal, new_spike_train_ind[k])
            for k in range(signal.sweep_number)
        }
        for k in range(signal.sweep_number):
            rise_rate[k] = np.array(rise_rate[k])[0]

        return rise_rate

    def calculate_statistics(self):
        corrected_signal_sweeps = self.sweeps_slices_dict_with_offcet(
            signal=self.signal, array=self.corrected_signal, offset=(0, 0)
        )
        corrected_signal_sweeps_no_res_part = self.sweeps_slices_dict_with_offcet(
            signal=self.signal,
            array=self.corrected_signal,
            offset=(self.signal.step[2], 0),
        )
        self.r_s_dict = self.series_resistance(
            signal=self.signal, baseline=self.baseline
        )
        self.r_in_dict = self.input_resistance(
            signal=self.signal, baseline=self.baseline, r_s_dict=self.r_s_dict
        )

        baseline_sweeps_no_res_part = self.sweeps_slices_dict_with_offcet(
            signal=self.signal, array=self.baseline[0], offset=(self.signal.step[2], 0)
        )

        self.distribution_param_sweep = {
            k: [
                baseline_sweeps_no_res_part[k].mean(),
                baseline_sweeps_no_res_part[k].std(),
            ]
            for k in range(self.signal.sweep_number)
        }

        self.tonic_cur = self.tonic_current_calc(
            signal=self.signal,
            colors=self.colors,
            distribution_param_sweep=self.distribution_param_sweep,
        )
        filtered_signal, filtered_signal_dict = self.low_pass_filter(
            corrected_signal=self.corrected_signal,
            signal=self.signal,
            corrected_signal_sweeps_no_res_part=corrected_signal_sweeps_no_res_part,
        )
        spike_train_ind = self.event_detection(
            signal=self.signal,
            corrected_signal=self.corrected_signal,
            filtered_signal=filtered_signal,
            baseline=self.baseline,
        )
        self.new_spike_train_ind = self.spike_selection(
            corrected_signal=self.corrected_signal,
            signal=self.signal,
            spike_train_ind=spike_train_ind,
        )
        self.freqs, self.spont_freq = self.frequency(
            signal=self.signal,
            new_spike_train_ind=self.new_spike_train_ind,
            corrected_signal_sweeps_no_res_part=corrected_signal_sweeps_no_res_part,
        )
        self.ampls_spks_pos, self.SpikeTrainTimes, ampls_spks = self.amplitude(
            signal=self.signal,
            corrected_signal=self.corrected_signal,
            new_spike_train_ind=self.new_spike_train_ind,
        )

        self.tau_decay = self.calculate_tau_decay(
            signal=self.signal, corrected_signal=self.corrected_signal, new_spike_train_ind=self.new_spike_train_ind, ampls_spks=ampls_spks
        )
        self.rise_rate = self.calculate_rise_rate(
            signal=self.signal, new_spike_train_ind=self.new_spike_train_ind
        )
