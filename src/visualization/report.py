from src.tools.report_functions import (
    electrical_param_table,
    baseline_report,
    make_result_feature_table,
    tonic_cur_report,
make_result_table_spont_freq,
)
from src.visualization.plots import (
    plot_electric_params,
    plt_I_hold,
    plt_event_param,
    box_plot,
    plt_box_plot,
)


def make_report(signal, signal_features):
    plot_electric_params(
        signal_features.r_s_dict,
        signal,
        signal_features,
        **{"xlabel": "Time, s", "ylabel": "R_s, MOhm", "plotname": "rs"}
    )
    plot_electric_params(
        signal_features.r_in_dict,
        signal,
        signal_features,
        **{
            "xlabel": "Time, s",
            "ylabel": "R_in, MOhm",
            "ylim": (0, 1200),
            "plotname": "rin",
        }
    )

    plt_I_hold(signal, signal_features, is_noise=False, **{"plotname": "Ihold"})
    plt_I_hold(signal, signal_features, is_noise=True, **{"plotname": "Inoise"})

    plt_event_param(
        signal_features.ampls_spks_pos,
        signal,
        signal_features,
        **{"xlabel": "Time, s", "ylabel": "Amplitude, pA", "plotname": "ampl"}
    )
    plt_event_param(
        signal_features.tau_decay,
        signal,
        signal_features,
        **{"xlabel": "Time, s", "ylabel": "Tau decay, ms", "plotname": "taudec"}
    )
    plt_event_param(
        signal_features.rise_rate,
        signal,
        signal_features,
        **{"xlabel": "Time, s", "ylabel": "Rise rate, pA/ms", "plotname": "riserate"}
    )

    tonic_cur_report(signal=signal, signal_features=signal_features)
    electrical_param_table(signal=signal, signal_features=signal_features)
    baseline_report(signal=signal, signal_features=signal_features)

    result = make_result_feature_table(signal=signal, signal_features=signal_features)
    spont_freq = make_result_table_spont_freq(signal=signal, signal_features=signal_features)

    box_plot(
        signal=signal,
        data=spont_freq,
        y="spont_freq",
        **{"title": "Spont freqs", "ylabel": "Spont freq, Hz", "plotname": "spont freq"}
    )

    plt_box_plot(
        parametres=["amplitudes", "tau_decay", "rise_rate_10_90"],
        lparametres=["amplitudes, pA", "tau_decay, sec", "rise_rate_10_90, pA/ms"],
        res=result,
        signal=signal,
    )
