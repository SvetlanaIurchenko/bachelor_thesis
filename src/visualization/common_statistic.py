from src.tools.build_common_statistic_function import *

def make_common_statistic_plots(report_folder, r_s_final, r_s_new_final, y1, y2, units):
    create_box_plot_with_hue(report_folder, r_s_new_final[r_s_new_final.mice_age == '8m'].loc[r_s_new_final[r_s_new_final.mice_age == '8m']['series'] != 'CTRL']
                             , 'series', y1, 'mice_type', **{'plotname': y1 + ' 7-8 m.o.', "ylabel": y1 + ", %"})
    create_box_plot_with_hue(report_folder, r_s_new_final[r_s_new_final.mice_age == '4m'].loc[r_s_new_final[r_s_new_final.mice_age == '4m']['series'] != 'CTRL']
                             , 'series', y1, 'mice_type', **{'plotname': y1 + ' 3-4 m.o.', "ylabel": y1 + ", %"})
    create_box_plot_with_hue(report_folder, r_s_final.loc[r_s_final['series'] == 'CTRL']
                             ,'mice_type', y2, 'mice_age', **{'plotname': y2 + ' CTRL series', "ylabel": y2 + ',' + units})

def electric_parameters_group_common_statistic(report_folder, files, colors):
    dict_df_r_s, dict_df_r_in = read_electric_cell_parameters(files, colors, 'electric_cell_param','R_s, Mohm', 'R_in, Mohm')
    r_s_common_stat = common_statistic(colors, dict_df_r_s)
    r_s_final, r_s_new_final = create_table_for_boxplot(r_s_common_stat, colors, 'R_s')
    make_common_statistic_plots(report_folder, r_s_final, r_s_new_final, 'delta_R_s', 'R_s', ' MOhm')
    create_result_table_statistic_tests(report_folder, r_s_final, 'R_s', **{'plotname': 'stat_tests_for_series_resistance'})

    r_in_common_stat = common_statistic(colors, dict_df_r_in)
    r_in_final, r_in_new_final = create_table_for_boxplot(r_in_common_stat, colors, 'R_in')
    make_common_statistic_plots(report_folder, r_in_final, r_in_new_final, 'delta_R_in', 'R_in', ' MOhm')
    create_result_table_statistic_tests(report_folder, r_in_final, 'R_in', **{'plotname': 'stat_tests_for_input_resistance'})

def baseline_analysis_group_common_statistic(report_folder, files, colors):
    dict_df_i_h, dict_df_i_noise = read_electric_cell_parameters(files, colors, 'baseline_param', 'Ihold, pA', 'Inoise, pA')
    i_h_common_stat = common_statistic(colors, dict_df_i_h)
    i_h_final, i_h_new_final = create_table_for_boxplot(i_h_common_stat, colors, 'I_h')
    make_common_statistic_plots(report_folder, i_h_final, i_h_new_final, 'delta_I_h', 'I_h', ' pA')
    create_result_table_statistic_tests(report_folder, i_h_final, 'I_h', **{'plotname': 'stat_tests_for_i_hold'})

    i_noise_common_stat = common_statistic(colors, dict_df_i_noise)
    i_noise_final, i_noise_new_final = create_table_for_boxplot(i_noise_common_stat, colors, 'I_noise')
    make_common_statistic_plots(report_folder, i_noise_final, i_noise_new_final, 'delta_I_noise', 'I_noise', ' pA')
    create_result_table_statistic_tests(report_folder, i_noise_final, 'I_noise', **{'plotname': 'stat_tests_for_i_noise'})

def tonic_cur_group_common_statistic(report_folder, files, colors):
    tonic_cur_final = tonic_current_common_statistic(files, '_tonic_cur', colors)
    make_common_statistic_plots(report_folder, tonic_cur_final, 'tonic_cur', ' pA')
    create_result_table_statistic_tests(report_folder, tonic_cur_final, 'tonic_cur', **{'plotname': 'stat_tests_for_tonic_cur'})

def event_parameters_group_common_statistic(report_folder, files, colors):
    event_param = read_event_parameters(files, 'event_param')
    amplitudes = event_param_processing('amplitudes', event_param, colors)
    tau_decay = event_param_processing('tau_decay', event_param, colors)
    rise_rate_10_90 = event_param_processing('rise_rate_10_90', event_param, colors)
    amplitudes_final, amplitudes_new_final = create_table_for_boxplot_event_params(amplitudes, colors, 'amplitudes')
    tau_decay_final, tau_decay_new_final = create_table_for_boxplot_event_params(tau_decay, colors, 'tau_decay')
    rise_rate_10_90_final, rise_rate_10_90_new_final = create_table_for_boxplot_event_params(rise_rate_10_90, colors, 'rise_rate_10_90')

    make_common_statistic_plots(report_folder, amplitudes_final, amplitudes_new_final, 'delta_amplitudes', 'amplitudes', ' pA')
    make_common_statistic_plots(report_folder, tau_decay_final, tau_decay_new_final, 'delta_tau_decay', 'tau_decay', ' ms')
    make_common_statistic_plots(report_folder, rise_rate_10_90_final, rise_rate_10_90_new_final, 'delta_rise_rate_10_90', 'rise_rate_10_90', ' pA/ms')
    create_result_table_statistic_tests(report_folder, amplitudes_final, 'amplitudes', **{'plotname': 'stat_tests_for_amplitudes'})
    create_result_table_statistic_tests(report_folder, tau_decay_final, 'tau_decay', **{'plotname': 'stat_tests_for_tau_decay'})
    create_result_table_statistic_tests(report_folder, rise_rate_10_90_final, 'rise_rate_10_90', **{'plotname': 'stat_tests_for_rise_rate_10_90'})

    spont_freqs = read_spontaneous_frequency(files, 'spont_freqs')
    spont_freqs_mean = common_statistic(colors, spont_freqs)
    spont_freqs_final, spont_freqs_new_final = create_table_for_boxplot_event_params(spont_freqs_mean, colors, 'spont_freqs')
    make_common_statistic_plots(report_folder, spont_freqs_final, spont_freqs_new_final, 'delta_spont_freqs', 'spont_freqs', ' Hz')
    create_result_table_statistic_tests(report_folder, spont_freqs_final, 'spont_freqs', **{'plotname': 'stat_tests_for_spont_freqs'})













































































































