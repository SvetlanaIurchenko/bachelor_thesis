from src.tools.build_common_statistic_function import *

def make_common_statistic_plots(report_folder, r_s_final, y):
    create_box_plot_with_hue(report_folder, r_s_final[r_s_final.mice_age == '8m'].loc[r_s_final[r_s_final.mice_age == '8m']['series'] != 'CTRL']
                             , 'series', y, 'mice_type', **{'plotname': y + ' 7-8 m.o.'})
    create_box_plot_with_hue(report_folder, r_s_final[r_s_final.mice_age == '4m'].loc[r_s_final[r_s_final.mice_age == '4m']['series'] != 'CTRL']
                             , 'series', y, 'mice_type', **{'plotname': y + ' 3-4 m.o.'})
    create_box_plot_with_hue(report_folder, r_s_final.loc[r_s_final['series'] == 'CTRL']
                             ,'mice_type', y, 'mice_age', **{'plotname': y + ' CTRL series'})

def electric_parameters_group_common_statistic(report_folder, files, colors):
    dict_df_r_s, dict_df_r_in = read_electric_cell_parameters(files, colors, 'electric_cell_param','R_s, Mohm', 'R_in, Mohm')
    r_s_common_stat = common_statistic(colors, dict_df_r_s)
    r_s_final = create_table_for_boxplot(r_s_common_stat, colors, 'R_s')
    make_common_statistic_plots(report_folder, r_s_final, 'R_s')

    r_in_common_stat = common_statistic(colors, dict_df_r_in)
    r_in_final = create_table_for_boxplot(r_in_common_stat, colors, 'R_in')
    make_common_statistic_plots(report_folder, r_in_final, 'R_in')

def baseline_analysis_group_common_statistic(report_folder, files, colors):
    dict_df_i_h, dict_df_i_noise = read_electric_cell_parameters(files, colors, '_baseline_param', 'Ihold, pA', 'Inoise, pA')
    i_h_common_stat = common_statistic(colors, dict_df_i_h)
    i_h_final = create_table_for_boxplot(i_h_common_stat, colors, 'I_h')
    make_common_statistic_plots(report_folder, i_h_final, 'I_h')
    i_noise_common_stat = common_statistic(colors, dict_df_i_noise)
    i_noise_final = create_table_for_boxplot(i_noise_common_stat, colors, 'I_noise')
    make_common_statistic_plots(report_folder, i_noise_final, 'I_noise')

def tonic_cur_group_common_statistic(report_folder, files, colors):
    tonic_cur_final = tonic_current_common_statistic(files, '_tonic_cur', colors)
    make_common_statistic_plots(report_folder, tonic_cur_final, 'tonic_cur')

def event_parameters_group_common_statistic(report_folder, files, colors):
    event_param = read_event_parameters(files, 'event_param')
    amplitudes = event_param_processing('amplitudes', event_param, colors)
    tau_decay = event_param_processing('tau_decay', event_param, colors)
    rise_rate_10_90 = event_param_processing('rise_rate_10_90', event_param, colors)
    amplitudes_final = create_table_for_boxplot_event_params(amplitudes, colors, 'amplitudes')
    tau_decay_final = create_table_for_boxplot_event_params(tau_decay, colors, 'tau_decay')
    rise_rate_10_90_final = create_table_for_boxplot_event_params(rise_rate_10_90, colors, 'rise_rate_10_90')

    make_common_statistic_plots(report_folder, amplitudes_final, 'amplitudes')
    make_common_statistic_plots(report_folder, tau_decay_final, 'tau_decay')
    make_common_statistic_plots(report_folder, rise_rate_10_90_final, 'rise_rate_10_90')

    spont_freqs = read_spontaneous_frequency(files, 'spont_freqs')
    spont_freqs_mean = common_statistic(colors, spont_freqs)
    spont_freqs_final = create_table_for_boxplot_event_params(spont_freqs_mean, colors, 'spont_freqs')
    make_common_statistic_plots(report_folder, spont_freqs_final, 'spont_freqs')














































































































