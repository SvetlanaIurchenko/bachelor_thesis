import pandas as pd
from functools import partial
from src import settings
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.tools.report_functions import sweep_seria_correspondence


def read_electric_cell_parameters(files, colors, report_table_name, col_1, col_2):
    electric_param = {}
    r_s = {}
    r_in = {}
    for f in files:
        electric_param[f] = pd.read_csv(settings.REPORTS_DIR_RECORDING / f / f"{report_table_name}.csv")
        r_s[f] = np.array(electric_param[f][col_1]) # 'R_s, Mohm'
        r_in[f] = np.array(electric_param[f][col_2]) # 'R_in, Mohm'

    dict_df_r_s = {}
    for k in r_s.keys():
        dict_df_r_s[k] = pd.DataFrame.from_dict(r_s[k])

    dict_df_r_in = {}
    for k in r_in.keys():
        dict_df_r_in[k] = pd.DataFrame.from_dict(r_in[k])

    sweep_seria_correspondence_fun = partial(
        sweep_seria_correspondence, colors=colors
    )

    for k, v in dict_df_r_s.items():
        v['sweep_number'] = [i for i in range(len(r_s[k]))]
        v['series'] = v.sweep_number.apply(sweep_seria_correspondence_fun)

    for k, v in dict_df_r_in.items():
        v['sweep_number'] = [i for i in range(len(r_in[k]))]
        v['series'] = v.sweep_number.apply(sweep_seria_correspondence_fun)

    return dict_df_r_s, dict_df_r_in

def common_statistic(colors, dict_df_r_s):
    r_s_mean = {}
    for k in dict_df_r_s.keys():
        r_s_mean_temp = []
        for c in colors.keys():
            r_s_mean_temp.append(dict_df_r_s[k].loc[dict_df_r_s[k]['series'] == c][0].mean())
        r_s_mean[k] = r_s_mean_temp
    return r_s_mean

def create_table_for_boxplot(r_s_common_stat, colors, value_name):
    r_s_common_stats_res_df = pd.DataFrame(columns=['filename', 'mice_type', 'mice_age'] + list(colors.keys()))
    for i, (key, val) in enumerate(r_s_common_stat.items()):
        mice_type = key.split()[2].upper()
        mice_age = key.split()[3]
        r_s_common_stats_res_df.loc[i] = [key, mice_type, mice_age] + val
    r_s_final = r_s_common_stats_res_df.melt(id_vars=['filename', 'mice_type', 'mice_age'],
                            value_vars=list(colors.keys()),
                            var_name='series',
                            value_name=value_name)
    return r_s_final

def create_box_plot_with_hue(report_folder, data, x,  y, hue, **kwargs):
    plt.figure(figsize=(10, 10))
    ax = sns.boxplot(data=data, x=x, y=y, hue=hue, linewidth=1.3, showfliers=False)
    ax.set_xlabel('series', fontsize=18)
    ax.set_ylabel(y, fontsize=18)
    plt.xticks(fontsize=12, rotation=0)
    plt.yticks(fontsize=12, rotation=0)
    plt.legend(fontsize=18)
    plt.title('plotname', fontsize=18)
    plt.savefig(report_folder / f"{kwargs.get('plotname')}.png")
    plt.savefig(report_folder / f"{kwargs.get('plotname')}.png")
    plt.close()

def tonic_current_common_statistic(files, report_table_name, colors):
    tonic_cur = {}
    for f in files:
        tonic_cur[f] = pd.read_csv(settings.REPORTS_DIR_RECORDING / f / f"{report_table_name}.csv")

    dict_df_tonic_cur = {}
    for k in tonic_cur.keys():
        dict_df_tonic_cur[k] = list(tonic_cur[k]['Iton, pA'])

    tonic_cur_common_stats_res_df = pd.DataFrame(
        columns=['filename', 'mice_type', 'mice_age'] + list(colors.keys())[1:])
    for i, (key, val) in enumerate(dict_df_tonic_cur.items()):
        mice_type = key.split()[2].upper()
        mice_age = key.split()[3]
        tonic_cur_common_stats_res_df.loc[i] = [key, mice_type, mice_age] + val

    tonic_cur_final = tonic_cur_common_stats_res_df.melt(id_vars=['filename', 'mice_type', 'mice_age'],
                                                         value_vars=list(colors.keys())[1:],
                                                         var_name='series',
                                                         value_name='Iton, pA')
    return tonic_cur_final

def read_event_parameters(files, report_table_name):
    event_param = {}
    for f in files:
        event_param[f] = pd.read_csv(settings.REPORTS_DIR_RECORDING / f / f"{report_table_name}.csv")
    return event_param

def event_param_processing(column, event_param, colors):
    amplitudes = {}
    for k in event_param.keys():
        ampl_temp = []
        for c in list(colors.keys())[:3]:
            ampl_temp.append(np.array(event_param[k].loc[event_param[k]['series'] == c][column]).mean())
        amplitudes[k] = ampl_temp
    return amplitudes

def create_table_for_boxplot_event_params(r_s_common_stat, colors, value_name):
    r_s_common_stats_res_df = pd.DataFrame(columns=['filename', 'mice_type', 'mice_age'] + list(colors.keys())[:3])
    for i, (key, val) in enumerate(r_s_common_stat.items()):
        mice_type = key.split()[2].upper()
        mice_age = key.split()[3]
        r_s_common_stats_res_df.loc[i] = [key, mice_type, mice_age] + val[:3]
    r_s_final = r_s_common_stats_res_df.melt(id_vars=['filename', 'mice_type', 'mice_age'],
                            value_vars=list(colors.keys())[:3],
                            var_name='series',
                            value_name=value_name)
    return r_s_final

def read_spontaneous_frequency(files, report_table_name):
    spont_freqs = {}
    for f in files:
        spont_freqs[f] = pd.read_csv(settings.REPORTS_DIR_RECORDING / f / f"{report_table_name}.csv").rename(columns={'spont_freq': 0})
    return spont_freqs


