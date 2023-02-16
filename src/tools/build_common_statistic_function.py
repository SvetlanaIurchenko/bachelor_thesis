import pandas as pd
from functools import partial
from src import settings
import matplotlib.patheffects as path_effects
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.tools.report_functions import sweep_seria_correspondence
from scipy import stats
from itertools import product


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
    r_s_ctrl = r_s_final[r_s_final.series == 'CTRL'].rename(columns={value_name: value_name + '_ctrl'}).drop(columns='series')
    r_s_new_final = r_s_final.merge(r_s_ctrl, on=['filename', 'mice_type', 'mice_age'])
    r_s_new_final['delta_' + value_name] = (r_s_new_final[value_name] - r_s_new_final[value_name + '_ctrl']) * 100 / r_s_new_final[value_name + '_ctrl']

    return r_s_final, r_s_new_final

def add_median_labels(ax, fmt='.1f'):
    lines = ax.get_lines()
    boxes = [c for c in ax.get_children() if type(c).__name__ == 'PathPatch']
    lines_per_box = int(len(lines) / len(boxes))
    for median in lines[4:len(lines):lines_per_box]:
        x, y = (data.mean() for data in median.get_data())
        # choose value depending on horizontal or vertical plot orientation
        value = x if (median.get_xdata()[1] - median.get_xdata()[0]) == 0 else y
        text = ax.text(x, y, f'{value:{fmt}}', ha='center', va='center',
                       fontweight='bold', color='white')
        # create median-colored border around white text for contrast
        text.set_path_effects([
            path_effects.Stroke(linewidth=3, foreground=median.get_color()),
            path_effects.Normal(),
        ])

def create_box_plot_with_hue(report_folder, data, x,  y, hue, **kwargs):
    plt.figure(figsize=(10, 10))
    ax = sns.boxplot(data=data, x=x, y=y, hue=hue, linewidth=1.3, showfliers=False)
    add_median_labels(ax)
    ax.set_xlabel('series', fontsize=18)
    ax.set_ylabel(y, fontsize=18)
    plt.xticks(fontsize=12, rotation=0)
    plt.yticks(fontsize=12, rotation=0)
    plt.legend(fontsize=18)
    plt.title(f"{kwargs.get('plotname')}", fontsize=18 )
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
    r_s_ctrl = r_s_final[r_s_final.series == 'CTRL'].rename(columns={value_name: value_name + '_ctrl'}).drop(columns='series')
    r_s_new_final = r_s_final.merge(r_s_ctrl, on=['filename', 'mice_type', 'mice_age'])
    r_s_new_final['delta_' + value_name] = (r_s_new_final[value_name] - r_s_new_final[value_name + '_ctrl']) * 100 / r_s_new_final[value_name + '_ctrl']
    return r_s_final, r_s_new_final

def read_spontaneous_frequency(files, report_table_name):
    spont_freqs = {}
    for f in files:
        spont_freqs[f] = pd.read_csv(settings.REPORTS_DIR_RECORDING / f / f"{report_table_name}.csv").rename(columns={'spont_freq': 0})
    return spont_freqs


def create_result_table_statistic_tests(report_folder, r_s_final, y, **kwargs):
    df = r_s_final.groupby(['mice_type', 'mice_age', 'series'])[y].apply(list).reset_index()
    itteration_array = list(product(df.index, df.index))
    itteration_array_names, t_tests_res, wilcoxon_res, mann_whitney_res, k_s_res = [], [], [], [], []
    for el in itteration_array:
        temp = df.mice_type[el[0]] + df.mice_age[el[0]] + df.series[el[0]] + ' vs ' + df.mice_type[el[1]] + df.mice_age[
            el[1]] + df.series[el[1]]
        itteration_array_names.append(temp)
        a = df[y][el[0]]
        b = df[y][el[1]]
        t_tests_res.append(stats.ttest_ind(a, b)[1])
        mann_whitney_res.append(stats.mannwhitneyu(a, b)[1])
        k_s_res.append(stats.kstest(a, b)[1])
        if el[0] != el[1] and len(a) == len(b):
            wilcoxon_res.append(stats.wilcoxon(a, b)[1])
        else:
            wilcoxon_res.append('Not calculated')
    statistic_tests_results = pd.DataFrame(itteration_array_names, columns=['pairs'])
    statistic_tests_results['student_t_test, p'] = t_tests_res
    statistic_tests_results['wilcoxon_test, p'] = wilcoxon_res
    statistic_tests_results['mann_whitneyu_test, p'] = mann_whitney_res
    statistic_tests_results['k_s_test, p'] = k_s_res
    statistic_tests_results.to_csv(report_folder / f"{kwargs.get('plotname')}.csv", index=True)

