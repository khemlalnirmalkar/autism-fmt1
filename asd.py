import seaborn as sns
import scipy.stats
import bisect
import numpy as np
import pandas as pd

from skbio.stats.distance import MissingIDError

palette = {'Responder': '#68228B', 'Non-responder': '#807dba',
           'neurotypical': '#228B22', 'donor': '#d95f02',
           'Rectal': '#68228B', 'Oral': '#807dba'}

def get_donor_sids(sample_md):
    donor_sids = sample_md[sample_md['SampleType'] == 'donor-stool']['SubjectID']
    donor_sids = dict((v, k) for k, v in donor_sids.items())
    return donor_sids

def control_distance_to_donors(sample_md, dm, sample_type):
    donor_sids = get_donor_sids(sample_md).values()
    group_sids = filter_sample_md(sample_md,
                                  [('Group', 'neurotypical'), ('time_point', "1"),
                                   ('SampleType', sample_type)]).index
    results = []
    for donor_sid in donor_sids:
        for group_sid in group_sids:
            try:
                d = dm[donor_sid, group_sid]
            except MissingIDError:
                continue
            results.append(d)
    return results

def donor_metric(sample_md, metric):
    group_sids = filter_sample_md(sample_md,
                                  [('Group', 'donor-initial'), ('SampleType', 'donor-stool')]).index
    return group_metric(sample_md, group_sids, metric)

def control_metric(sample_md, sample_type, metric):
    group_sids = filter_sample_md(sample_md,
                                  [('Group', 'neurotypical'), ('time_point', "1"), ('SampleType', sample_type)]).index
    return group_metric(sample_md, group_sids, metric)


def group_metric(sample_md, group_sids, metric):
    results = []
    for group_sid in group_sids:
        try:
            v = float(sample_md[metric][group_sid])
        except ValueError:
            continue
        if not np.isnan(v):
            results.append(v)
    return results

def _distance_to_donor(sid, donor, dm, donor_sids):
    try:
        donor_sid = donor_sids[donor]
    except KeyError:
        return np.nan

    try:
        return dm[donor_sid, sid]
    except MissingIDError:
        return np.nan

def add_distance_to_donor(sample_md, dm, metric_name, verbose=False):
    distance_to_initial_donor_data = []
    distance_to_most_relevant_donor_data = []
    donor_sids = get_donor_sids(sample_md)
    for sid in sample_md.index:
        if sample_md['Group'][sid] == 'autism':
            # compile distance to initial donor
            initial_donor = sample_md['BBT_Donor_ID'][sid]
            distance_to_initial_donor_data.append(_distance_to_donor(sid, initial_donor, dm, donor_sids))

            # compile distance to most relevant donor (either the most recent, or the initial donor
            # if pre-FMT, which helps to evaluate initial change)
            most_relevant_donor = sample_md['MostRecentDonorID'][sid]
            if most_relevant_donor == 'pre-treatment':
                most_relevant_donor = sample_md['BBT_Donor_ID'][sid]
            distance_to_most_relevant_donor_data.append(_distance_to_donor(sid, most_relevant_donor, dm, donor_sids))
        elif sample_md['Group'][sid] == 'neurotypical' and sample_md['time_point'][sid] == "1":
            # compile median distance from timepoint 0 of control to all donor samples
            distances_to_donors = []
            for donor_sid in donor_sids.values():
                try:
                    distances_to_donors.append(dm[donor_sid, sid])
                except MissingIDError:
                    continue
            median_distance_to_donors = np.median(distances_to_donors)
            distance_to_initial_donor_data.append(np.median(median_distance_to_donors))
            distance_to_most_relevant_donor_data.append(np.median(median_distance_to_donors))
        else:
            # not a distance that we're compiling
            distance_to_initial_donor_data.append(np.nan)
            distance_to_most_relevant_donor_data.append(np.nan)

    sample_md['%s to initial donor' % metric_name] = distance_to_initial_donor_data
    sample_md['%s to most relevant donor' % metric_name] = distance_to_most_relevant_donor_data
    return sample_md

def filter_sample_md(sample_md, includes):
    result = sample_md.copy()
    for column, value in includes:
        result = result[result[column] == value]
    return result

def plot_week_data(df, sample_type, metric, hue=None, hide_donor_baseline=False, hide_control_baseline=False):
    df['week'] = pd.to_numeric(df['week'], errors='coerce')
    df[metric] = pd.to_numeric(df[metric], errors='coerce')
    asd_data = filter_sample_md(df, [('SampleType', sample_type), ('Group', 'autism')])
    asd_data = asd_data.sort_values(by='week')
    ax = sns.boxplot(data=asd_data, x='week', y=metric, color='white')
    ax = sns.swarmplot(data=asd_data, x='week', y=metric, hue=hue, palette=palette)
    control_y = np.median(control_metric(df, sample_type, metric=metric))
    donor_y = np.median(donor_metric(df, metric=metric))
    x0 = np.min(df['week']) - 1
    x1 = np.max(df['week']) + 1
    if not hide_control_baseline:
        ax.plot([x0, x1], [control_y, control_y],
                color=palette['neurotypical'], linestyle='--', label='neurotypical (median)')
    if not hide_donor_baseline:
        ax.plot([x0, x1], [donor_y, donor_y],
            color=palette['donor'], linestyle='--', label='donor (median)')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    return ax

alphas = [(0.001, '***'), (0.01, '**'), (0.05, '*')]

def get_sig_text(p, alphas, null_text=""):
    if np.isnan(p):
        return null_text
    alphas.sort()
    if p >= alphas[-1][0]:
        return 'ns'
    sorted_location = bisect.bisect([e[0] for e in alphas], p)
    return alphas[sorted_location][1]

def plot_week_data_with_stats(sample_md, sample_type, metric, hue=None, alphas=alphas, one_tailed=False,
                              hide_donor_baseline=False, hide_control_baseline=False):
    ax = plot_week_data(sample_md, sample_type, metric, hue, hide_donor_baseline=hide_donor_baseline,
                        hide_control_baseline=hide_control_baseline)
    stats = tabulate_week_to_week0_paired_stats(sample_md, sample_type, metric)
    ymax = ax.get_ylim()[1]
    stats.sort_index()
    for i, w in enumerate(stats.index):
        t, p = stats['t'][w], stats['p-value'][w]
        if one_tailed:
            if t < 0.0:
                p = p/2.
            else:
                p = 1.0
        sig_text = get_sig_text(p, alphas)
        ax.text(i, 1.02*ymax, sig_text, ha='center',va='center')
    return ax

# Paired t-test: change in distance to donor from time zero is different than zero
# (positive t means more different than donor, negative t means more similar to donor)
def tabulate_week_to_week0_paired_stats(df, sample_type, metric):
    asd_data = filter_sample_md(df, [('SampleType', sample_type), ('Group', 'autism')])
    results = []
    asd_data['week'] = pd.to_numeric(asd_data['week'], errors='coerce')
    asd_data = asd_data.sort_values(by='week')
    weeks = asd_data['week'].unique()
    for i in weeks:
        g = asd_data[np.logical_or(asd_data['week'] == 0, asd_data['week'] == i)].groupby('SubjectID')
        paired_diffs = g.diff()[metric].dropna()
        t, p = scipy.stats.ttest_1samp(paired_diffs, popmean=0)
        results.append((len(paired_diffs), np.median(paired_diffs), t, p))
    return pd.DataFrame(results, index=pd.Index(weeks, name='week'),
                        columns=['n', metric, 't', 'p-value'])

# two sample t-test: does ASD data at each week differ significantly from
# the control group?
def tabulate_week_to_control_stats(df, sample_type, metric):
    control_week0 = control_metric(df, sample_type, metric=metric)
    asd_data = filter_sample_md(df, [('SampleType', sample_type), ('Group', 'autism')])
    results = []
    asd_data['week'] = pd.to_numeric(asd_data['week'], errors='coerce')
    asd_data = asd_data.sort_values(by='week')
    weeks = asd_data['week'].unique()
    for i in weeks:
        weeki = asd_data[metric][asd_data['week'] == i].dropna()
        t, p = scipy.stats.ttest_ind(weeki, control_week0, equal_var=False)
        results.append((len(weeki), np.median(weeki), t, p))
    return pd.DataFrame(results, index=pd.Index(weeks, name='week'),
                        columns=['n', metric, 't', 'p-value'])
