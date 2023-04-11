import numpy as np
import pandas as pd
from scipy import stats


def get_t_scores(df):
    """
    Perform a t test to check for difference in means between variants of concern and non-concern
    :param df:
    :return:
    """
    # set option to displa all columns
    pd.set_option('display.max_columns', None)
    df = df.copy()
    # get the means
    means = df.groupby('variant_of_concern').mean()
    # get the standard deviations
    stds = df.groupby('variant_of_concern').std()
    # get the counts
    counts = df.groupby('variant_of_concern').count()
    # get the degrees of freedom
    print("\nMEANS ")
    print(means)
    print("\nSTDS")
    print("\nSCORES")
    print(stds)
    d_f = counts - 1
    # get the t scores
    t_scores = (means.iloc[0] - means.iloc[1]) / (stds.iloc[0] / np.sqrt(counts.iloc[1]))
    # get the p values as a series based on the t scores
    # make a new dataframe stats with t_scores as a column
    stats_ = pd.DataFrame({'t_scores': t_scores})
    stats_['p_values'] = 2 * (1 - stats.t.cdf(stats_['t_scores'], df=d_f.iloc[0]))

    # # get the p values
    # p_values = 2 * (1 - stats.t.cdf(t_scores, df=d_f))
    # get the confidence intervals
    # ci = stats.t.interval(0.95, df=d_f, loc=means['mean'], scale=stds['mean'] / np.sqrt(counts['mean']))
    # create a dataframe to return
    # output = pd.DataFrame({'t_scores': t_scores, 'p_values': p_values})
    # return the dataframe
    return stats_

