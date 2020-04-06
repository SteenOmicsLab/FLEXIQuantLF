#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FLEXIQuant LF CLI version

This script:
- fits multiple robust linear regresion models (using RANSAC) for each sample with median reference sample peptide
intensities as independent variable and sample peptide intensities as dependent variable and selects the best model
- calculates vertical distances to the regression line for each peptide datapoint
- calculates raw scores by dividing the distances by the median reference sample intensity times the slope
of the regression line and subtracting the result from 1
- for each sample, removes peptides resulting in raw scores higher than the median of all raw scores of a sample
plus three times the median absolute deviation of all raw scores of a sample
- calculates RM score by applying TOP3 normalization to filtered raw scores

Output:
    - FQ-LF-output_raw_scores.csv: contains raw scores for all peptides
    - FQ-LF-output_RM_scores.csv: contains RM scores for filtered peptides plus Slope, R2 model (based on
        RANSAC inliers), R2 data (based on all data points), Reproducibility factor (fraction of the most frequently
        resulting slope of all RANSAC initiations)
    - FQ-LF-output_diff_modified.csv: boolean matrix of same dimention as output_RM_scores.
        True: peptide is differentially modified (RM score smaller than modification cutoff)
        False: peptide is not differentially modified (RM score equal or larger than modification cutoff)
    - FQ-LF-output_removed_peptides.csv: contains all peptides that were removed during raw score filtering

Parameters:
    - input_file: path to .csv file containing peptide intensities (str)
    - output_folder: path of the folder into which the output files should be written (str)
    - reference: identifier of the reference sample(s) (str)
    - num_init: number of ransac initiations (int)
    - mod_cutoff: RM score cutoff used to classify peptides as differentially modified (float between 0 and 1)
    - create_plots: if True, plots are created, if False no plots are created

"""

# Author: Konstantin Kahnert
# Date: 2019_10_15
# Python version: 3.7.3


# # import required libraries
import click
from sys import exit
from tqdm import tqdm
# used to manipulate the data
from pandas import read_csv, Series, DataFrame
from numpy import nan, array, square, sqrt
# used for regular expressions
from re import findall
# used for calculating statistics
from scipy.stats import f
from scipy.stats import median_absolute_deviation
# used for ransac linear regression
from sklearn import linear_model
# used for plotting
from seaborn import distplot, scatterplot, lineplot
from matplotlib import gridspec, use
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
use("pdf")


def add_filename_to_path(path, filename, ending):
    """
    Takes path of input file, removes file extension and adds new filename
    Parameters:
        path: path to input file (string)
        filename: name of new file (string)
        ending: file type e.g. .csv (string)

    Returns:
        inputpath_filename (string)
    """
    # match everything until "." in the path
    filename_start = findall("^(.*?)\..*", filename)[0]

    # add new filename with "_" in between
    new_filename = filename_start + "_" + ending

    # add new_filename to path_output
    path_out = path + "/" + new_filename

    return path_out


def calculate_confidence_band(slope, median_int, dataframe_train, X, y, row, idx, matrix_distance_RL, alpha):
    """ Calculates confidence bands arround the regression line"""

    # calculate predicted intensity with Reference intensity of a peptide and slope of the sample (Y hat)
    Y_pred = slope * median_int

    # calculate W
    N = len(dataframe_train)
    F = f.ppf(q=1 - alpha, dfn=2, dfd=N - 1)
    W = sqrt(2 * F)

    # calculate standard deviation (s(Y hat))
    # calculate prediction error
    error = y - Y_pred
    error.dropna(inplace=True)

    # calculate mean squared error
    MSE = sum(error.apply(square)) / (N - 1)

    # calculate mean X intensity
    X_bar = dataframe_train["Reference intensity"].mean()

    # iterate over all peptides of the sample
    CB_low = []
    CB_high = []

    # iterate through median peptide intensities
    for idx_2, elm in dataframe_train["Reference intensity"].iteritems():
        # calculate squared distance to mean X (numerator)
        dist_X_bar = square(elm - X_bar)

        # calculate sum of squared distances to mean X(denominator)
        sum_dist_X_bar = sum(square(X - X_bar))

        # calculate standard deviation
        s = float(sqrt(MSE * ((1 / N) + (dist_X_bar / sum_dist_X_bar))))

        # calculate predicted intensity for given X
        Y_hat = slope * elm

        # calculate high and low CB values and append to list
        cb_low = Y_hat - W * s
        cb_high = Y_hat + W * s

        CB_low.append(cb_low)
        CB_high.append(cb_high)

    # calculate predicted intensities
    pred_ints = median_int * slope

    # calculate distance to regression line
    distance_RL = pred_ints - row

    # save distances in matrix_distance
    matrix_distance_RL.loc[idx] = distance_RL

    # add CBs as columns to dataframe_train
    dataframe_train["CB low"] = CB_low
    dataframe_train["CB high"] = CB_high

    return matrix_distance_RL, dataframe_train



def create_regression_plots(dataframe_train, idx, r2_score_model, r2_score_data, slope, alpha):
    """ Creates a scatter plot with regression line and confidence bands"""

    # create new figure with two subplots
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 6])
    ax1 = plt.subplot(gs[1])
    ax0 = plt.subplot(gs[0], sharex=ax1)

    # set space between subplots
    gs.update(hspace=0.05)

    # plot histogram in upper subplot
    plt.sca(ax0)

    # add title
    plt.title("RANSAC Linear Regression of Sample " + str(idx))

    # plot histogram
    distplot(a=dataframe_train["Reference intensity"], bins=150, kde=False)

    # remove axis and tick labels
    plt.xlabel("")
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=True,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off

    # plot scatter plot
    plt.sca(ax1)
    scatterplot(x="Reference intensity", y="Sample intensity", data=dataframe_train, hue=dataframe_train.index)

    # draw regression line
    line_label = "R2 model: " + str(r2_score_model) + "\nR2 data: " + str(r2_score_data)
    max_int = dataframe_train["Reference intensity"].max()
    X = [0, max_int]
    y = [0, slope * max_int]
    plt.plot(X, y, color="darkblue", linestyle="-", label=line_label)

    # draw confidence band
    lineplot("Reference intensity", "CB low", data=dataframe_train, color="darkgreen",
             label="CB, alpha=" + str(alpha))
    lineplot("Reference intensity", "CB high", data=dataframe_train, color="darkgreen")

    # set line style of CB lines to dashed
    for i in [1, 2]:
        ax1.lines[i].set_linestyle("--")

    # create legend if sample has 20 peptides or less otherwise don't create a legend
    if len(dataframe_train) <= 20:
        # set right x axis limit
        plt.gca().set_xlim(right=1.4 * max_int)
        plt.legend()
    else:
        plt.gca().get_legend().remove()

    # set y axis label
    plt.ylabel("Intensity sample " + str(idx))
    plt.xlabel("Reference intensity")

    return fig


def calc_raw_scores(df_distance, median_int):
    """
    Parameters:
         df_distance: Pandas DataFrame containing the vertical distances to the regression line
                      rows: samples, columns: peptides
         median_int: Pandas Series containing the median intensities for each peptide of the reference samples
    Returns:
        Pandas DataFrame of same dimension as df_distance containing raw scores
    """
    # copy df_distance
    df_rs = df_distance.copy()

    # iterate through rows of df_distance (samples)
    for idx, row in df_distance.iterrows():
            # extract slope
            slope = row["Slope"]

            # delete slope from row
            row.drop("Slope", inplace=True)

            # calculate raw scores
            raw_scores = 1 - row / (slope * median_int)

            # add slope to raw scores
            raw_scores["Slope"] = slope

            # replace idx row in df_RM_score with calculated raw scores
            df_rs.loc[idx] = raw_scores

    return df_rs


def normalize_t3median(dataframe):
    """
    Applies Top3 median normalization to dataframe
    Determines the median of the three highest values in each row and divides every value in the row by it

    Parameter:
        dataframe: Pandas dataframe of datatype float or integer

    Returns:
        normalized pandas dataframe of same dimensions as input dataframe
    """
    # copy dataframe
    dataframe_t3med = dataframe.copy()

    # for each row, normalize values by dividing each value by the median
    # of the three highest values of the row
    # iterate over rows of dataframe
    for idx, row in dataframe.iterrows():
        # calculate the median of the three highest values
        median_top3 = row.nlargest(3).median()

        # normalize each value of row by dividing by median_top3
        row_norm = row / median_top3

        # update row in dataframe_norm with row_norm
        dataframe_t3med.loc[idx] = row_norm

    return dataframe_t3med


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '-i',
    '--input_file',
    type=click.Path(exists=True),
    required=True,
    help='Path of the comma separated input file (.csv)'
)
@click.option(
    '-o',
    '--output_folder',
    type=click.Path(exists=True),
    required=True,
    help='Path to the folder in which the output files should be placed'
)
@click.option(
    '-r',
    '--reference',
    required=True,
    help='Reference sample identifier. This must match the Group column value of the reference samples.'
)
@click.option(
    '-n',
    '--num_init',
    default=30,
    help='Number of times FLEXIQuant LF fits a new RANSAC linear regression model to each sample to choose the best '
         'model. Should be an integer between 5 and 100. The more initiations, the higher the reproducibility '
         'and the probability that FLEXIQuant LF finds the optimal model. However, choosing a high number can '
         'significantly increases the run time and more than 50 initiations rarely provide additional benefit. '
         'The default value is 30'
)
@click.option(
    '-mc',
    '--mod_cutoff',
    default=0.5,
    help='RM score cutoff used to classify peptides as differentially modified. Should be a a float between 0 and 1. '
         'The default value is 0.5'
)
@click.option('-p',
              '--create_plots',
              is_flag=True,
              help='If selected a linear regression plot for each sample is created. All plots will be saved in '
                   'one pdf file')
def main(input_file, output_folder, reference, num_init, mod_cutoff, create_plots):
    """
    FLEXIQuant LF CLI

    Usage example:

    python FLEXIQuant_LF_CLI.py -i path/to/input_file.csv -o path/to/output_folder/ -r "Control" -n 30 -mc 0.5 -p


    """
    # get input file name
    input_file_name = findall("[^\\\\,/]*$", input_file)[0]

    # load data in pandas data frames
    dataframe_int_matrix = read_csv(input_file, sep=",", index_col=0)

    # check if Group column exists
    try:
        tmp = dataframe_int_matrix["Group"]

    except KeyError:
        # print error message
        print("Error: Incorrect input format!\n" + "No column named \"Group\" found in input file.")

        # terminate process
        exit()

    # check given reference identifier exists in Group column
    if str(reference) not in set(dataframe_int_matrix["Group"].astype(str)):

        # print error message
        print("Error: Given reference sample identifier \"" + str(reference) + "\" not found in \"Group\" column")

        # terminate process
        exit()

    else:
        # filter dataframe for controls
        dataframe_int_matrix_control = dataframe_int_matrix[dataframe_int_matrix["Group"].astype(str) ==
                                                            str(reference)]

    # get modified peptides
    modified = []
    for elm in dataframe_int_matrix_control.columns:
        if "ph" in elm or "ac" in elm or "gly" in elm:
            modified.append(elm)

    # delete modified peptides
    dataframe_int_matrix_control = dataframe_int_matrix_control.drop(modified, axis=1)
    dataframe_int_matrix.drop(modified, inplace=True, axis=1)

    # delete group column
    group_column = dataframe_int_matrix["Group"]
    dataframe_int_matrix_control.drop("Group", inplace=True, axis=1)
    dataframe_int_matrix.drop("Group", inplace=True, axis=1)

    # calculate median intensities for unmodified peptides of controls
    median_int = dataframe_int_matrix_control.median(axis=0)

    # delete columns where all entries are nan
    dataframe_int_matrix.dropna(how="all", axis=1)

    # initiate empty lists to save results of linear regressions
    slope_list = []
    r2_score_list_model = []
    r2_score_list_data = []
    reproducibility_list = []

    # copy dataframe_int_matrix
    matrix_distance_RL = dataframe_int_matrix.copy()

    # check if plots should be created or not
    if create_plots:
        # create output file path
        path_plots = add_filename_to_path(output_folder, input_file_name, "FQ-LF_regression_plots.pdf")

        # try to excess file and raise permission error if not possible
        try:
            plots_pdf = PdfPages(path_plots)
        except PermissionError:
            # print error message
            print("Error: Permission denied!\n" + "Please close " + path_plots)

            # terminate process
            exit()
    else:
        pass

    # determine progress step size and set progress to 0
    num_samples = len(dataframe_int_matrix)
    progress_step = (100 / num_samples) * 0.99
    progress = 0

    # iterate through rows (samples) of dataframe_int_matrix
    for idx, row in tqdm(dataframe_int_matrix.iterrows(), total=dataframe_int_matrix.shape[0]):

        # create data frame with sample intensities and median intensities
        dataframe_train = DataFrame({"Sample intensity": row, "Reference intensity": median_int})

        # remove nan values
        dataframe_train.dropna(inplace=True)

        # sort dataframe
        dataframe_train.sort_index(inplace=True, axis=0)

        # if number of peptides is smaller than 5, skip sample and continue with next interation
        if len(dataframe_train) < 5:

            # set all metrics to nan
            matrix_distance_RL.loc[idx] = nan
            dataframe_train["CB low"] = nan
            dataframe_train["CB high"] = nan
            slope_list.append(nan)
            r2_score_list_model.append(nan)
            r2_score_list_data.append(nan)
            reproducibility_list.append(nan)

            # continue with next iteration
            continue

        else:
            pass

        # initiate empty list to save results of linear regression
        list_model_slopes = []
        list_r2_model = []
        list_r2_data = []
        all_iter_slopes = []
        all_iter_r2_data = []
        all_iter_r2_model = []


        # set training data
        X = array(dataframe_train["Reference intensity"]).reshape(-1, 1)
        y = dataframe_train["Sample intensity"]

        # calculate squared MAD
        sq_mad = square(dataframe_train["Sample intensity"].mad())

        # run ransac linear regression num_ransac_init times to select best fitting model
        for i in range(num_init):
            # initiate linear regression model with ransac regressor
            ransac_model = linear_model.RANSACRegressor(base_estimator=
                                                        linear_model.LinearRegression(fit_intercept=False, n_jobs=-2),
                                                        max_trials=1000,
                                                        stop_probability=1,
                                                        min_samples=0.5,
                                                        loss="squared_loss",
                                                        residual_threshold=sq_mad)
            # fit model
            ransac_model.fit(X, y)

            # get coefficient
            slope = (float(ransac_model.estimator_.coef_))

            # get inlier and outlier
            inlier_mask = ransac_model.inlier_mask_

            # add as column to dataframe_train
            dataframe_train["Outlier"] = ~inlier_mask.astype(bool)

            # calculate r2 score based on inliers
            dataframe_train_inlier = dataframe_train[dataframe_train["Outlier"] == False]
            X_inlier = array(dataframe_train_inlier["Reference intensity"]).reshape(-1, 1)
            y_inlier = dataframe_train_inlier["Sample intensity"]
            r2_score_model = round(ransac_model.score(X_inlier, y_inlier), 4)
            r2_score_data = round(ransac_model.score(X, y), 4)

            # save model and r2 score to corresponding lists
            list_model_slopes.append(slope)
            list_r2_model.append(r2_score_model)
            list_r2_data.append(r2_score_data)

        # save all slopes and score to corresponding lists
        all_iter_slopes.append(list_model_slopes)
        all_iter_r2_model.append(list_r2_model)
        all_iter_r2_data.append(list_r2_data)

        # determine best model based on r2 scores
        best_model = list_r2_model.index(max(list_r2_model))

        # save slope of best model to slope_list
        slope_list.append(list_model_slopes[best_model])
        slope = list_model_slopes[best_model]

        # calculate reproducibility factor and save to list
        series_slopes = Series(list_model_slopes)
        reproducibility_factor = max(series_slopes.value_counts()) / num_init
        reproducibility_list.append(reproducibility_factor)

        # get r2 scores of best model
        r2_score_model = list_r2_model[best_model]
        r2_score_data = list_r2_data[best_model]

        # save best r2 score to lists
        r2_score_list_model.append(r2_score_model)
        r2_score_list_data.append(r2_score_data)

        # check if plots need to be created
        if create_plots:
            # calculate confidence band
            alpha = 0.3
            matrix_distance_RL, dataframe_train = calculate_confidence_band(slope, median_int, dataframe_train, X, y,
                                                                            row, idx, matrix_distance_RL, alpha)

            # plot scatter plot with regression line
            fig = create_regression_plots(dataframe_train, idx, r2_score_model, r2_score_data, slope, alpha)

            # save plot
            plots_pdf.savefig(figure=fig, bbox_inches='tight')
            plt.close()

        else:
            # calculate predicted intensities
            pred_ints = median_int * slope

            # calculate distance to regression line
            distance_RL = pred_ints - row

            # save distances in matrix_distance
            matrix_distance_RL.loc[idx] = distance_RL


    # close pdf file if plots were created
    if create_plots:
        # close pdf file
        plots_pdf.close()
    else:
        pass

    # add slope to dataframe
    matrix_distance_RL["Slope"] = slope_list

    # calculate raw scores
    dataframe_raw_scores = calc_raw_scores(matrix_distance_RL, median_int)

    # calculate MAD per sample
    dataframe_raw_scores.drop("Slope", axis=1, inplace=True)
    dataframe_raw_scores_T = dataframe_raw_scores.T
    mad = median_absolute_deviation(dataframe_raw_scores_T, scale=1)
    median = dataframe_raw_scores_T.median(axis=0)

    # calculate cutoff value for each time point (> 3*MAD)
    cutoff = median + 3 * mad

    # remove peptides with raw scores > cutoff for each sample
    dataframe_raw_scores_T_cutoff = dataframe_raw_scores_T[round(dataframe_raw_scores_T, 5) <= round(cutoff, 5)]
    removed = Series(dataframe_raw_scores_T_cutoff.index[dataframe_raw_scores_T_cutoff.isna().any(axis=1)])
    dataframe_raw_scores_T_cutoff.dropna(axis=0, how="any", inplace=True)
    dataframe_raw_scores_cutoff = dataframe_raw_scores_T_cutoff.T

    # apply t3median normalization to calculate RM scores
    dataframe_RM = normalize_t3median(dataframe_raw_scores_cutoff)

    # check if peptides are modified (RM score below modification cutoff)
    dataframe_RM_mod = dataframe_RM < mod_cutoff

    # add metrics of regline to matrix_distance_RL
    dataframe_raw_scores["Slope"] = slope_list
    dataframe_raw_scores["R2 model"] = r2_score_list_model
    dataframe_raw_scores["R2 data"] = r2_score_list_data
    dataframe_raw_scores["Reproducibility factor"] = reproducibility_list

    dataframe_RM["Slope"] = slope_list
    dataframe_RM["R2 model"] = r2_score_list_model
    dataframe_RM["R2 data"] = r2_score_list_data
    dataframe_RM["Reproducibility factor"] = reproducibility_list

    # add Group column again
    dataframe_raw_scores["Group"] = group_column
    dataframe_RM["Group"] = group_column
    dataframe_RM_mod["Group"] = group_column

    # save raw scores as csv file, raise permission error if file can't be accessed
    path_out = add_filename_to_path(output_folder, input_file_name, "FQ-LF-output_raw_scores.csv")
    try:
        dataframe_raw_scores.to_csv(path_out, sep=",", index=True)
    except PermissionError:
        # print error message
        print("Error: Permission denied!\n" + "Please close " + path_out)

        # terminate process
        exit()

    # save RM scores as csv file, raise permission error if file can't be accessed
    path_out = add_filename_to_path(output_folder, input_file_name, "FQ-LF-output_RM_scores.csv")
    try:
        dataframe_RM.to_csv(path_out, sep=",", index=True)
    except PermissionError:
        # print error message
        print("Error: Permission denied!\n" + "Please close " + path_out)

        # terminate process
        exit()

    # save differentially modified dataframe as csv file, raise permission error if file can't be accessed
    path_out = add_filename_to_path(output_folder, input_file_name, "FQ-LF-output_diff_modified.csv")
    try:
        dataframe_RM_mod.to_csv(path_out, sep=",", index=True)
    except PermissionError:
        # print error message
        print("Error: Permission denied!\n" + "Please close " + path_out)

        # terminate process
        exit()

    # save removed peptides as csv file, raise permission error if file can't be accessed
    path_out = add_filename_to_path(output_folder, input_file_name, "FQ-LF-output_removed_peptides.csv")
    try:
        removed.to_csv(path_out, sep=",", index=False, header=False)
    except PermissionError:
        # print error message
        print("Error: Permission denied!\n" + "Please close " + path_out)

        # terminate process
        exit()


if __name__ == '__main__':
    main()
