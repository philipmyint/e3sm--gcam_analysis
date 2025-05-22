import fileinput
import numpy as np
import pandas as pd
from scipy import stats
import sys

# These options will format floating-point values in scientific notation as specified below and will display the 
# complete contents of Pandas DataFrames without any kind of truncation. 
pd.set_option("display.float_format","{:+.8e}".format)
pd.set_option("display.max_rows",None)
pd.set_option("display.max_columns",None)
pd.set_option("display.width",None)
pd.set_option("display.max_colwidth",None)

def clean_up_dataframe(df, print_to_console=False):
    """ 
    Cleans up a Pandas DataFrame from junk columns that contain only NaN values and returns this cleaned-up DataFrame.
    The junk columns may be a result of pd.read_fwf() separating out a single column with a long label into two columns where the second column 
    contains the junk values. For example, the column "rho A+B+C (g/cm^3)" may get separated to "rho A+B+C" and "(g/cm^3)", in which the latter
    column, which will contain only NaN entries, should be deleted.

    Parameters:
        df: DataFrame to be cleaned.
        print_to_console: Indicates if the contents of the data should be printed to the console to monitor the progress.

    Returns:
        df: DataFrame after cleaning operations have been applied.
    """
    if print_to_console:
        print(f'\nThe original columns of the DataFrame are:\n{df.columns}\n')

    # Number of columns.
    num_columns = len(df.columns)

    # List containing the labels of junk columns.
    columns_to_delete = []

    # List containing the labels of non-junk columns.
    columns_to_keep = []

    # Loop over the columns in the DataFrame.
    for i in np.arange(num_columns):

        # Extract the label of the current column.
        label = df.columns[i]

        # If this is a junk column, record its label and go to the next column.
        if df[label].isnull().values.all():
            columns_to_delete.append(label)
            if print_to_console:
                print(f'Column {i} contains junk; its label is {label}\n')
            continue

        # If the next column is a junk column and we are not in the last column, append the label of that 
        # next column (with some modifications) to the label of the current one.
        if i < num_columns - 1 and df[df.columns[i+1]].isnull().values.all():
            # If the label of the junk column indicates the units of the quantity in the previous column 
            # placed between parentheses, extract only the parentheses and the unit label between them.
            start = df.columns[i+1].find("(")
            if start != -1:
                end = df.columns[i+1].find(")") + 1
                label += " " + df.columns[i+1][start:end]
            # If the label of the junk column does not contain parentheses, then append that label to that of the current column.
            else:
                label += " " + df.columns[i+1]
            if print_to_console:
                print(f'The new label of column {i} is {label}')

        # Record the labels of all non-junk columns in the appropriate list.
        columns_to_keep.append(label)

    # Delete all junk columns.
    for label in columns_to_delete:
        df.drop(label, axis=1, inplace=True)
    
    # Reassign the column labels to those that we want to keep; print these labels and return the DataFrame.
    df.columns = columns_to_keep
    if print_to_console:
        print(f'The new columns of the modified DataFrame are:\n{df.columns}\n')
    return df

def get_columns_without_units_in_dataframe(df):
    """
    Returns all column names in a Pandas DataFrame without the units. For example, return 'Mass' if the colum name is 'Mass (kg)'.

    Parameters:
        df: The input DataFrame.

    Returns:
        A list of all column names, with the units removed.
    """
    columns_without_units = []
    for column in df.columns:
        stop_index = column.find(' (')
        # If the column does not have unit associated with it (so that find returns a -1), just return that column unmodified.
        if stop_index == -1:
            columns_without_units.append(column)
        else:
            columns_without_units.append(column[:stop_index])
    return columns_without_units

def get_matching_column_in_dataframe(df, variable, all_matches=False):
    """
    Returns the column name(s) in a Pandas DataFrame that matches the variable.

    Parameters:
        df: The input DataFrame.
        variable: Variable name to search for in the column names.
        all_matches: If True, all columns that matches the variable will be returned. If False, only the first matching column will be returned.

    Returns:
        The matching column name(s), or None if no match is found.
    """
    all_matching_columns = []
    for column in df.columns:
        if variable + ' (' in column:
            if not all_matches:
                return column
            else:
                all_matching_columns.append(column)
    if all_matching_columns:
        return all_matching_columns
    else:
        return None

def move_columns_next_to_each_other_in_dataframe(df, column1, column2):
    """
    Rearranges a Pandas DataFrame by moving column2 next to column1.

    Parameters:
        column1: Column next to which column2 will be moved. The order of all columns prior to column1 will be unaffected.
        column2: Column to move next to column1.
        df: DataFrame whose columns will be rearranged.

    Returns:
        A new DataFrame with column1 and column2 moved next to each other.
    """
    index = df.columns.get_loc(column1)
    new_columns = list(df.columns[:index]) + [column1, column2] + list(df.columns[index+1:])
    df = df[new_columns]
    # Remove the duplicate column2, keeping the first occurrence.
    df = df.loc[:, ~df.columns.duplicated()]
    return df

def perform_ttest(df, columns_set_1, columns_set_2):
    """ 
    Performs a t-test for the means of two (presumed) independent data sets.

    Parameters:
        df: Pandas DataFrame containing the columns for both data sets.
        columns_set_1: List of columns in the DataFrame for the first data set.
        columsn_set_2: List of columns in the DataFrame for the second data set.
        
    Returns:
        The p-value produced by the t-test.
    """
    set_1 = df[columns_set_1]
    set_2 = df[columns_set_2]
    ttest = stats.ttest_ind(set_1, set_2)
    return ttest.pvalue

def read_file_into_dataframe(file_name, clean_up_df=False):
    """ 
    Reads a csv or fixed-width-format file, puts the contents into a Pandas DataFrame, and returns the DataFrame.

    Parameters:
        file_name: Complete path and name of the output file.
        clean_up_df: Boolean that specifies if we want to call clean_up_dataframe() on the DataFrame before returning it.

    Returns:
        DataFrame containing the contents of the file.
    """
    if file_name.endswith('.csv'):
        df = pd.read_csv(file_name)
    else:
        df = pd.read_fwf(file_name)
    if clean_up_df:
        df = clean_up_dataframe(df)
    return df

def write_data_and_labels_to_csv(file_name, data, column_labels=None, transpose_data=False, 
                                 keep_index_column=False, keep_header=True, separation=',', format='%12.8e'):
    """ 
    Stores the given data and column labels in a Pandas DataFrame and writes the contents of the DataFrame to a csv file.

    Parameters:
        file_name: Complete path and name of the output file.
        data: Contains one or more columns of data.
        column_labels: List of strings representing the column labels.
        transpose_data: Indicates if it is necessary to transpose the data (which it would be if the data are stored in a list of lists, for example)
                        so that the data appear along the columns of the output file rather along the rows.
        keep_index_column: Specifies whether to print the index column in the output file.
        keep_header: Indicates if the output file should contain headers for column labels.
        separation: Indicates the character(s) to be used as column separators.
        format: Format for floats in the output file.

    Returns:
        N/A.
    """
    df = pd.DataFrame(data)
    if transpose_data:
        df = df.T
    if column_labels:
        df.columns = column_labels
    df.to_csv(file_name, index=keep_index_column, header=keep_header, sep=separation, float_format=format)

def write_data_and_labels_to_fwf(file_name, data, column_labels, transpose_data=False, keep_index_column=False, print_to_console=False):
    """ 
    Stores the given data and column labels in a Pandas DataFrame and writes the contents of the DataFrame to an output file in fixed-width format.

    Parameters:
        file_name: Complete path and name of the output file.
        data: Contains one or more columns of data.
        column_labels: List of strings representing the column labels.
        transpose_data: Indicates if it is necessary to transpose the data (which it would be if the data are stored in a list of lists, for example)
                        so that the data appear along the columns of the output file rather along the rows.
        keep_index_column: Specifies whether to print the index column in the output file.
        print_to_console: Indicates if the contents of the data should be printed to the console.

    Returns:
        N/A.
    """
    df = pd.DataFrame(data)
    if transpose_data:
        df = df.T
    df.columns = column_labels
    if print_to_console:
        print(df)
    write_dataframe_to_fwf(file_name, df, keep_index_column)

def write_dataframe_to_fwf(file_name, df, keep_index_column=False, width_index_column=None):
    """ 
    Writes the contents of a Pandas DataFrame to an output file in fixed-width format (fwf), omitting the index column if specified to do so.

    Parameters:
        file_name: Complete path and name of the output file.
        df: DataFrame in which the data to be written to the output file are stored.
        keep_index_column: Specifies whether to print the index column in the output file.
        width_index_column: Specifies the number of characters that span the width of the index column.

    Returns:
        N/A.
    """
    # Write the contents to the file.
    with open(file_name, 'w') as file:
        file.write(df.__repr__())
    
    # Remove the index column by deleting the number of characters equal to the width of the column.
    # This includes the space between it and the next column.
    if (keep_index_column == False):

        # If the width of the index column has not been specified in the call to this function, set it based on the number of lines in the file.
        if not width_index_column:
            num_lines = len(df)
            if num_lines > 1:
                width_index_column = int(np.log10(num_lines-1) + 2)
            else:
                width_index_column = 2

        for line in fileinput.input(files=(file_name), inplace=True):
            sys.stdout.write(line[width_index_column:])