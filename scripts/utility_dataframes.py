import numpy as np
import pandas as pd
import sys, fileinput

# These options will format floating-point values in scientific notation as specified below and will display the 
# complete contents of Pandas DataFrames without any kind of truncation. 
pd.set_option("display.float_format","{:+.8e}".format)
pd.set_option("display.max_rows",None)
pd.set_option("display.max_columns",None)
pd.set_option("display.width",None)
pd.set_option("display.max_colwidth",None)

def write_dataframe_to_fwf(file_name, df, keep_index_column=False, width_index_column=None):
    """ 
    Writes the contents of a Pandas DataFrame to an output file in fixed-width format (fwf), omitting the index column if specified to do so.

    Parameters:
        1) file_name (string): Directory location and name of the output file.
        2) df (Pandas DataFrame): DataFrame in which the data to be written to the output file are stored.
        3) keep_index_column (bool): Specifies whether to print the index column in the output file (default = False).
        4) width_index_column (int): Specifies the number of characters that span the width of the index column (default = None).

    Returns:
        N/A
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

def write_data_and_labels_to_fwf(file_name, data, column_labels, transpose_data=False, 
                                 keep_index_column=False, print_to_console=False):
    """ 
    Stores the given data and column labels in a Pandas DataFrame and writes the contents of the DataFrame to an output file in fixed-width format.

    Parameters:
        1) file_name (string): Directory location and name of the output file.
        2) data (Array or list): Contains one or more columns of data.
        3) column_labels (list): List of strings representing the column labels.
        4) transpose_data (bool): Indicates if it is necessary to transpose the data (which it would be if the data are stored in a list of lists,
            for example) so that the data appear along the columns of the output file rather along the rows (default = False).
        5) keep_index_column (bool): Specifies whether to print the index column in the output file (default = False).
        6) print_to_console (bool): Indicates if the contents of the data should be printed to the console (default = False).

    Returns:
        N/A
    """
    df = pd.DataFrame(data)
    if transpose_data:
        df = df.T
    df.columns = column_labels
    if print_to_console:
        print(df)
    write_dataframe_to_fwf(file_name, df, keep_index_column)

def write_data_and_labels_to_csv(file_name, data, column_labels=None, transpose_data=False, 
                                 keep_index_column=False, keep_header=True, separation='\t', format='%12.8e'):
    """ 
    Stores the given data and column labels in a Pandas DataFrame and writes the contents of the DataFrame to a csv file.

    Parameters:
        1) file_name (string): Directory location and name of the output file.
        2) data (Array or list): Contains one or more columns of data.
        3) column_labels (list): List of strings representing the column labels.
        4) transpose_data (bool): Indicates if it is necessary to transpose the data (which it would be if the data are stored in a list of lists,
            for example) so that the data appear along the columns of the output file rather along the rows (default = False).
        5) keep_index_column (bool): Specifies whether to print the index column in the output file (default = False).
        6) keep_header (bool): Indicates if the output file should contain headers for column labels (default = True).
        7) separation (string): Indicates the character(s) to be used as column separators (default = tabs).
        8) format (string): Format for floats in the output file.

    Returns:
        N/A
    """
    df = pd.DataFrame(data)
    if transpose_data:
        df = df.T
    if column_labels:
        df.columns = column_labels
    df.to_csv(file_name, index=keep_index_column, header=keep_header, sep=separation, float_format=format)

def clean_up_dataframe(df, print_to_console=False):
    """ 
    Cleans up a Pandas DataFrame from junk columns that contain only NaN values and returns this cleaned-up DataFrame.
    The junk columns may be a result of pd.read_fwf() separating out a single column with a long label into two columns where the second column 
    contains the junk values. For example, the column "rho A+B+C (g/cm^3)" may get separated to "rho A+B+C" and "(g/cm^3)", in which the latter
    column, which will contain only NaN entries, should be deleted.

    Parameters:
        1) df (Pandas DataFrame): DataFrame to be cleaned.
        2) print_to_console (bool): Indicates if the contents of the data should be printed to the console to monitor the progress (default = False).

    Returns:
        df (Pandas DataFrame): DataFrame after cleaning operations have been applied.
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