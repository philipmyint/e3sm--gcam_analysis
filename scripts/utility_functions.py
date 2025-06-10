import numpy as np
import os
import re

def add_lists_elementwise(list1, list2, list2_are_units=False):
    """
    Performs elementwise addition of two lists.

    Parameters:
        list1: The first list.
        list2: The second list.
        list2_are_units: Boolean that specifies whether list2 represents the units corresponding to list1. This would be useful when forming 
                         column headers, where list1 specifies the quantities and list2 specifies the corresponding units for those quantities.

    Returns:
        A new list containing the element-wise sums of list1 and list2. If list2_are_units is True, then the elements of list2 will be in parentheses 
        so that each element of the combined list will be of the form 'a (b)', where a is from list1 and b is from list2.
        Returns None if the lists are not of the same length.
    """
    if len(list1) != len(list2):
        return None
    if list2_are_units:
        return [a + ' (' + b + ')' for a, b in zip(list1, list2)]
    else:
        return [a + b for a, b in zip(list1, list2)]

def check_is_list_of_lists(data):
    """
    Checks if the given iterable is a list of lists.

    Parameters:
        data: Iterable we want to check is a list of lists.

    Returns:
        True if data is a list of lists, False otherwise.
    """
    if not isinstance(data, list):
        return False
    return all(isinstance(item, list) for item in data)

def check_substrings_in_list(substrings, list, all_or_any='all'):
    """
    Checks if either all or any of the elements of the substrings list are substrings of at least one element in list.

    Parameters:
        substring: A list of strings (substrings to search for).
        list: A list of strings (strings to search within).
        all_or_any: String whose value should be either 'all' or 'any'.

    Returns:
        True if either all or any of the elements of substrings are substrings of at least one element in list, False otherwise.
    """
    if all_or_any == 'all':
        return all(any(substring in element for element in list) for substring in substrings)
    else:
        return any(any(substring in element for element in list) for substring in substrings)

def check_substrings_in_string(substrings, string, all_or_any='all'):
    """
    Checks if either all or any of the elements of the substrings list are substrings of the string.

    Parameters:
        substring: A list of strings (substrings to search for).
        string: String to search within.
        all_or_any: String whose value should be either 'all' or 'any'.

    Returns:
        True if either all or any of the elements of the substrings list are substrings of the string, False otherwise.
    """
    if all_or_any == 'all':
        return all(substring in string for substring in substrings)
    else:
        return any(substring in string for substring in substrings)

def convert_month_numbers_to_(substrings, string, all_or_any='all'):
    """
    Checks if either all or any of the elements of the substrings list are substrings of the string.

    Parameters:
        substring: A list of strings (substrings to search for).
        string: String to search within.
        all_or_any: String whose value should be either 'all' or 'any'.

    Returns:
        True if either all or any of the elements of substrings are substrings of the string, False otherwise.
    """
    if all_or_any == 'all':
        return all(substring in string for substring in substrings)
    else:
        return any(substring in string for substring in substrings)

def create_numpy_array_from_ds(ds, variables, fill_nan_values):
    """
    Creates a list of NumPy arrays from the specified variables of an xarray Dataset. 

    Parameters:
        ds: xarray Dataset.
        variables: List of variables.
        fill_nan_values: List that indicates what to set NaN values to for each variable.

    Returns:
        List of arrays, one for each ds variable. If only one variable is specified, then a single array (not a list with this array) is returned.
    """
    np_arrays = []
    for index, column in enumerate(variables):
        array = ds[column].to_numpy().reshape(-1,1)
        array = np.nan_to_num(array, nan=fill_nan_values[index])
        np_arrays.append(array)
    if len(np_arrays) == 1:
        return array
    else:
        return np_arrays

def get_all_files_in_path(path, file_name_substrings=None, file_extension=None):
    """
    Gets a list of complete paths for all files that are in a particular directory.

    Parameters:
        path: Path of directory where we want to search for files.
        file_name_substrings: A list of all substrings that must be in the file names.
        file_extension: File extension that should be in all files. 

    Returns:
        A list with complete paths to all files in the directory. 
        If both file_name_substrings and file_extension are None, then all files in the directory will be included in the list.
    """
    file_paths = []
    for root, _, files in os.walk(path):
        for file in files:
            if not file_name_substrings or all([substring in file for substring in file_name_substrings]):
                if not file_extension or file.endswith(file_extension):
                    file_path = os.path.join(root, file)
                    file_paths.append(file_path)
    return file_paths

def modify_list_based_on_condition(original_list, condition, new_value_function):
    """
    Modifies a list by applying a condition and a function to generate new values.

    Parameters:
        original_list: The list to be modified.
        condition: A function that takes an element as input and returns True if the condition is met, and False otherwise.
        new_value_function: A function that takes an element as input and returns the new value to replace the original element.

    Returns:
        A new list with the modified elements.
    """
    return [new_value_function(element) if condition(element) else element for element in original_list]

def print_p_values(ttest, variable, p_value_threshold, p_value_file, output_file_or_label, p_value_file_print_only_if_below_threshold):
    """
    Prints the p-values from a t-test to the console and optionally prints to an output file.

    Parameters:
        ttest: t-test object.
        variable: Variable for which the t-test was performed.
        p_value_threshold: Threshold for the p-value. The message to the console will indicate if the p-value falls below this threshold.
        p_value_file: Path and name for the file where the p-value result will be printed.
        output_file_or_label: Output file or label from where the data used to perform the t-test were obtained.
        p_value_file_print_only_if_below_threshold: If true, the p-value gets printed to the file only if it falls below the threshold.

    Returns:
        N/A.
    """
    if ttest.pvalue < p_value_threshold:
        print(f'p-value of {variable} in {output_file_or_label}: {ttest.pvalue:.4e}, which is less than {p_value_threshold}')
        if p_value_file:
            with open(p_value_file, 'a+') as f:
                f.write(f'{variable} in {output_file_or_label}: {ttest.pvalue:.4e}\n')
    else:
        print(f'p-value of {variable} in {output_file_or_label}: {ttest.pvalue:.4e}')
        if p_value_file and not p_value_file_print_only_if_below_threshold:
            with open(p_value_file, 'a+') as f:
                f.write(f'{variable} in {output_file_or_label}: {ttest.pvalue:.4e})\n')

def replace_inside_parentheses(text, replacement):
    """
    Replaces the string inside parentheses with a given replacement string.

    Parameters:
        text: The input string.
        replacement: The string with which we will replace the content inside parentheses.

    Returns:
        The modified string with the content inside parentheses replaced.
    """
    return re.sub(r"\([^)]*\)", replacement, text)

def sort_file(file_path):
    """
    Sorts the lines of a file alphabetically.

    Parameters:
        file_path: The file to be sorted.

    Returns:
        N/A.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        return f"Error: File not found: {file_path}"
    lines.sort()
    with open(file_path, 'w') as file:
        file.writelines(lines)
    return f"File sorted successfully: {file_path}"