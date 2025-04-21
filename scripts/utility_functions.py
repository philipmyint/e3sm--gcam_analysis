import numpy as np
import os

def add_lists_elementwise(list1, list2, list2_are_units=False):
    """
    Performs elementwise addition of two lists.

    Parameters:
        1) list1: The first list.
        2) list2: The second list.
        3) list2_are_units: Boolean that specifies whether list2 represents the units corresponding to list1,
        the latter of which are presumably column headers.

    Returns:
        A new list containing the element-wise sums of list1 and list2.
        Returns None if the lists are not of the same length.
    """
    if len(list1) != len(list2):
        return None
    if list2_are_units:
        return [a + ' (' + b + ')' for a, b in zip(list1, list2)]
    else:
        return [a + b for a, b in zip(list1, list2)]
    
def check_substrings_in_list(substrings, list, all_or_any='all'):
    """
    Checks if either all or any of the elements of the substrings list are substrings of at least one element in list.

    Parameters:
        substring: A list of strings (substrings to search for).
        list: A list of strings (strings to search within).

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

    Returns:
        True if either all or any of the elements of substrings are substrings of the string, False otherwise.
    """
    if all_or_any == 'all':
        return all(substring in string for substring in substrings)
    else:
        return any(substring in string for substring in substrings)

def create_numpy_array_from_xrds_columns(xrds, columns, fill_nan_values):
    np_arrays = []
    for index, column in enumerate(columns):
        array = xrds[column].to_numpy().reshape(-1,1)
        array = np.nan_to_num(array, nan=fill_nan_values[index])
        np_arrays.append(array)
    if len(np_arrays) == 1:
        return array
    else:
        return np_arrays

def get_all_files_in_path(path, file_name_substrings=None, file_extension=None):
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
        condition: A function that takes an element as input and returns True if the
                    condition is met, and False otherwise.
        new_value_function: A function that takes an element as input and returns
                            the new value to replace the original element.

    Returns:
        A new list with the modified elements.
    """
    return [new_value_function(element) if condition(element) else element for element in original_list]