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
    
def find_between_chars(text, start_char, end_char):
    """
    Finds and returns the substring located between the first occurrence of start_char and the first occurrence of end_char after start_char.

    Parameters:
        text: The string to search within.
        start_char: The character marking the beginning of the desired substring.
        end_char: The character marking the end of the desired substring.

    Returns:
        str or None: The extracted substring if both characters are found in order, otherwise None.
    """
    start_index = text.find(start_char)
    if start_index == -1:
        return None  # Start character not found
    end_index = text.find(end_char, start_index + 1)
    if end_index == -1:
        return None  # End character not found after start character
    return text[start_index + 1:end_index]

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
                f.write(f'{variable} in {output_file_or_label}: {ttest.pvalue:.4e}\n')

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

def transpose_scenarios_if_needed(scenarios, scenario_sets=None):
    """
    Detects the format of the scenarios list of lists representing collections of ensembles of scenarios, 
    and transposes the list of lists if necessary to match the internal format expected by the plotting functions.
    
    The plotting functions expect scenarios organized by ensemble member (rows = members, columns = scenario sets):
        [["Control", "Full feedback"],           # Member 1
         ["Control_2", "Full feedback_2"],       # Member 2
         ["Control_3", "Full feedback_3"]]       # Member 3
    
    This function also supports an alternative user-friendly format organized by scenario set (rows = sets, columns = members):
        [["Control", "Control_2", "Control_3"],           # Control set
         ["Full feedback", "Full feedback_2", "Full feedback_3"]]  # Full feedback set
    
    The function uses the following heuristics to detect the format:
    1. If scenario_sets is provided, its length should match the number of scenario sets (columns in the original format,
       rows in the alternative format). This is the primary detection method and works regardless of scenario naming conventions.
    2. If scenario_sets is not provided and the matrix is not square, the function assumes the format where rows > columns
       is the original format (more ensemble members than scenario sets is typical).
    3. For square matrices without scenario_sets, the function cannot determine the format and returns as-is.
    
    Parameters:
        scenarios: A list of lists containing scenario names. Can be in either format. Scenario names can follow any 
                   naming convention - they do not need to share common base names or follow patterns like "Name_1", "Name_2", etc.
        scenario_sets: Optional list of scenario set names (e.g., ["Control", "Full feedback"]). This is the recommended
                       way to specify the format, as it allows unambiguous detection regardless of scenario naming conventions.
                       When provided, the function compares len(scenario_sets) to the dimensions of the scenarios matrix
                       to determine which format is being used.
    
    Returns:
        A tuple of (transposed_scenarios, was_transposed):
        - transposed_scenarios: The scenarios in the internal format (organized by ensemble member).
        - was_transposed: Boolean indicating whether the input was transposed (True if the input was in the 
          alternative format and was converted to the internal format).
    
    Examples:
        # Example 1: Using scenario_sets for unambiguous detection (recommended)
        scenarios = [["run_A", "run_B", "run_C"], ["exp_X", "exp_Y", "exp_Z"]]
        scenario_sets = ["Control", "Treatment"]
        # len(scenario_sets) = 2 matches num_rows = 2, so this is the alternative format
        # Result: [["run_A", "exp_X"], ["run_B", "exp_Y"], ["run_C", "exp_Z"]], True
        
        # Example 2: Without scenario_sets, using row/column count heuristic
        scenarios = [["A1", "B1"], ["A2", "B2"], ["A3", "B3"]]  # 3 rows, 2 cols
        # More rows than columns suggests original format (3 members, 2 sets)
        # Result: scenarios (unchanged), False
        
        # Example 3: Square matrix without scenario_sets (ambiguous)
        scenarios = [["A", "B"], ["C", "D"]]  # 2x2 - could be either format
        # Result: scenarios (unchanged), False
    """
    if not scenarios or not isinstance(scenarios[0], list):
        # Not a list of lists; return as-is.
        return scenarios, False
    
    num_rows = len(scenarios)
    num_cols = len(scenarios[0]) if scenarios else 0
    
    # If all inner lists don't have the same length, we cannot reliably detect format
    # abort with message
    if not all(len(row) == num_cols for row in scenarios):
        print("In json input file: netdcf_files requires equal length lists of file names")
        os._exit(1)

    # Heuristic 1 (Primary): If scenario_sets is provided, use its length to determine the format.
    # This is the most reliable method and works regardless of scenario naming conventions.
    if scenario_sets is not None:
        num_sets = len(scenario_sets)
        if num_sets == num_rows and num_sets != num_cols:
            # scenario_sets length matches rows, so input is in alternative format (organized by scenario set).
            # Transpose to convert to internal format (organized by ensemble member).
            transposed = [[scenarios[row][col] for row in range(num_rows)] for col in range(num_cols)]
            return transposed, True
        elif num_sets == num_cols and num_sets != num_rows:
            # scenario_sets length matches columns, so input is already in internal format.
            return scenarios, False
        elif num_sets == num_rows == num_cols:
            # Square matrix and scenario_sets matches both dimensions - ambiguous.
            # In this case, we assume the user is using the new format (organized by scenario set) since that's
            # the more intuitive format for users. Each row in scenarios corresponds to a scenario set.
            transposed = [[scenarios[row][col] for row in range(num_rows)] for col in range(num_cols)]
            return transposed, True
        # If num_sets matches neither, fall through to heuristic 2.
    
    # Heuristic 2: Without scenario_sets, use the shape of the matrix.
    # Typically, there are more ensemble members than scenario sets (e.g., 5 members x 2 sets).
    # So if rows > cols, assume original format; if cols > rows, assume alternative format and transpose.
    if num_cols > num_rows:
        # More columns than rows suggests alternative format (few sets, many members per set).
        # Transpose to convert to internal format.
        transposed = [[scenarios[row][col] for row in range(num_rows)] for col in range(num_cols)]
        return transposed, True
    else:
        # rows >= cols: assume original format or ambiguous square matrix
        # abort with error message
        print("In json input file: netdcf_files has ambiguous format. Specify netcdf_file_sets for two separate ensemble inputs and set netcdf_files accordingly (one ensemble in each sublist)")
        os._exit(1)
