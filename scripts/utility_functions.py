import os

def get_all_files_in_paths(paths, file_name_substrings=None, file_extension=None):
    file_paths = []
    for path in paths:
        for root, _, files in os.walk(path):
            for file in files:
                if not file_name_substrings or all([substring in file for substring in file_name_substrings]):
                    if not file_extension or file.endswith(file_extension):
                        file_path = os.path.join(root, file)
                        file_paths.append(file_path)
    return file_paths

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
