import json
import os
import pandas as pd
import numpy as np
import copy

# Set up directory paths
if '__file__' in globals():
    script_dir = os.path.dirname(os.path.realpath(__file__))
else:
    # Manually define the script_dir if '__file__' is not present in the globals
    script_dir = os.path.abspath("C:/Users/Computer/Documents/GitHub/Amino-Acid-Project")

# Define the directory for the original IF data
data_dir = os.path.join(script_dir, 'IF_Data')

# Define the directory for the groupings data
groupings_dir = os.path.join(script_dir, 'IF_Groupings')

# Define the directory for the acid data
acid_data_dir = os.path.join(script_dir, 'Acid_Base_Data')

def load_if_data(if_values):
    """
    Load the original IF data file.

    Args:
        if_values (tuple): A tuple representing the IF values.

    Returns:
        dict: The data loaded from the IF JSON file.
    """
    # Calculate the project root directory
    if '__file__' in globals():
        project_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    else:
        project_root = os.path.abspath("C:/Users/Computer/Documents/GitHub/Amino-Acid-Project")

    data_dir = os.path.join(project_root, 'IF_Data')
    file_name = f"IF_({', '.join(map(str, if_values))})_data.json"
    file_path = os.path.join(data_dir, file_name)

    with open(file_path, 'r') as file:
        return json.load(file)

def load_groupings_json(if_values, selected_time_dict):
    """
    Load the groupings JSON file for a specified IF and filter for the selected time_dict.

    Args:
        if_values (tuple): A tuple representing the IF values.
        selected_time_dict (str): The selected time_dict key as a string.

    Returns:
        pd.DataFrame: The DataFrame loaded from the groupings JSON file, filtered for the selected time_dict.
    """

    # Check if __file__ is defined, if not, set it manually
    if '__file__' in globals():
        # Calculate the project root directory
        project_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    else:
        # Manually define the project root directory if __file__ isn't defined
        project_root = os.path.abspath("C:/Users/Computer/Documents/GitHub/Amino-Acid-Project")

    groupings_dir = os.path.join(project_root, 'IF_Groupings')
    file_name = f"IF_({', '.join(map(str, if_values))})_groupings.json"
    file_path = os.path.join(groupings_dir, file_name)

    groupings_df = pd.read_json(file_path, orient='records', lines=True)
    
    # Filter the DataFrame for the selected time_dict
    groupings_df = groupings_df[groupings_df["Time Dict"] == selected_time_dict]
    
    return groupings_df

def load_acid_data(filepath):
    """
    Load acid data from the provided path.

    Parameters:
    - filepath: str
        The path to the data directory or file.

    Returns:
    - dict
        Dictionary with acid names as keys and corresponding dataframes as values.
    """
    acids = [
        "Alanine", "Arginine", "Asparagine", "Aspartic Acid", "Cysteine", "Glutamic Acid",
        "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
        "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"
    ]

    acid_dfs = {}
    for acid in acids:
        acid_file = os.path.join(filepath, f"{acid}.xlsx")
        if not os.path.exists(acid_file):
            acid_file = os.path.join(filepath, "New", f"{acid}.xlsx")
        acid_dfs[acid] = pd.read_excel(acid_file)

    return acid_dfs

def transform_acid_data(acid_processing_df_dict):
    """
    Convert loaded acid data into a desired format.

    Parameters:
    - acid_processing_df_dict: dict
        Dictionary with raw acid dataframes.

    Returns:
    - dict
        Dictionary with transformed acid dataframes.
    """
    acid_df_dict = {}
    for acid_name, df in acid_processing_df_dict.items():
        conversion_dict = {
            str(node_value): [j.strip(',') for j in df['Edges'][i].split()]
            for i, node_value in enumerate(range(1, len(df) + 1))
        }

        new_df = pd.DataFrame(columns=['Node_1', 'Node_2', 'Node_1_type', 'Node_2_type'])
        for entry, list_of_values in conversion_dict.items():
            for node2 in list_of_values:
                new_df.loc[len(new_df.index)] = [entry, node2, entry, node2]

        type_dict = {f'{i + 1}': df['Type'].loc[df.index[i]] for i in range(len(df))}
        new_df = new_df.replace({'Node_1_type': type_dict, 'Node_2_type': type_dict})
        new_df[['Node_1', 'Node_2']] = new_df[['Node_1', 'Node_2']].astype(int)
        # This is probably unneeded, since we can just use the basic columns (but one change at a time...)
        new_df['Normalized_node_1'] = new_df['Node_1']
        new_df['Normalized_node_2'] = new_df['Node_2']
        # Drop 'Node_1' and 'Node_2' columns as they are no longer needed
        new_df.drop(columns=['Node_1', 'Node_2'], inplace=True)
        acid_df_dict[acid_name] = copy.deepcopy(new_df)

    return acid_df_dict

def derive_solution_sets(l, m):
    """
    Derive all possible solution sets for given range of values.
    
    This function generates and returns solution sets and acid length 
    associations for a specified range of values (l to m). Each solution 
    set is a list containing lists of bucketed integers derived from 
    the original list of integers ranging from 1 to a given value.
    
    Parameters:
    - l (int): The starting value of the range.
    - m (int): The ending value of the range (representing max graph size).
    
    Returns:
    tuple: (solution_set_derivation, acid_length_dict)
    - solution_set_derivation (dict): Dictionary where keys are strings 
      in the format "list_n[i]_k[j]" (where [i] is the length of the 
      original list of integers and [j] is the total number of derived 
      solution sets) and values are lists of solution sets.
    - acid_length_dict (dict): Dictionary where keys are integers 
      representing the length of the original list and values are 
      strings in the format "list_n[i]_k[j]".
    """
    
    solution_set_derivation = {}
    acid_length_dict = {}

    for i in range(l, m + 1):
        hold_list = list(range(1, i + 1))
        solution_set_list = []

        ticker = 0
        while ticker <= i:
            bucket_1, bucket_2, bucket_3 = hold_list[:], [], []

            if ticker > 0:
                bucket_1 = hold_list[:-ticker]
                bucket_2 = hold_list[-ticker:]

            solution_set = [bucket_1, bucket_2, bucket_3]
            solution_set_list.append(copy.deepcopy(solution_set))

            while bucket_2:
                bucket_3.insert(0, bucket_2.pop())
                solution_set = [bucket_1, bucket_2, bucket_3]
                solution_set_list.append(copy.deepcopy(solution_set))

            ticker += 1

        solution_set_derivation[f"list_n{i}_k{len(solution_set_list)}"] = solution_set_list
        acid_length_dict[i] = f"list_n{i}_k{len(solution_set_list)}"

    return solution_set_derivation, acid_length_dict

import os
import json

def load_genetic_code(code_name="Standard"):
    """
    Loads a specific genetic code from the JSON file.
    See: https://en.wikipedia.org/wiki/List_of_genetic_codes

    Args:
        code_name (str): The name of the genetic code to load.

    Returns:
        tuple: A tuple containing the sequence dictionary and the full acid list.
    """

    # Check if __file__ is defined, if not, set it manually
    if '__file__' in globals():
        # If __file__ is defined, calculate the project root directory
        # This assumes the util folder is one level down from the project root
        project_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    else:
        # Manually define the project root directory if __file__ isn't defined
        project_root = os.path.abspath("C:/Users/Computer/Documents/GitHub/Amino-Acid-Project")
    
    # Construct the full path to the JSON file
    json_path = os.path.join(project_root, 'genetic_codes.json')

    with open(json_path, 'r') as file:
        data = json.load(file)

    if code_name not in data:
        raise ValueError(f"Genetic code '{code_name}' not found in the JSON file.")

    sequence_dict = data[code_name]["Sequences"]
    full_acid_list = data[code_name]["FullAcidList"]
    solution_set_mappings = data[code_name]["SolutionSets"]
    solution_set_mappings_reduced = data[code_name]["ReducedSolutionSets"]

    return sequence_dict, full_acid_list, solution_set_mappings, solution_set_mappings_reduced

def remove_acids(full_acid_list, acids_to_remove):
    """
    Removes specific amino acids from the full acid list.

    Args:
        full_acid_list (list): The original list of amino acids.
        acids_to_remove (list): A list of amino acids to remove from the full_acid_list.

    Returns:
        list: A new list with specified amino acids removed.
    """
    return [acid for acid in full_acid_list if acid not in acids_to_remove]

def filter_acid_data(acid_df_dict, threshold):
    """
    Filter out rows from acid dataframes based on a threshold and re-index the remaining nodes.

    Args:
    - acid_df_dict: dict
        Dictionary of acid dataframes.
    - threshold: int
        The threshold value for 'Normalized_node_1' to determine which rows to keep.

    Returns:
    - dict
        Updated dictionary with filtered dataframes.
    """
    filtered_acid_df_dict = {}
    for acid, df in acid_df_dict.items():
        # Filter based on threshold
        filtered_df = df[df['Normalized_node_1'] >= threshold].copy()

        # Re-index remaining nodes
        unique_nodes = pd.unique(filtered_df[['Normalized_node_1', 'Normalized_node_2']].values.ravel('K'))
        node_mapping = {node: i+1 for i, node in enumerate(sorted(unique_nodes))}
        filtered_df.replace({'Normalized_node_1': node_mapping, 'Normalized_node_2': node_mapping}, inplace=True)

        filtered_acid_df_dict[acid] = filtered_df
    return filtered_acid_df_dict

def validate_data(acid_dfs):
    """
    Validate loaded acid data.

    Parameters:
    - acid_dfs: dict
        Dictionary with acid dataframes.

    Raises:
    - ValueError: If required columns or values are missing.
    - TypeError: If values are of unexpected data type.
    """
    for acid_name, df in acid_dfs.items():
        # Check for required columns
        required_columns = ["Edges", "Type"]
        for column in required_columns:
            if column not in df.columns:
                raise ValueError(f"'{column}' column not found in {acid_name} dataframe.")
            
        # Check for data types and missing values
        if not df["Edges"].apply(isinstance, args=(str,)).all():
            raise TypeError(f"Unexpected data type in 'Edges' column of {acid_name} dataframe.")
        if df["Edges"].isna().any():
            raise ValueError(f"Missing values detected in 'Edges' column of {acid_name} dataframe.")
        if df["Type"].isna().any():
            raise ValueError(f"Missing values detected in 'Type' column of {acid_name} dataframe.")