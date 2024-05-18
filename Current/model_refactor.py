# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 15:40:33 2023

@author: Computer
"""

### 1) Data import
## - Add data validation
import json
import pandas as pd
import numpy as np
import os
import copy
import time
import csv
from itertools import product
import multiprocessing.pool
from joblib import Parallel, delayed
from functools import partial
from datetime import datetime
from util.sim_data_prep import filter_acid_data, remove_acids, derive_solution_sets, load_genetic_code, load_acid_data, transform_acid_data, validate_data
start_time = time.time()
# ----------------------------------------------
# INITIAL SETUP
# ----------------------------------------------

# Set up directory paths
# Check if __file__ is defined, if not, set it manually
if '__file__' in globals():
    script_dir = os.path.dirname(os.path.realpath(__file__))
else:
    # Manually define the __file__ if it's not present in the globals
   script_dir = os.path.abspath("C:/Users/Computer/Onedrive/Documents/GitHub/Amino-Acid-Project")

data_dir = os.path.join(script_dir, 'Acid_Base_Data')

# ----------------------------------------------
# DATA LOADING FUNCTIONS
# ----------------------------------------------

# ----------------------------------------------
# DATA VALIDATION
# ----------------------------------------------

# ----------------------------------------------
# MAIN EXECUTION
# ----------------------------------------------

# Load and validate data
acid_processing_df_dict = load_acid_data(data_dir)
validate_data(acid_processing_df_dict)

# Transform data
acid_df_dict = transform_acid_data(acid_processing_df_dict)

print("Data loading and transformation complete!")

### 3) Solution set derivation
## Could store results but that eliminates customization
## Real issue is that this is tied to model assumptions (max graph size (nodes))
# Define the function for deriving solution sets
# Compute and store the solution sets if not already done
if not os.path.exists("solution_sets.json"):
    l, m = 0, 30

    solution_set_derivation, acid_length_dict = derive_solution_sets(l, m)

    # Store the derived solution sets
    with open("solution_sets.json", "w") as file:
        json.dump({
            'solution_set_derivation': solution_set_derivation,
            'acid_length_dict': acid_length_dict
        }, file)

# Load the solution sets when you run the algorithm
with open("solution_sets.json", "r") as file:
    """
    Load the solution sets from a JSON file.
    
    Read the 'solution_sets.json' file to load the precomputed 
    solution sets and acid length associations.
    """
    
    data = json.load(file)
    solution_set_derivation = data['solution_set_derivation']
    acid_length_dict = data['acid_length_dict']

### 4) Choose genetic codes to test
# Load the codes to test
## Old model could do one at a time, we want to generalize...?
sequence_dict, full_acid_list, solution_set_mappings, solution_set_mappings_reduced = load_genetic_code('Standard')
# Building the solution set dictionary
solution_set_dict = {}
for acid in full_acid_list:
    solution_set_name = solution_set_mappings[acid]
    solution_set_dict[acid] = solution_set_derivation[solution_set_name]
### 6) Select range to test over
## Unconstrained
## non-decreasing (i.e. y >= x, r >= q, and z >= w)
## y > x, r > q, or z > w.
## Anchoring specific pairs, e.g. (0, 3, *)
def generate_permutations(criteria_pair=(0, 3)):
    """
    Generate valid permutations for Interpretation Framework (IF) based on specific criteria.
    
    This function generates 6-tuple permutations where each tuple represents slicing indices.
    A reverse reading is implied when the first number of a pair is larger than the second.
    
    :param criteria_pair: tuple, default=(0, 3)
        The pair which the first elements of the permutation should match.
        
    :return: list
        A list of valid 6-tuple permutations that satisfy the specified conditions.
    """
    valid_permutations = []
    
    for permutation in product(range(4), repeat=6):
        if (permutation[1] >= permutation[0] and  # Ensure second element is not less than the first in each pair
            permutation[3] >= permutation[2] and
            permutation[5] >= permutation[4] and
            permutation[:2] == criteria_pair and  # The first pair must match the specified criteria_pair
            permutation[2] != permutation[3] and  # Additional condition: pairs must be non-identical
            permutation[4] != permutation[5]):
            
            # Additional check: If any pair is (0, 3), at least one of the other two pairs must be non-empty
            if ((permutation[0] == 0 and permutation[1] == 3 and
                 (permutation[2] != 0 or permutation[3] != 0 or permutation[4] != 0 or permutation[5] != 0)) or
                (permutation[2] == 0 and permutation[3] == 3 and
                 (permutation[0] != 0 or permutation[1] != 0 or permutation[4] != 0 or permutation[5] != 0)) or
                (permutation[4] == 0 and permutation[5] == 3 and
                 (permutation[0] != 0 or permutation[1] != 0 or permutation[2] != 0 or permutation[3] != 0))):
                valid_permutations.append(permutation)

    return valid_permutations

# Example usage:
permutation_results = generate_permutations(criteria_pair=(0, 3))
answer_dict = {}

# Order in which amino acids will be tested
gate_order = [
    'Glycine', 'Proline', 'Alanine', "Arginine", 'Aspartic Acid', 'Asparagine', 
    'Cysteine', 'Glutamic Acid', 'Glutamine', 'Histidine', 'Isoleucine', 'Leucine', 
    'Lysine', 'Methionine', 'Phenylalanine', #'Selenocysteine',
    'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'
]
"""
gate_order is a predefined order list of amino acids.

This list determines the sequence in which the acids will be processed or tested in the simulation. 
The order can be essential because some acids might have dependencies or effects on subsequent acids. 

Notes:
- 'Selenocysteine' is commented out, which means it's currently excluded from the simulation but can 
  be included if desired by uncommenting it.
- The sequence of acids in this list can be modified based on experimental or scientific requirements.
"""

# ----------------------------------------------
# Modifications to the Run
# ----------------------------------------------
useReducedRun = True

if useReducedRun:
    # Remove specific amino acids from the full acid list for a reduced run.
    # In this case, Glycine and Proline are removed.
    # This step is crucial for focusing the simulation on a subset of amino acids,
    # potentially reducing computational complexity and focusing on more relevant cases.
    acids_to_remove = ["Glycine", "Proline"]
    full_acid_list = remove_acids(full_acid_list, acids_to_remove)

    # Rebuild the solution set dictionary based on the reduced list of amino acids.
    # This step ensures that the solution sets align with the updated list of acids,
    # reflecting the removal of Glycine and Proline.
    solution_set_dict = {}
    for acid in full_acid_list:
        solution_set_name = solution_set_mappings_reduced[acid]
        solution_set_dict[acid] = solution_set_derivation[solution_set_name]
    

if useReducedRun:
    # Determine the threshold value based on Glycine's properties, specifically the number of atoms.
    # Here, a threshold of 10 is chosen based on the understanding that Glycine contains 10 atoms at physiological pH.
    # This threshold is used to filter out rows in the acid dataframes that correspond to or are influenced by Glycine.
    # Atoms started at 1, so we use 10 here instead of 9 (so no zero-indexing)
    ## PRECEDING IS ALL WRONG
    ## ONLY THE FIRST C=O and C-O-H are identical in a meaningful fashion
    ## Since we run into Glycine's second C leading to H instead of C like most other acids do...
    ## TODO: This means we also have to update the genetic_code values for which solution_set to use...
    glycine_threshold = 9

    # Apply the filtering to the dataframes.
    # This step removes rows from each amino acid's dataframe where 'Normalized_node_1' is below the glycine_threshold.
    # It helps in refining the data to exclude influences or characteristics of Glycine.
    acid_df_dict = filter_acid_data(acid_df_dict, glycine_threshold)
    
if useReducedRun:
    # Update the gate order for the reduced run.
    # This step excludes Glycine and Proline from the processing sequence, aligning with the earlier removal of these acids.
    # The updated gate order reflects the focus on a specific subset of amino acids,
    # potentially leading to more targeted and efficient processing in the simulation.
    gate_order = [
    'Alanine', "Arginine", 'Aspartic Acid', 'Asparagine', 
    'Cysteine', 'Glutamic Acid', 'Glutamine', 'Histidine', 'Isoleucine', 'Leucine', 
    'Lysine', 'Methionine', 'Phenylalanine', #'Selenocysteine',
    'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'
    ]

## Main algorithm
def test_framework(possible_IF, gate_order, full_acid_list, result_folder):
    # This shouldn't do much, but I'm curious, so leaving it in.
    run_counter = 0
    ## TODO: Only for local test
    possible_IF = (3, 1, 3, 0, 2, 0)
    # Extract the variables from the possible_IF tuple
    x, y, q, r, w, z = possible_IF
    # Create the slices list
    slices = [(x, y), (q, r), (w, z)]
    print('Current IF:', x, y, q, r, w, z)
    IF_name = f'{x, y, q, r, w, z}'
    # Initialize dictionaries for storing time_dict_X data
    time_dict_1, time_dict_2, time_dict_3 = {}, {}, {}

    # Function to generate key based on forward or reverse reading
    def generate_key(sequence, start, end):
        """
        Generate a key from a sequence based on start and end indices. Handles both forward and reverse readings.

        If the start index is less than or equal to the end index, the function slices the sequence normally.
        If the start index is greater than the end index, it indicates a reverse reading. In this case,
        the function slices the sequence in reverse order.

        This function is designed to handle indices that may extend beyond the current sequence, 
        allowing for referencing neighboring codons in future expansions of the Interpretation Framework.

        Parameters:
        - sequence: str
            The sequence from which to generate the key.
        - start: int
            The starting index for slicing the sequence.
        - end: int
            The ending index for slicing the sequence.

        Returns:
        - str
            A substring of the sequence, sliced according to the specified indices, handling reverse reading if required.
        """
        return sequence[start:end] if start <= end else sequence[end:start][::-1]

    # Generate time_dict_X dictionaries for the selected IF
    for acid in full_acid_list:
        acid_of_interest = sequence_dict[acid]

        # Iterate through the sequences of the selected acid
        for sequence in acid_of_interest:
             # Generate keys for each slice pair in slices.
            # The generate_key function is used to handle both normal and reverse readings
            # based on the order of the start and end indices in each slice pair.
            keys = [generate_key(sequence, *slice_pair) for slice_pair in slices]
            # Pair each generated key with its corresponding time dictionary.
            # This loop iterates over the generated keys and the time dictionaries (time_dict_1, time_dict_2, time_dict_3).
            # For each key, if it is not an empty string, a new DataFrame is assigned to that key in the respective time dictionary.
            for key, time_dict in zip(keys, [time_dict_1, time_dict_2, time_dict_3]):
                if key:  # Check if the key is not an empty string
                    time_dict[key] = pd.DataFrame()

    # Generate save_state_dict for each acid
    save_state_dict = {}
    for acid in full_acid_list:
        save_state_dict[acid] = [
            copy.deepcopy(time_dict_1),
            copy.deepcopy(time_dict_2),
            copy.deepcopy(time_dict_3),
        ]

    # Create an empty_save_state_dict as a reference point
    empty_save_state_dict = copy.deepcopy(save_state_dict)

    # Initialize the save_state_dict with the time_dict_Xs for each amino acid
    save_state_dict = {
        acid: [time_dict_1, time_dict_2, time_dict_3]
        for acid in full_acid_list
    }
    
    ### Function definition ####
    # Tracks gate position
    gate_position = 0
    # Tracks position within solution_sets (THIS USED TO BE "o")
    solution_set_position = 0
    # Working values for the current run
    acid_solution_set = []
    # Where we store the values for all solutions
    answer_list = []
    # Global variables
    current_acid = gate_order[gate_position]
    #acid_df = acid_df_dict[current_acid]
    bases_current_acid = sequence_dict[current_acid]
    testing_sets = solution_set_dict[current_acid]
    testing_set = testing_sets[solution_set_position]
    # Global dictionary to store the history of empty entries
    empty_entries_history = {}

    def gate_change(gate_order, gate_position, acid_df_dict, sequence_dict, solution_set_dict):
        """
        Updates the testing parameters based on the current amino acid.

        Parameters:
        - gate_order (list): A list of amino acids included in the simulation.
        - gate_position (int): Index to access the current acid being evaluated.
        - acid_df_dict (dict): Dictionary storing acid-related data (currently commented out).
        - sequence_dict (dict): Dictionary storing associated codons for each amino acid.
        - solution_set_dict (dict): Dictionary storing all possible subdivisions of the hypergraph representing each amino acid.

        Returns:
        - tuple: 
            * current_acid (str): The amino acid being evaluated.
            * bases_current_acid (list): Codons associated with the current acid.
            * testing_sets (list of lists): Different ways to subdivide the hypergraph for the current acid.
            
        Note:
        - `gate_position` only changes by +1 or -1 to navigate through the acids.
        - `sequence_dict` holds the codons for the currently tested genetic code (current_acid is the key).
        - Each `testing_sets` is a list of lists, where each child list holds three lists representing potential subdivisions of the amino acid's hypergraph.
        """
        current_acid = gate_order[gate_position]
        #acid_df = acid_df_dict[current_acid]
        bases_current_acid = sequence_dict[current_acid]
        testing_sets = solution_set_dict[current_acid]
        return current_acid, bases_current_acid, testing_sets
    
    def testing_set_update(testing_sets, solution_set_position):
        """
        Fetches the current subdivision of the hypergraph for evaluation.

        Parameters:
        - testing_sets (list of lists): All possible ways to subdivide the hypergraph for a particular amino acid.
        - solution_set_position (int): Index to retrieve the current subdivision for evaluation.

        Returns:
        - list: A single testing_set representing a subdivision of the hypergraph for the current stage of computation.
        """
        testing_set = testing_sets[solution_set_position]
        return testing_set
    
    # Function to reset time_dicts based on the history of empty entries
    def reset_time_dicts_to_history(time_dict_1, time_dict_2, time_dict_3, gate_position):
        """
        Resets the time_dicts to their state at a given gate position based on the history of empty entries.

        Parameters:
        - time_dict_1, time_dict_2, time_dict_3 (dict): The time dictionaries to be reset.
        - gate_position (int): The gate position to which the time_dicts will be reset.
        """
        nonlocal empty_entries_history
        if gate_position in empty_entries_history:
            empty_entries_1, empty_entries_2, empty_entries_3 = empty_entries_history[gate_position]
            for key in empty_entries_1:
                time_dict_1[key] = pd.DataFrame()
            for key in empty_entries_2:
                time_dict_2[key] = pd.DataFrame()
            for key in empty_entries_3:
                time_dict_3[key] = pd.DataFrame()

    def process_time_dict(time_dict, time_dict_name, testing_set, testing_set_pos, slice_1, slice_2, bases_current_acid, current_acid, solution_set_position, run_counter):
        """
        Iteratively process each time dictionary and update positions based on returned values.
        
        This function calls the process_time_dict function for each time dictionary, passing 
        relevant parameters. If all the time dictionaries are successfully processed without 
        conflicts, the gate_change function is invoked to update certain values for the next iteration.
        
        Parameters:
        - solution_set_position: int
            Current position within the solution set.
        - gate_position: int
            Position or index of the current gate.
        - acid_solution_set: list
            Set of solutions corresponding to the acids.
        - time_dict_1, time_dict_2, time_dict_3: dict
            The time dictionaries to be processed.
        - bases_current_acid: list
            List of bases corresponding to the current acid.
        - current_acid: str
            The current acid being processed.
        - run_counter: int
            Counter for the number of runs or iterations.
        - testing_sets: list
            List of testing sets.
        - x, y, q, r, w, z: int
            Parameters used for slicing and processing within the process_time_dict function.
        
        Returns:
        tuple: (solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3)
            Updated values based on the processing of the time dictionaries and the gate change.
        """
        norm_dict = create_norm_dict(testing_set[testing_set_pos])
        # TODO: Check work on: Consider moving this outsie of the for loop to avoid recreating it.
        df_selection = create_df_selection(current_acid, testing_set[testing_set_pos], norm_dict, acid_df_dict)
        for base in bases_current_acid:
                #key = base[slice_1:slice_2]  # The key for the current segment of the base
                key = generate_key(base, slice_1, slice_2)
                # Happens when IF does not use the time dict
                if key == "":
                    return solution_set_position, None, True
                    #current_df = time_dict.get(key)  # Safely get the value from the dictionary
                    #if current_df is None:
                    #    return solution_set_position, None, True
                
                current_df = time_dict[key]  # The DataFrame corresponding to the current segment of the base, reaching here means it's not an empty step.

                #print(current_df, key, current_acid, bases_current_acid, time_dict, time_dict_name, IF_name, testing_set[testing_set_pos], testing_set, testing_set_pos)
                if testing_set[testing_set_pos] == []:
                    if not current_df.empty:
                        print(f"Overlap error in t_d_{testing_set_pos + 1} for {current_acid} at solution_set_position = {solution_set_position}. Run counter 1: {run_counter}. acid_solution_set = {acid_solution_set}")
                        solution_set_position += 1
                        return solution_set_position, None, False
                    else:
                        #handle the scenario where current_df is empty (move to the next base)
                        continue
                else:
                    if current_df.empty:
                        ## TODO: Figure out if we have to assign value here for consistency. (Only seems to matter if we're changing the df_selection, but see rework for more on that debate)
                        pass  # This condition is okay; no error here
                    elif current_df.shape != df_selection.shape or not current_df.equals(df_selection):
                        print(f"Overlap error in t_d_{testing_set_pos + 1} for {current_acid} at solution_set_position = {solution_set_position}. Run counter 2: {run_counter}. acid_solution_set = {acid_solution_set}")
                        solution_set_position += 1
                        return solution_set_position, None, False

        return solution_set_position, df_selection, True

    def create_norm_dict(testing_set_step):
        """
        Creates a normalized dictionary for the given testing set step.

        This function creates a dictionary with keys based on the elements of testing_set_step and assigns normalized values 
        (incremented integers starting from 1). Additional entries are created for out-of-bounds lower and upper values to 
        assist in the comparison of graph subsets.

        Parameters:
        - testing_set_step (list): A list representing a testing set step.

        Returns:
        - dict: Normalized dictionary for the given testing set step.

        Example:
        Given `testing_set_step` = [2, 3, 4], the returned `norm_dict` would be:
        {1: 'Out of bounds lower', 2: 1, 3: 2, 4: 3, 5: 'Out of bounds upper'}
        """
        norm_dict = {}
        norm = 1
        for i in range(len(testing_set_step)):
            norm_dict[testing_set_step[i]] = norm
            norm += 1
        if len(testing_set_step) > 0:
            for i in range(1, min(norm_dict)):
                norm_dict[i] = "Out of bounds lower"
            for i in range(max(norm_dict) + 1, testing_set_step[-1] + 1):
                norm_dict[i] = "Out of bounds upper"
        return norm_dict

    ## This drops the unnecessary rows and columns, leaving only the targeted and normalized selection
    def create_df_selection(current_acid, testing_set_step, norm_dict, acid_df_dict):
        """
        Creates a DataFrame selection based on the specified acid and testing set step, with normalization applied.
        
        Parameters:
        - current_acid (str): The name of the current acid being processed.
        - testing_set_step (list): The current testing set step containing node range to be selected.
        - norm_dict (dict): Dictionary mapping original node numbers to their normalized values.
        - acid_df_dict (dict): Dictionary containing data for different acids, with acid names as keys.
        
        Returns:
        - DataFrame: The selected and normalized data for the given acid.
        """
        #Seeing if this is why we're getting nonetype (since we're trying to create normalization for empty steps)
        if (testing_set_step == []):
            return pd.DataFrame()
        # Extract the relevant part of the dataframe from acid_df_dict based on testing_set_step range
        # This used to filter on node_1
        df_selection = acid_df_dict[current_acid][acid_df_dict[current_acid].Normalized_node_1.isin(range(min(testing_set_step), max(testing_set_step) + 1))]
        
        # Drop the 'Node_1' and 'Node_2' columns
        #df_selection = df_selection.drop(columns=['Node_1', 'Node_2'])
        
        # Replace 'Normalized_node_1' and 'Normalized_node_2' values using the norm_dict for normalization
        df_selection = df_selection.replace({'Normalized_node_1': norm_dict, 'Normalized_node_2': norm_dict})
    
        return df_selection
    
    def record_empty_entries(time_dict):
        """
        Records the keys of empty entries in the given time dictionary.
        
        Parameters:
        - time_dict (dict): The time dictionary to check for empty entries.
        
        Returns:
        - set: A set of keys corresponding to empty entries in the time_dict.
        """
        empty_keys = set()
        for key, value in time_dict.items():
            if value.empty:
                empty_keys.add(key)
        return empty_keys

    def actual_run(solution_set_position, gate_position, acid_solution_set, time_dict_1, time_dict_2, time_dict_3, bases_current_acid, current_acid, run_counter, testing_sets, x, y, q, r, w, z):
        """
        Iteratively process each time dictionary and update positions based on returned values.
        
        This function calls the process_time_dict function for each time dictionary, passing 
        relevant parameters. If all the time dictionaries are successfully processed without 
        conflicts, the gate_change function is invoked to update certain values for the next iteration.
        
        Parameters:
        - solution_set_position: int
            Current position within the solution set.
        - gate_position: int
            Position or index of the current gate.
        - acid_solution_set: list
            Set of solutions corresponding to the acids.
        - time_dict_1, time_dict_2, time_dict_3: dict
            The time dictionaries to be processed.
        - bases_current_acid: list
            List of bases corresponding to the current acid.
        - current_acid: str
            The current acid being processed.
        - run_counter: int
            Counter for the number of runs or iterations.
        - testing_sets: list
            List of testing sets.
        - x, y, q, r, w, z: int
            Parameters used for slicing and processing within the process_time_dict function.
        
        Returns:
        tuple: (solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3)
            Updated values based on the processing of the time dictionaries and the gate change.
        """
        nonlocal empty_entries_history
        solution_set_position, selection_1, success = process_time_dict(time_dict_1, "time_dict_1", testing_set, 0, x, y, bases_current_acid, current_acid, solution_set_position, run_counter)
        if not success:
            return solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3
    
        solution_set_position, selection_2, success = process_time_dict(time_dict_2, "time_dict_2", testing_set, 1, q, r, bases_current_acid, current_acid, solution_set_position, run_counter)
        if not success:
            return solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3
    
        solution_set_position, selection_3, success = process_time_dict(time_dict_3, "time_dict_3", testing_set, 2, w, z, bases_current_acid, current_acid, solution_set_position, run_counter)
        if not success:
            return solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3
    
        #Before we modify the time_dicts, we need to establish which entries are empty at this point.
        #If the entry is empty here, then when we step back to this point, the entries need to be returned to that state.
        # Record the state of empty entries in each time_dict before updating them
        empty_entries_1 = record_empty_entries(time_dict_1)
        empty_entries_2 = record_empty_entries(time_dict_2)
        empty_entries_3 = record_empty_entries(time_dict_3)

        # Now, we can modify the time_dicts with the validated selections:
        # TODO: This will unnecessarily rewrite already full key-values...
        # Modify the time_dicts with the validated selections from process_time_dict
        for base in bases_current_acid:
            # Generate keys using the generate_key function
            key_1 = generate_key(base, x, y)
            key_2 = generate_key(base, q, r)
            key_3 = generate_key(base, w, z)
            # Only check if the string slice (e.g., base[x:y]) is not empty. 
            # In the process_time_dict function, when the key (derived from the slice like base[x:y]) 
            # is an empty string, it returns None for the selection (e.g., selection_1). 
            # Hence, if the slice is not empty, it guarantees that the corresponding selection 
            # is not None and is a valid DataFrame to be assigned to the time_dict. 
            # This approach avoids redundant checks and keeps the code concise and efficient.
            # Assign the validated selections to the time_dicts if the keys are not empty
            if key_1:
                time_dict_1[key_1] = selection_1
            if key_2:
                time_dict_2[key_2] = selection_2
            if key_3:
                time_dict_3[key_3] = selection_3

        # Update the history dictionary with the current gate's empty entries
        empty_entries_history[gate_position] = (empty_entries_1, empty_entries_2, empty_entries_3)

        acid_solution_set.append(solution_set_position)
        print(f"Value found for {current_acid} at gate_position = {gate_position}; solution_set_position = {solution_set_position}")
        #save_state_dict[current_acid] = copy.deepcopy([time_dict_1, time_dict_2, time_dict_3])
        #TODO: Ensure that we need to set solution_set_position to 0 after moving to the next gate. (Need to set at 0, so correct value is used since it's a new gate and so new statistical combination)
        solution_set_position = 0
        gate_position += 1
    
        if gate_position < len(gate_order):
            current_acid, bases_current_acid, testing_sets = gate_change(gate_order, gate_position, acid_df_dict, sequence_dict, solution_set_dict)
        elif gate_position > len(gate_order) - 1:
            print(f"Solution found b/c gate_order exceeds max length. gate_position: {gate_position} len(gate_order - 1): {len(gate_order) - 1}")
            testing_sets = "MADE UP VALUE"
    
        return solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3
  
    ########## The code that executes the search ##########
    counter = 0
    #Our running condition
    while gate_position >= 0:
        print(f'Counter iteration: {counter}')
        #This is the first condition run, it is only run once per Interpretation Framework.
        if gate_position == 0 and solution_set_position == 0:
            solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3 = actual_run(solution_set_position, gate_position, acid_solution_set, time_dict_1, time_dict_2, time_dict_3, bases_current_acid, current_acid, counter, testing_sets, x, y, q, r, w, z)
            counter += 1
    
        #This is the condition that runs the most. First expression (below upper bound), second expression (exist testing set: need to check - 1) 
        elif gate_position <= len(gate_order) - 1 and (solution_set_position < len(testing_sets) and solution_set_position > 0):
            #Likely modified a t_d_x in previous run, so need to load the last working 
            ## TODO: Should be able to eliminate this in some scenarios (unless a gate change has just happened)?
            ## If we move the load logic to right after gate_change, then we should be able to eliminate
            #load_state_dict = load_save_state(gate_order, gate_position, solution_set_position, empty_save_state_dict, save_state_dict)
            #time_dict_1, time_dict_2, time_dict_3 = load_state_dict["time_dict_1"], load_state_dict["time_dict_2"], load_state_dict["time_dict_3"]
            # Unpack time_dicts from the loaded state dictionary for better clarity
            testing_set = testing_set_update(testing_sets, solution_set_position)
            print(f"({IF_name})acid: {current_acid}. gate_position: {gate_position}. solution_set_position: {solution_set_position}.")
            solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3 = actual_run(solution_set_position, gate_position, acid_solution_set, time_dict_1, time_dict_2, time_dict_3, bases_current_acid, current_acid, counter, testing_sets, x, y, q, r, w, z)
            counter += 1 
    
        #This condition occurs right after a gate_change() and so there is no need to call load_save_state
        elif gate_position <= len(gate_order) - 1 and solution_set_position == 0:
            testing_set = testing_set_update(testing_sets, solution_set_position)
            print(f"<THIS IS WILD> solution_set_position: {solution_set_position}")
            print(f"acid: {current_acid}. gate_position: {gate_position}. solution_set_position: {solution_set_position}.")
            solution_set_position, gate_position, current_acid, bases_current_acid, testing_sets, time_dict_1, time_dict_2, time_dict_3 = actual_run(solution_set_position, gate_position, acid_solution_set, time_dict_1, time_dict_2, time_dict_3, bases_current_acid, current_acid, counter, testing_sets, x, y, q, r, w, z)
            counter += 1 
    
        #If this condition is met, we've found a solution (tested every gate and found no conflicts, so the last gate finding a working value pushes us past len(gate_order) - 1)    
        elif gate_position > len(gate_order) - 1:
            print("An answer has been found using acid_solution_set = ", acid_solution_set)
            answer = copy.deepcopy(acid_solution_set)
            answer_list.append(answer)
            print("Number of answers is now: ", len(answer_list))
            ## Adding this code to safely test all permutations
            # Write the result to a unique file
            result_file_path = os.path.join("run_data", result_folder, f"IF_{possible_IF}.csv")
            with open(result_file_path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([IF_name, answer_list, counter])
        
            return [IF_name, answer_list, counter]
            #Return to previous gate to check for more solutions
            gate_position = gate_position - 1
            #We've changed gates, so we update our values
            current_acid, bases_current_acid, testing_sets = gate_change(gate_order, gate_position, acid_df_dict, sequence_dict, solution_set_dict)
            #Pick up from where we left off in previous acid
            solution_set_position = acid_solution_set[-1]
            #Remove the last value to search for another solution
            acid_solution_set.pop()
            #Increment by 1, so we're checking a new combination
            solution_set_position = solution_set_position + 1
            #Now we load the save_state
            ## Empty out the values from the time_dict for current acid because we need to look for more solutions
            reset_time_dicts_to_history(time_dict_1, time_dict_2, time_dict_3, gate_position)
            continue
    
        #This is our ending condition for the current Interpretation framework.
        #It stores information for later analysis.
        #Should this be len(testing_sets) - 1? (Pretty sure it shouldn't be, but will check again)
        elif gate_position == 0 and solution_set_position >= len(testing_sets):
            print(f"Program has finished running. There were {len(answer_list)} answers found. Counter: {counter}")
            #This is where we put code to move between overlap_groups if we're testing acid combinations separately.
            answer_dict[IF_name] = answer_list
            ## results_df.loc[len(results_df.index)] = [IF_name, answer_list, counter]
            run_counter += 1
            return [IF_name, answer_list, counter]
        #All testing_set values for given gate have been exhausted. Move back to previous gate and increment by 1.
        elif solution_set_position >= len(testing_sets) :
            print(f"No solution found for {current_acid}. Moving to previous gate.")
            gate_position = gate_position - 1
            #The gate has changed, so we must once again update
            ## If we remove the load from the earlier elif, then we need to put it here
            current_acid, bases_current_acid, testing_sets = gate_change(gate_order, gate_position, acid_df_dict, sequence_dict, solution_set_dict)
            #Set value for testing set
            solution_set_position = acid_solution_set[-1]
            #Remove old value
            acid_solution_set.pop()
            # increment solution_set_position by 1 to avoid rerunning the answer, then we need to start this again, or check to see if 
            # we're > len(solution_sets)
            solution_set_position = solution_set_position + 1
            #Load_save_state() here with reset_time_dict.
            reset_time_dicts_to_history(time_dict_1, time_dict_2, time_dict_3, gate_position)
            
        else:
            print("An unaccounted for scenario occurred. Exiting program.")
            print(f"current_acid: {current_acid} gate_position: {gate_position} solution_set_position: {solution_set_position}")
            break

def create_unique_result_folder(timestamp, run_type_label='default_run'):
    """
    Create a unique folder for storing results of the current run with a specified run type label.

    Args:
    - run_type_label (str): A label describing the type of the run, e.g., 'full_run', 'reduced_run'.

    Returns:
    - str: The path of the created folder, used to store results from the current run.
    """
    folder_name = f"{run_type_label}_{timestamp}"
    result_folder_path = os.path.join("run_data", folder_name)
    os.makedirs(result_folder_path, exist_ok=True)
    return result_folder_path

# Create a partial function with gate_order and full_acid_list pre-filled
partial_test_framework = partial(test_framework, gate_order=gate_order, full_acid_list=full_acid_list)

def run_multiprocessing(permutation_results, result_folder):
    # Determine the number of cores
    num_cores = multiprocessing.cpu_count()
    
    # Execute the jobs in parallel
    results = Parallel(n_jobs=num_cores - 1)(delayed(partial_test_framework)(possible_IF=i, result_folder=result_folder) for i in permutation_results)
    #results = Parallel(n_jobs=1)(delayed(partial_test_framework)(possible_IF=i) for i in permutation_results)

    # Execute the jobs in parallel
    #results = Parallel(n_jobs=num_cores - 1)(
    #results = Parallel(n_jobs=1)(
    #    delayed(partial_test_framework)(possible_IF=i) for i in permutation_results
    #)
    
    # Convert results into a DataFrame
    results_df = pd.DataFrame(results, columns=['IF', 'Solution_array', 'Counter'])
    
    return results_df

def save_to_json(data, filename, folder):
    with open(os.path.join(folder, filename), 'w') as file:
        json.dump(data, file, indent=4)

def save_to_csv(data, filename, folder):
    if isinstance(data, pd.DataFrame):
        data.to_csv(os.path.join(folder, filename), index=False)
    else:
        pd.DataFrame(data).to_csv(os.path.join(folder, filename), index=False)

if __name__ == "__main__":
    #run_type_label = "reduced_run_full_exit_on_success"  # Or "full_run", "test_run", etc. 
    run_type_label = "reduced_run_03_exit_on_success"  # Or "full_run", "test_run", etc. 
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    result_folder = create_unique_result_folder(run_type_label)
    # Run the multiprocessing
    results_df = run_multiprocessing(permutation_results, result_folder)

    # Save key configuration data
    save_to_json(gate_order, f"gate_order_{current_time}.json", result_folder)
    save_to_json(sequence_dict, f"sequence_dict_{current_time}.json", result_folder)
    save_to_json(full_acid_list, f"full_acid_list_{current_time}.json", result_folder)
    save_to_json(solution_set_dict, f"solution_set_dict_{current_time}.json", result_folder)

    # Save the results DataFrame
    save_to_csv(results_df, f"results_{current_time}.csv", result_folder)

    # Calculate and print the computation time
    end_time = time.time()
    difference = end_time - start_time
    hours, minutes, seconds = difference // 3600, (difference % 3600) // 60, difference % 60
    print(f"Computation took {hours} hours, {minutes} minutes, and {seconds} seconds.")