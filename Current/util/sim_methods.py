import pandas as pd
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
    - DataFrame: The selected and normalized data for the given acid. Raises error if selection is unexpectedly empty.
    """
    # Check if the testing set step is empty, and return an empty DataFrame if so
    if not testing_set_step:
        return pd.DataFrame()

    # Get the DataFrame for the current acid
    acid_data = acid_df_dict[current_acid]

    # Determine the range of values to filter on
    filter_min, filter_max = min(testing_set_step, default=0), max(testing_set_step, default=0)
    filter_range = range(filter_min, filter_max + 1)

    # Filter the DataFrame based on the range
    df_selection = acid_data[acid_data['Normalized_node_1'].isin(filter_range)]

    # Apply normalization to the selected DataFrame
    df_selection = df_selection.replace({'Normalized_node_1': norm_dict, 'Normalized_node_2': norm_dict})

    # Check if the resulting DataFrame is empty and it shouldn't be (i.e., testing_set_step isn't empty)
    if df_selection.empty and testing_set_step:
        # Handling the error - you might want to log this or raise an exception
        if current_acid == "Alanine":
            print(current_acid, testing_set_step, "acid data: ", acid_data)
        raise ValueError(f"Expected data for {current_acid} within the range {filter_range} but got an empty DataFrame.")
    
    # print(filter_range, df_selection, current_acid, testing_set_step)

    return df_selection

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
        for i in testing_set_step:
            norm_dict[i] = norm
            norm += 1

        if testing_set_step:
            min_val, max_val = min(testing_set_step), max(testing_set_step)
            # Add out-of-bounds entries
            for i in range(1, min_val):
                norm_dict[i] = "Out of bounds lower"
            for i in range(max_val + 1, max_val + 1 + len(testing_set_step)):  # Assuming an equal range on the upper side
                norm_dict[i] = "Out of bounds upper"

        return norm_dict