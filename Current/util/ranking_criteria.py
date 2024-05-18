import pandas as pd

def rank_by_atom_proportion(acid_df_dict, atom_type='Carbon'):
    """
    Rank amino acids based on the proportion of a specific atom type.

    Parameters:
    - acid_df_dict: dict
        Dictionary of amino acid dataframes.
    - atom_type: str
        The type of atom to base the ranking on (default: 'Carbon').

    Returns:
    - DataFrame
        DataFrame with amino acids ranked by the specified atom proportion.
    """
    ranking_data = []

    for acid_name, df in acid_df_dict.items():
        total_atoms = len(df)
        specific_atom_count = len(df[(df['Node_1_type'] == atom_type) | (df['Node_2_type'] == atom_type)])
        proportion = specific_atom_count / total_atoms if total_atoms > 0 else 0
        ranking_data.append((acid_name, proportion))

    ranking_df = pd.DataFrame(ranking_data, columns=['Amino Acid', 'Proportion of ' + atom_type])
    ranking_df.sort_values(by='Proportion of ' + atom_type, ascending=False, inplace=True)

    return ranking_df