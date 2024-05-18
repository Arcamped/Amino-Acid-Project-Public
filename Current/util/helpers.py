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

def print_test(string):
    """
    Prints and returns the given string
    """
    print(string)
    return string