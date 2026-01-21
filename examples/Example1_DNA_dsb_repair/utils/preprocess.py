#Packages
import pandas as pd

#Auxiliary functions
def get_variation_type(sequence,len_ref):
    """
    Returns the type of variation in a sequence (possibly containing '-' gaps)
    compared to a reference length.

    Args:
        sequence (str): Input sequence, possibly containing '-'.
        len_ref (int): Length of the reference sequence.

    Returns:
        str: One of 'none', 'deletion', or 'insertion'.
    """
    len_seq=len(sequence.replace('-',''))
    if len_seq==len_ref:
        return 'none'
    elif len_seq<len_ref:
        return 'deletion'
    else:
        return 'insertion'
    
def inserted_sequence(original,with_insertion):
    """
    Extract the inserted subsequence in `with_insertion`, comparing it with `original`. 
    Both sequences are considered to have the same length, and `original` has '-' 
    at the locations of the inserted subsequence. 

    Args:
        original (str): Reference sequence with placeholders '-'.
        with_insertion (str): Modified sequence with an inserted subsequence.

    Returns:
        str: Inserted subsequence.
    """
    sequence=''
    for l in (with_insertion[i] for i in range(len(original)) if original[i]!=with_insertion[i]):
        sequence+=l
    return sequence

#Main Function to re-format the raw data
def format_experiment_table(file_root, file_name): 
    """
    Load and reformat a experimental table into a clean, analysis-ready format.

    Args:
        file_root (str): Path to the directory containing the raw data file.
        file_name (str): Name of the raw data file (csv).

    Returns:
        pandas.DataFrame: A dataframe with columns:
            - `library` (str): Replicate/library identifier.
            - `sequence` (str): Inserted sequence (empty string for no insertion).
            - `freq` (float): Observed frequency.
    """
    #Load Table
    Raw_Table=pd.read_table(file_root+file_name)
    
    #Format names of replicates
    Raw_Table=Raw_Table.melt(id_vars=['ref_align','read_align'],var_name='library',value_name='freq')
    Raw_Table.library=Raw_Table.library.apply(lambda w: w.split('_')[1])
    
    #Get reference sequence
    ref_sequence=Raw_Table.ref_align.apply(lambda x: x.replace('-','')).unique()[0]
    len_ref=len(ref_sequence)
    
    #Get variation types
    Raw_Table['variation_type']=Raw_Table.read_align.apply(lambda x: get_variation_type(x,len_ref))
    
    #Focus only on insertions and None
    Raw_Table=Raw_Table[Raw_Table['variation_type'].isin(['none','insertion'])].reset_index(drop=True)
    
    #Isolate Inserted Sequence
    Raw_Table.loc[:,'sequence']=\
        [inserted_sequence(Raw_Table['ref_align'][i],Raw_Table['read_align'][i])\
        for i in range(len(Raw_Table))]
    
    #Keep only relevant columns
    Raw_Table=Raw_Table[['library','sequence','freq']]

    return Raw_Table