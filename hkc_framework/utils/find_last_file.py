import os
from string import ascii_letters
def findLargestIndex(path):
    """ Purpose: Find largest file index
        Args: path = path to file 
        Returns : max file index 
    """
    files = os.listdir(path)
    extracted_indicies = []
    for x in files:
        if x == ".keep":
            continue
        extracted_indicies.append(int(x.strip(ascii_letters + '_.')))    
    #files.rstrip(ascii_letters)
    #files.rstrip('_.')
    #files = int(files)
    return(max(extracted_indicies))
