#!/usr/bin/env python 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import DivergingNorm

__aa_seq__ = ("LRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGETGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKEPEERPTFEYLQAFL")
__cd_pos__ = np.arange(270, 520)
__aa_key__ = {
    'His' : 'H',
    'Lys' : 'K',
    'Arg' : 'R',
    'Asp' : 'D',
    'Glu' : 'E',
    'Cys' : 'C',
    'Met' : 'M',
    'Asn' : 'N',
    'Gln' : 'Q',
    'Ser' : 'S',
    'Thr' : 'T',
    'Ala' : 'A',
    'Ile' : 'I',
    'Leu' : 'L',
    'Val' : 'V',
    'Phe' : 'F',
    'Trp' : 'W',
    'Tyr' : 'Y',
    'Gly' : 'G',
    'Pro' : 'P',
    'Ter' : '*'
    
}
__pos_idx_dict__ = dict(zip(__cd_pos__, np.arange(len(__cd_pos__))))
__aa_idx_dict__ = dict(zip(list(__aa_key__.values()), np.arange(len(__aa_key__))))


def create_variant_index(df):
    """
    Create variant index from DataFrame with index in hgvs format.
    
    Args:
        df: DataFrame with index in hgvs format
    Returns:
        variant_index: variant in [1 letter WTAA][position][1 letter MutAA]
                       format, maintaining original order (list)
        
    """
    
    # Grab original hgvs undex
    variant_index = df.index.tolist()
    
    # Loop through all variants
    for idx, variant in enumerate(variant_index):
        
        # Skip _sy and _wt variants
        if 'p.' not in variant:
            continue

        # Get WT, mutant identities and position
        wt_aa = variant[2:5]
        pos = variant[5:-3]
        mut_aa = variant[-3:]

        # Create shortened form
        variant_shortened = __aa_key__[wt_aa] + str(pos) + __aa_key__[mut_aa]

        # Reassign value
        variant_index[idx] = variant_shortened
    
    return(variant_index)


def create_heatmap_arr(scores, variant_index):
    """
    Create heatmap from given scores.
    
    Args:
        scores: scores in order of variant_index (np array)
        variant_index: order of variants in the format 
                        [1 letter WTAA][position][1 letter MutAA] (list)
    Returns:
        heatmap_arr: 2D array where positions=rows, aas=cols
    
    """
    
    # Create empty array for storing values
    heatmap_arr = np.empty((len(__pos_idx_dict__), len(__aa_idx_dict__)))
    heatmap_arr[:] = np.nan
    
    # Loop through individual variants
    for variant, score in zip(variant_index, scores):

        # Skip "_wt" and "_sy" rows
        if len(variant) < 5:
            continue

        # Collect wt, mutant, and position data from variant_index
        wt_aa = variant[0]
        pos = variant[1:-1]
        mut_aa = variant[-1]

        # Assign score to heatmap array
        row_idx = __pos_idx_dict__[int(pos)]
        col_idx = __aa_idx_dict__[mut_aa]
        heatmap_arr[row_idx, col_idx] = score
        
    return(heatmap_arr)


def plot_heatmap(heatmap, title):
    """
    Plot heatmap.
    
    Args:
        heatmap: 2D array where positions=rows, aas=cols
        title: title of plot
    Returns:
        fig: matplotlib.pyplot figure object
        ax: matplotlib.pyplot axis object

    """
    fig, ax = plt.subplots(figsize=(50,300))
    resid_map = plt.imshow(heatmap.T, cmap='bwr', norm=DivergingNorm(0.0))

    # Set tick locations
    ax.set_yticks(np.arange(heatmap.shape[1]))
    ax.set_xticks(np.arange(heatmap.shape[0]))

    # Set tick labels
    ax.set_yticklabels(__aa_idx_dict__.keys())
    ax.set_xticklabels(__pos_idx_dict__.keys())
    
    # Set title
    plt.title(title)
    
    # Show figure
    plt.show()
    
    return(fig, ax)


def heatmap_to_pymol(heatmap):
    """Calculate values for PyMOL recoloring."""
    
    pymol_vals = np.nanmean(heatmap, axis=1)
    return(pymol_vals)