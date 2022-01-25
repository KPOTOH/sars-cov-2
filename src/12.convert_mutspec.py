"""
format mutspec table to such style for SigProfilerExtractor

Mutation Types      Sample1 Sample2
A[C>A]A             58      74
A[C>A]C             36      66
A[C>A]G             13      12
A[C>A]T             37      64

C>A, C>G, C>T, T>A, T>C, T>G
"""

import numpy as np
import pandas as pd

PATH_TO_MUTSPEC_192 = "./data/share/07.All_MutSpec192_ForFullGenome.csv"
nucs_that_mutated = {"C", "U"}

translator = str.maketrans("AUGC", "TACG")

mutspec192 = pd.read_csv(PATH_TO_MUTSPEC_192)
print(mutspec192)


def collapse_mutspec(nuc_subst: str):
    """ return substitution in A[A>C]A style

    nuc_subst: AAA>ACA format"""
    assert isinstance(nuc_subst, str) and len(nuc_subst) == 7

    mutated_nuc = nuc_subst[1]
    if mutated_nuc not in nucs_that_mutated:
        nuc_subst = nuc_subst.translate(translator)
    else:
        nuc_subst = nuc_subst.replace("U", "T")

    down_nuc = nuc_subst[0]
    up_nuc = nuc_subst[-1]
    mut = nuc_subst[1] + ">" + nuc_subst[-2]
    new_style_subst = down_nuc + "[" + mut + "]" + up_nuc
    return new_style_subst


mut_types = mutspec192.NucSubst.apply(collapse_mutspec)
print(mut_types)


# mutspec96 = pd.DataFrame()
