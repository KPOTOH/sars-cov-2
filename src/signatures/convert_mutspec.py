"""
format mutspec table to such style for SigProfilerExtractor:

Mutation Types      Sample1 Sample2
A[C>A]A             58      74
A[C>A]C             36      66
A[C>A]G             13      12
A[C>A]T             37      64

C>A, C>G, C>T, T>A, T>C, T>G
"""

import random
import pandas as pd

PATH_TO_MUTSPEC_192 = "./data/share/07.All_MutSpec192_ForFullGenome.csv"
PATH_TO_OUT = "./data/share/mutspec_96_all.tsv"
COLS = ["NucSubst" ,"ObsFr" ,"Colour" ,"ExpFr" ,"ObsToExp" ,"Context"]
nucs_that_mutated = {"C", "U"}
translator = str.maketrans("AUGC", "TACG")


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


def process_one_mutspec(path: str, label=None, rounding=True):
    label = label or "Sample_{:02}".format(random.randint(1, 99))

    mutspec192 = pd.read_csv(PATH_TO_MUTSPEC_192, usecols=COLS)
    mutspec192["NucSubst"] = mutspec192.NucSubst.apply(collapse_mutspec)

    mutspec96 = mutspec192.groupby("NucSubst").ObsToExp.sum().reset_index()
    new_style_mutspec = mutspec96[["NucSubst", "ObsToExp"]]
    new_style_mutspec.columns = ["Mutation Types", label]
    if rounding:
        new_style_mutspec[label] = new_style_mutspec[label].round(0).astype(int)
    return new_style_mutspec


def write_new_style(data: pd.DataFrame, path):
    data.to_csv(path, index=None, sep="\t")


def main():
    m = process_one_mutspec(PATH_TO_MUTSPEC_192, "Sars_all")
    write_new_style(m, PATH_TO_OUT)


if __name__ == "__main__":
    main()
