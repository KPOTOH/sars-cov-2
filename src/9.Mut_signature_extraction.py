from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerExtractor import sigpro as sig


def main_function():
    # genInstall.install('GRCh37')

    data = sig.importdata("matrix")
    sig.sigProfilerExtractor(
        "matrix", "../data/example_sig_output", 
        data, minimum_signatures=1, maximum_signatures=3)


if __name__ == "__main__":
    main_function()
