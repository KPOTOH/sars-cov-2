from SigProfilerExtractor import sigpro as sig
import click


@click.command("extractor", help="Pass path to mutspec table in SigProfiler format and path to output directory. Script release signatures from mutspec and decompose it with COSMIC database")
@click.argument("input_mutspec", required=True, type=click.Path(exists=True))
@click.argument("outdir", required=True, type=click.Path(exists=False))
@click.option("-m", "--max_signatures", default=5, show_default=True, type=int, help="maximum number of signatures to release")
@click.option("-t", "--threads", default=-1, show_default=True, type=int, help="number of threads to use")
def main(input_mutspec, outdir, max_signatures, threads):
    sig.sigProfilerExtractor(
        "matrix", outdir,
        input_mutspec,
        minimum_signatures=1,
        maximum_signatures=max_signatures,
        cpu=threads,
        gpu=False,
    )


if __name__ == "__main__":
    main()
