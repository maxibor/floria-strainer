import rich_click as click
from floria_strainer.main import strainer
from floria_strainer import __version__, __author__


@click.command()
@click.version_option(__version__)
@click.argument("floria_outdir", type=click.Path(exists=True))
@click.option(
    "-n",
    "--nb-strains",
    help="Number of strains to keep. If 0, the number of strains will be determined by the mean floria average strain count with HAPQ > HAPQ treshold",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "-h",
    "--hapq-cut",
    help="Minimum HAPQ threshold",
    type=int,
    default=15,
    show_default=True,
)
@click.option(
    "-s",
    "--sp-cut",
    help="Minimum strain clustering probability threshold",
    type=float,
    default=0.5,
    show_default=True,
)
@click.option(
    "-b", "--bam", help="Input BAM file", type=click.Path(exists=True), required=True
)
@click.option(
    "-m",
    "--mode",
    help="BAM output mode. Tag: add ST (strain) tags to the reads. Split: split the reads in different BAM files per strain.",
    type=click.Choice(["tag", "split"]),
    default="tag",
    show_default=True,
)
@click.option(
    "-o",
    "--basename",
    help="Output file basaneme",
    type=str,
    required=True,
    default="floria_strained",
    show_default=True,
)
def cli(
    floria_outdir: str,
    nb_strains: int,
    hapq_cut: int,
    sp_cut: float,
    bam: str,
    mode: str,
    basename: str,
):
    """
    Strain the haplotypes in the floria output directory.

    Author: Maxime Borry
    """
    strainer(floria_outdir, nb_strains, hapq_cut, sp_cut, bam, mode, basename)
