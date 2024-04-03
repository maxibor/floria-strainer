import rich_click as click
from floria_strainer.main import strainer
from floria_strainer import __version__, __author__


@click.command()
@click.version_option(__version__)
@click.argument("floria_outdir", type=click.Path(exists=True))
@click.option(
    "-n",
    "--nb-strains",
    help="Number of strains to keep. If 0, the number of strains will be determined by the mean floria average strain count with HAPQ > 15.",
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
def cli(floria_outdir: str, nb_strains: int, hapq_cut: int, sp_cut: float):
    """
    Strain the haplotypes in the floria output directory.

    Author: Maxime Borry
    """
    strainer(floria_outdir, nb_strains, hapq_cut, sp_cut)
