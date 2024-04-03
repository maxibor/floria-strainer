import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from collections import ChainMap
import logging
import pysam
import os

from floria_strainer.parser import (
    parse_vartigs,
    parse_vartig_info,
    parse_floria_contig_ploidy,
)


def parse_files(floria_outdir: str) -> tuple[pd.DataFrame, int]:

    cp = parse_floria_contig_ploidy(
        os.path.join(floria_outdir, "contig_ploidy_info.tsv")
    )

    contigs = cp["contig"].unique().tolist()
    nb_strains = round(cp["average_straincount_min15hapq"].mean())

    vartigs_files = [
        os.path.join(floria_outdir, contig, f"{contig}.vartigs") for contig in contigs
    ]
    vartig_info_files = [
        os.path.join(floria_outdir, contig, "vartig_info.txt") for contig in contigs
    ]

    parsed_vartigs = [parse_vartigs(vartig_file) for vartig_file in vartigs_files]

    parsed_vartig_info = [
        parse_vartig_info(vartig_info_file) for vartig_info_file in vartig_info_files
    ]

    vartigs_df = pd.DataFrame.from_dict(dict(ChainMap(*parsed_vartigs)))
    vartig_info_df = pd.DataFrame.from_dict(dict(ChainMap(*parsed_vartig_info)))

    df = pd.merge(vartigs_df, vartig_info_df, on=["contig", "haploset"])

    return df, nb_strains


def compute_gmm(obs: np.array, n_components: int) -> tuple[np.array, np.array]:
    """
    Compute the Gaussian Mixture Model for the given observations.

    Parameters
    ----------
    obs : np.array
        The observations to compute the GMM for.
    n_components : int
        The number of components to use in the GMM. If not specified, the
        number of components will be selected using the Silhouette score.

    Returns
    -------
    List(np.array, np.array)
        The labels and the maximum probability for each observation.
    """
    if n_components == 0:
        logging.warning(
            "The number of components is 0. The optimal number of components will be selected using the Silhouette score."
        )

        best_comp = 2
        old_score = 0
        for n in range(2, 5):
            gmm = GaussianMixture(n_components=n)
            gmm.fit(obs)
            labels = gmm.predict(obs)
            score = silhouette_score(obs, labels)
            if score > old_score:
                best_comp = n
                old_score = score
        n_components = best_comp

    elif n_components < 2:
        raise ValueError("The number of components must be at least 2.")
    gmm = GaussianMixture(n_components=n_components)
    gmm.fit(obs)
    labels = gmm.predict(obs)
    max_proba = np.max(gmm.predict_proba(obs), axis=1)
    return labels, max_proba


def process_df(df: pd.DataFrame, hapq_cut: int, sp_cut: float, nb_strains: int) -> dict:
    """
    Process the DataFrame to get the strains.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to process.
    hapq_cut : int
        The hapq cut-off to use.
    sp_cut : float
        The clustering probability support cut-off to use.
    nb_strains : int
        The number of strains to use. Use 0 to automatically select the best number of strains.

    Returns
    -------
    dict: {str: {str: int}}
        {
            contig(str) : {
                haploset(str) : strain(int)
            }
        }
    list: [int]
        List of strains
    """

    df = df[df["HAPQ"] >= hapq_cut]
    support = df.groupby(["contig", "pos", "allele"])["support"].sum().reset_index()
    coverage = (
        support.groupby(["contig", "pos"])["support"]
        .sum()
        .reset_index()
        .rename(columns={"support": "coverage"})
    )
    all_freq = support.merge(
        coverage, left_on=["contig", "pos"], right_on=["contig", "pos"]
    )

    all_freq["freq"] = all_freq["support"] / all_freq["coverage"]
    all_freq["MAF"] = np.where(all_freq["freq"] > 0.5, 1, 0)
    keys = ["contig", "pos", "allele"]
    df2 = all_freq.merge(df[["contig", "pos", "allele", "haploset"]], on=keys)

    obs = df2["freq"].values.reshape(-1, 1)

    labels, max_proba = compute_gmm(obs=obs, n_components=nb_strains)

    df2["strain"] = pd.Categorical(labels)
    df2["strain_proba"] = max_proba

    df2 = df2[df2["strain_proba"] >= sp_cut]
    strains = df2["strain"].unique().tolist()
    df2["strain"] = df2["strain"].astype(int)

    df3 = (
        df2.groupby(["contig", "haploset"])["strain"]
        .mean()
        .reset_index()
        .assign(strain=lambda df: df["strain"].apply(round))
    )

    return (
        df3.set_index("haploset")
        .groupby("contig")
        .apply(lambda x: x["strain"].to_dict())
        .to_dict(),
        strains,
    )


def write_bam_split(
    inbam: str, basename: str, haplostrain: dict, strains: list
) -> None:
    """
    Write the BAM files for each strain.

    Parameters
    ----------
    inbam : str
        The input BAM file.
    basename : str
        The basename to use for the output BAM files.
    haplostrain : dict
        The haplostrain dictionary.
    strains : list
        The strains to write the BAM files for.
    """

    bam = pysam.AlignmentFile(inbam, "rb")
    for strain in strains:
        with pysam.AlignmentFile(
            f"{basename}.{strain}.bam", "wb", template=bam
        ) as outbam:
            for read in bam:
                refname = read.reference_name
                if refname in haplostrain:
                    try:
                        h = read.get_tag("HP")
                        if haplostrain[refname][h] == strain:
                            outbam.write(read)
                    except KeyError:
                        pass
    bam.close()


def write_bam(inbam: str, outbam: str, haplostrain: dict) -> None:
    """
    Write the BAM files for each strain.

    Parameters
    ----------
    inbam : str
        The input BAM file.
    outbam : str
        The output BAM file.
    haplostrain : dict
        The haplostrain dictionary.
    """

    with pysam.AlignmentFile(inbam, "rb") as bam:
        with pysam.AlignmentFile(outbam, "wb", template=bam) as bamout:
            for read in bam:
                refname = read.reference_name
                if refname in haplostrain:
                    try:
                        h = read.get_tag("HP")
                        read.set_tag("ST", haplostrain[refname][h])
                    except KeyError:
                        pass
                bamout.write(read)


def strainer(floria_outdir, nb_strains: int, hapq_cut: int, sp_cut: float):
    fl_df, fl_nb_strains = parse_files(floria_outdir)
    if nb_strains == 0:
        nb_strains = fl_nb_strains
    fl_df_processed, strains = process_df(
        df=fl_df, hapq_cut=hapq_cut, sp_cut=sp_cut, nb_strains=nb_strains
    )

    logging.info(f"Strains: {strains}")
