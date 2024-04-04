import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from collections import ChainMap
from floria_strainer import logger
import pysam
import os

from floria_strainer.parser import (
    parse_haplosets,
    parse_vartig_info,
    parse_floria_contig_ploidy,
)


def parse_files(floria_outdir: str) -> tuple[pd.DataFrame, int, dict]:

    cp = parse_floria_contig_ploidy(
        os.path.join(floria_outdir, "contig_ploidy_info.tsv")
    )

    contigs = cp["contig"].unique().tolist()
    nb_strains = round(cp["average_straincount_min15hapq"].mean())

    haplosets_files = [
        os.path.join(floria_outdir, contig, f"{contig}.haplosets") for contig in contigs
    ]
    vartig_info_files = [
        os.path.join(floria_outdir, contig, "vartig_info.txt") for contig in contigs
    ]

    parsed_haplosets = [
        parse_haplosets(haplosets_file) for haplosets_file in haplosets_files
    ]

    haplosets, read_dicts = zip(*parsed_haplosets)

    parsed_vartig_info = [
        parse_vartig_info(vartig_info_file) for vartig_info_file in vartig_info_files
    ]

    haplosets_df = pd.DataFrame.from_dict(dict(ChainMap(*haplosets)))
    vartig_info_df = pd.DataFrame.from_dict(dict(ChainMap(*parsed_vartig_info)))
    read_dict = dict(ChainMap(*read_dicts))

    df = pd.merge(haplosets_df, vartig_info_df, on=["contig", "haploset"])

    return df, nb_strains, read_dict


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
        logger.warning(
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
    gmm = GaussianMixture(n_components=n_components, random_state=42)
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
    inbam: str, basename: str, haplostrain: dict, strains: list, read_dict: dict
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
    read_dict : dict
        The read dictionary with the references as keys for the reads as keys, and the haplotag as value.
    """
    no_haploset = False
    for strain in strains:
        logger.info(
            f"Writing the BAM files in split mode for strain {strain} to {basename}.{strain}.bam"
        )
        bam = pysam.AlignmentFile(inbam, "rb")
        with pysam.AlignmentFile(
            f"{basename}.{strain}.bam", "wb", template=bam
        ) as outbam:
            for read in bam:
                refname = read.reference_name
                if refname in read_dict:
                    try:
                        haplotag = read_dict[refname][read.query_name]
                        read.set_tag("HP", haplotag)
                        no_haploset = False
                    except KeyError:
                        no_haploset = True
                        pass
                if refname in haplostrain and not no_haploset:
                    try:
                        if haplostrain[refname][haplotag] == strain:
                            read.set_tag("ST", haplostrain[refname][haplotag])
                            outbam.write(read)
                    except KeyError:
                        pass
                elif no_haploset:
                    outbam.write(read)
    bam.close()


def write_bam(inbam: str, outbam: str, haplostrain: dict, read_dict: dict) -> None:
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
    read_dict : dict
        The read dictionary with the references as keys for the reads as keys, and the haplotag as value.
    """
    logger.info(f"Writing the BAM file in tag mode to {outbam}.")
    no_haploset = False
    with pysam.AlignmentFile(inbam, "rb") as bam:
        with pysam.AlignmentFile(outbam, "wb", template=bam) as outbam:
            for read in bam:
                refname = read.reference_name
                if refname in read_dict:
                    try:
                        haplotag = read_dict[refname][read.query_name]
                        read.set_tag("HP", haplotag)
                        no_haploset = False
                    except KeyError:
                        no_haploset = True
                        pass
                if refname in haplostrain and not no_haploset:
                    try:
                        read.set_tag("ST", haplostrain[refname][haplotag])
                    except KeyError:
                        pass
                outbam.write(read)


def strainer(
    floria_outdir,
    nb_strains: int,
    hapq_cut: int,
    sp_cut: float,
    bam: str,
    mode: str,
    basename: str,
):
    fl_df, fl_nb_strains, read_dict = parse_files(floria_outdir)
    if nb_strains == 0:
        nb_strains = fl_nb_strains
    fl_processed, strains = process_df(
        df=fl_df, hapq_cut=hapq_cut, sp_cut=sp_cut, nb_strains=nb_strains
    )

    logger.info(
        f"{len(strains)} strains were found: {', '.join([str(i) for i in strains])}"
    )
    if mode == "tag":
        write_bam(
            inbam=bam,
            outbam=f"{basename}.bam",
            haplostrain=fl_processed,
            read_dict=read_dict,
        )
    elif mode == "split":
        write_bam_split(
            inbam=bam,
            basename=basename,
            haplostrain=fl_processed,
            strains=strains,
            read_dict=read_dict,
        )
    else:
        raise ValueError("The mode must be either 'tag' or 'split'.")
