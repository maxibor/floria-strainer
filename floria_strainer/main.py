import pandas as pd
import numpy as np
from typing import List
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
import logging

from floria_strainer.parser import parse_vartigs, parse_vartig_info


def parse_files(vartigs: dict, vartig_info: dict) -> pd.DataFrame:

    vartigs = parse_vartigs(vartigs)
    vartig_info = parse_vartig_info(vartig_info)

    vartigs_df = pd.DataFrame.from_dict(vartigs)
    vartig_info_df = pd.DataFrame.from_dict(vartig_info)

    df = pd.merge(vartigs_df, vartig_info_df, on=["contig", "haploset"])

    return df


def compute_gmm(obs: np.array, n_components: int) -> List(np.array, np.array):
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
    if n_components < 2:
        logging.warning(
            "The number of components you gave less than 2. The optimal number of components will be selected using the Silhouette score."
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

    gmm = GaussianMixture(n_components=n_components)
    labels = gmm.predict(obs)
    max_proba = np.max(gmm.predict_proba(obs), axis=1)
    return labels, max_proba


def process_df(df: pd.DataFrame, hapq_cut: int, nb_strains: int) -> pd.DataFrame:
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
    all_freq["freq"] = all_freq["support"] / all_freq["cov"]
    all_freq["MAF"] = np.where(all_freq["freq"] > 0.5, 1, 0)
    keys = ["contig", "pos", "allele"]
    df2 = all_freq.merge(df[["contig", "pos", "allele", "haploset"]], on=keys)
    obs = df2["freq"].values.reshape(-1, 1)

    labels, max_proba = compute_gmm(obs=obs, n_components=nb_strains)

    df2["strain"] = pd.Categorical(labels)
    df2["strain_proba"] = max_proba
