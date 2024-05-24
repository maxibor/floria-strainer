import os
import pytest
import pandas as pd
from floria_strainer.main import parse_files, process_df

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
floria_out_dir = os.path.join(data_dir, "floria_out_dir")


@pytest.fixture(scope="module")
def parsed_files():
    df, nb_strains, read_dict = parse_files(floria_out_dir, hapq_cut=15)
    return df, nb_strains, read_dict


def test_parse_files(parsed_files):
    df, nb_strains, read_dict = parsed_files

    assert isinstance(df, pd.DataFrame)
    assert isinstance(nb_strains, int)
    assert isinstance(read_dict, dict)

    assert "contig" in df.columns
    assert "haploset" in df.columns
    assert "HAPQ" in df.columns
    assert "pos" in df.columns
    assert "cons" in df.columns
    assert "allele" in df.columns
    assert "support" in df.columns

    assert df.shape == (2522, 7)
    assert nb_strains == 2
    assert "NZ_CP081897.1" in read_dict
    assert len(read_dict["NZ_CP081897.1"]) == 21579
    assert "NZ_CP081897.1-865804" in read_dict["NZ_CP081897.1"]
    assert "NZ_CP081894.1-245850" in read_dict["NZ_CP081897.1"]
    assert "NZ_CP081896.1-334880" not in read_dict["NZ_CP081897.1"]


def test_process_df(parsed_files):
    df, nb_strains, read_dict = parsed_files

    dd, strains = process_df(df, hapq_cut=15, sp_cut=0.9, nb_strains=nb_strains, basename="test")

    assert isinstance(strains, list)
    assert len(strains) == 2
    assert isinstance(dd, dict)
    assert isinstance(dd["NZ_CP081897.1"], dict)
    assert len(dd["NZ_CP081897.1"]) == 118
    assert "HAP0" in dd["NZ_CP081897.1"]
    assert "HAP1" in dd["NZ_CP081897.1"]
    assert dd["NZ_CP081897.1"]["HAP0"] == 0
    assert dd["NZ_CP081897.1"]["HAP1"] == 1
