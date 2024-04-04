import os
import pytest
from floria_strainer.parser import parse_haplosets, parse_vartig_info

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
contig_dir = os.path.join(data_dir, "floria_out_dir", "NZ_CP081897.1")


def test_parse_haplosets():
    haplosets, read_dict = parse_haplosets(
        os.path.join(contig_dir, "NZ_CP081897.1.haplosets"), hapq_cut=15
    )

    # Checking haploset dictionary
    assert isinstance(haplosets, dict)
    assert "contig" in haplosets
    assert "haploset" in haplosets
    assert "HAPQ" in haplosets

    assert (
        len(haplosets["contig"]) == len(haplosets["haploset"]) == len(haplosets["HAPQ"])
    )

    assert haplosets["contig"][0] == "NZ_CP081897.1"
    assert haplosets["haploset"][0] == "HAP0"
    assert haplosets["HAPQ"][0] == 16

    # Checking reads dictionary
    print(haplosets)
    print(read_dict)

    assert isinstance(read_dict, dict)
    assert read_dict["NZ_CP081897.1"]["NZ_CP081897.1-865804"] == "HAP0"
    assert read_dict["NZ_CP081897.1"]["NZ_CP081894.1-245850"] == "HAP10"
    # Checking that haplosets with HAPQ < hapq_cut are not in the read_dict
    with pytest.raises(KeyError):
        assert read_dict["NZ_CP081897.1"]["NZ_CP081896.1-334880"] == "HAP325"


def test_parse_vartig_info():
    vartig_info = parse_vartig_info(os.path.join(contig_dir, "vartig_info_short.txt"))

    assert isinstance(vartig_info, dict)

    keys = set(["contig", "haploset", "pos", "cons", "allele", "support"])

    assert set(vartig_info.keys()) == keys

    for i, key in enumerate(vartig_info):
        if i != 0:
            assert length == len(vartig_info[key])
        length = len(vartig_info[key])

    assert vartig_info["contig"][0] == "NZ_CP081897.1"
    assert vartig_info["haploset"][0] == "HAP0"
    assert vartig_info["pos"][0] == 770
    assert vartig_info["pos"][1] == 1022
    assert vartig_info["pos"][2] == 770
    assert vartig_info["cons"][0] == 0
    assert vartig_info["cons"][4] == 1
    assert vartig_info["cons"][5] == 1
    assert vartig_info["allele"][0] == 0
    assert vartig_info["support"][0] == 35
    assert vartig_info["support"][4] == 9
    assert vartig_info["support"][5] == 4
    assert vartig_info["pos"][6] == 2034
    assert vartig_info["cons"][6] == 0
