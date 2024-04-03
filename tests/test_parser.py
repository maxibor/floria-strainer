import os
from floria_strainer.parser import parse_vartigs, parse_vartig_info

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")


def test_parse_vartigs():
    vartigs = parse_vartigs(os.path.join(data_dir, "NZ_CP081897.1.vartigs"))

    assert isinstance(vartigs, dict)
    assert "contig" in vartigs
    assert "haploset" in vartigs
    assert "HAPQ" in vartigs

    assert len(vartigs["contig"]) == len(vartigs["haploset"]) == len(vartigs["HAPQ"])

    assert vartigs["contig"][0] == "NZ_CP081897.1"
    assert vartigs["haploset"][0] == "HAP0"
    assert vartigs["HAPQ"][0] == 16


def test_parse_vartig_info():
    vartig_info = parse_vartig_info(os.path.join(data_dir, "vartig_info_short.txt"))

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
