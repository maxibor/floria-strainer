def parse_vartigs(filename: str) -> dict:
    """
    Read a vartigs file and return a dictionary with the values.

    Parameters
    ----------
    filename : str
        The name of the file to read.

    Returns
    -------
    dict
        A dictionary with the values in the file.
    """
    vartigs = {"contig": [], "haploset": [], "HAPQ": []}
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                lsplit = line.rstrip()[1:].split("\t")
                vartigs["contig"].append(lsplit[1].split(":")[1])
                vartigs["haploset"].append(lsplit[0].split(".")[0])
                vartigs["HAPQ"].append(int(lsplit[6].split(":")[1]))
            else:
                continue
    return vartigs


def parse_vartig_info(filename: str) -> dict:
    """
    Read a vartig_info.txt file and return a dictionary with the values.

    Parameters
    ----------
    filename : str
        The name of the file to read.

    Returns
    -------
    dict
        A dictionary with the values in the file.
    """

    vartig_all = {
        "contig": [],
        "haploset": [],
        "pos": [],
        "cons": [],
        "allele": [],
        "support": [],
    }
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                lsplit = line.rstrip()[1:].split("\t")
                hpset = lsplit[0].split(".")[0]
                contig = lsplit[0].split("/")[-1]
                continue
            else:
                try:
                    lsplit = line.split("\t")
                    pos = int(lsplit[0].split(":")[1])
                    cons = lsplit[1]
                    als = lsplit[2].split("|")
                    if cons != "?":
                        for all in als:
                            al, sup = all.split(":")
                            vartig_all["contig"].append(contig)
                            vartig_all["haploset"].append(hpset)
                            vartig_all["pos"].append(int(pos))
                            vartig_all["cons"].append(int(cons))
                            vartig_all["allele"].append(int(al))
                            vartig_all["support"].append(int(sup))
                except ValueError:
                    continue
                continue

    return vartig_all
