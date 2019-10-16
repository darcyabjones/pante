#!/usr/bin/env python3

import sys


def print_line(block):
    columns = ["ID", "Name", "Alias", "Ontology_term", "Dbxref", "Note"]
    line = [block.get(c, ".") for c in columns]
    print("\t".join(line))
    return


def main():

    if len(sys.argv) != 2:
        print(f"USAGE: {sys.argv[0]} infile.stk > out.tsv")
        sys.exit(1)

    print("\t".join(["ID", "Name", "Alias",
                     "Ontology_term", "Dbxref", "Note"]))

    try:
        infile_name = sys.argv[1]
        if infile_name == "-":
            infile_handle = sys.stdin
        else:
            infile_handle = open(infile_name, "r")

        block = {"Note": [], "Dbxref": [], "Ontology_term": []}

        for line in infile_handle:
            sline = line.strip()
            if line.startswith("//"):
                if len(block["Dbxref"]) > 0:
                    block["Dbxref"] = ",".join(block["Dbxref"])
                else:
                    block["Dbxref"] = "."

                if len(block["Note"]) > 0:
                    block["Note"] = " ".join(block["Note"]).replace("\t", "%09")
                else:
                    block["Note"] = "."

                if len(block["Ontology_term"]) > 0:
                    block["Ontology_term"] = ",".join(block["Ontology_term"])
                else:
                    block["Ontology_term"] = "."

                print_line(block)

                block = {"Note": [], "Dbxref": [], "Ontology_term": []}

            elif line.startswith("#=GF ID"):
                block["Name"] = sline.split("ID", maxsplit=1)[-1].strip()

            elif line.startswith("#=GF AC"):
                block["ID"] = sline.split("AC", maxsplit=1)[-1].strip()
                if block["ID"].startswith("PF"):
                    block["Dbxref"].append("Pfam:" + block["ID"])

            elif line.startswith("#=GF DE"):
                block["Alias"] = sline.split("DE", maxsplit=1)[-1].strip()

            elif line.startswith("#=GF DR"):
                acc_line = sline.split("DR", maxsplit=1)[-1].strip()
                sacc_line = acc_line.split(";")
                acc = sacc_line[0].strip() + ":" + sacc_line[1].strip()

                if acc.startswith("SO") or acc.startswith("GO"):
                    block["Ontology_term"].append(acc)
                else:
                    block["Dbxref"].append(acc)

            elif line.startswith("#=GF CC"):
                comment = sline.split("CC", maxsplit=1)[-1].strip()
                block["Note"].append(comment)

    finally:
        infile_handle.close()

    return


if __name__ == "__main__":
    main()
