from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
from dataclasses import dataclass


@dataclass
class Protein:
    id: str = ""
    name: str = ""
    organism: str = ""
    description: str = ""
    sequence: str = ""


def read_proteins(file_in: str) -> list[Protein]:
    file = open(file_in)
    records = SeqIO.parse(file, "fasta")

    proteins = []
    # filter records with a word "Fragment"
    stop_word = "Fragment"
    for record in records:
        if stop_word not in record.description:
            protein = Protein()
            protein.id = record.id
            protein.name = record.name
            protein.sequence = record.seq
            protein.description = record.description
            protein.organism = get_taxonomy_from_description(record.description)
            proteins.append(protein)

    proteins = sort_proteins_by_taxonomy(proteins)
    print(len(proteins))
    return proteins


# get list with all organism's names
def get_taxa(proteins: list[Protein]) -> list[str]:
    result = [p.organism for p in proteins]
    return list(set(result))


def print_taxa(taxa: list[str]) -> None:
    for t in taxa:
        print(t)


def sort_proteins_by_taxonomy(proteins:list[Protein]) -> list[Protein]:
    def sort_func(protein: Protein) -> str:
        return protein.organism
    result = sorted(proteins, key=sort_func)
    return result


def get_taxonomy_from_description(description:str) -> str:
    # description format:
    # tr|A0A663N598|A0A663N598_ATHCN Hexokinase OS=Athene cunicularia OX=194338 GN=LOC113481846 PE=3 SV=1
    # organism name:
    # OS=Athene cunicularia
    os = ""
    pattern = r"OS=\w+\s\w+"
    match = re.search(pattern, description)
    if match:
        os = match[0][3:]
    else:
        print("wrong os: ", description)
    return os


def write_proteins(proteins: list[Protein], file: str, file_format: str="fasta") -> None:
    def convert_to_record(protein: Protein) -> SeqRecord:
        # change the protein's id from format
        # "tr|A0A663N598|A0A663N598" to "A0A663N598"
        # change description from format
        # > tr|A0A663N598|A0A663N598_ATHCN Hexokinase OS=Athene cunicularia OX=194338 GN=LOC113481846 PE=3 SV=1
        # to format
        # > A0A663N598 Athene cunicularia
        ids = protein.id.split("|")
        record = SeqRecord(
            id=ids[1],
            seq=protein.sequence,
            name=protein.name,
            description=protein.organism  # + " " + ids[1]#protein.description
        )
        return record

    # convert Protein to SeqRecord
    sequences = [convert_to_record(p) for p in proteins]
    with open(file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, file_format)




def init() -> None:
    hex1 = "B20211126A084FC58F6BBA219896F365D15F2EB44020D0BZ.fasta"
    hex2 = "B20211126A084FC58F6BBA219896F365D15F2EB44020F41H.fasta"
    proteins_hex1 = read_proteins(hex1)
    proteins_hex2 = read_proteins(hex2)

    write_proteins(proteins_hex1, "hex1.fasta")
    write_proteins(proteins_hex2, "hex2.fasta")

    taxa_hex1 = get_taxa(proteins_hex1)
    taxa_hex2 = get_taxa(proteins_hex2)

    common_taxa = set(taxa_hex1).intersection(set(taxa_hex2))
    taxa_only_hex1 = set(taxa_hex1).difference(common_taxa)
    taxa_only_hex2 = set(taxa_hex2).difference(common_taxa)
    print(taxa_only_hex1)
    print(taxa_only_hex2)


if __name__ == '__main__':
    init()

