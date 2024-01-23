'''
This script was created by Zhang Yujian on January 22nd, 2024.
This script is aimed to convert different sequence infomation involved in bioinformation.
'''

class ConvKit:
    # codon
    dna2ProDict: dict[str, str] = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', \
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', \
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', \
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', \
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', \
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', \
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', \
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', \
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', \
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', \
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', \
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', \
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', \
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', \
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', \
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }


    rna2ProDict: dict[str, str] = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', \
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', \
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', \
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', \
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', \
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', \
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', \
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', \
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', \
        'ACU': 'U', 'ACC': 'U', 'ACA': 'U', 'ACG': 'U', \
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', \
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', \
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', \
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', \
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', \
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }


    # T/U
    dna2RnaDict: dict[str, str] = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'U'}

    rna2DnaDict: dict[str ,str] = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'T'}




    def __init__(self):
        pass
        

# transcription
def dna2Rna(dnaSeq: str) -> str:
    convKit: ConvKit = ConvKit()
    output = dnaSeq.upper().translate(str.maketrans(convKit.dna2RnaDict))
    return(output)


# reverse transcription
def rna2Dna(rnaSeq: str) -> str:
    convKit: ConvKit = ConvKit()
    output = rnaSeq.upper().translate(str.maketrans(convKit.rna2DnaDict))
    return(output)


# DNA translation
def dna2Pro(dnaSeq: str, start: int = 0, end: int = -1) -> str:
    if (end == -1):
        end = len(dnaSeq) - 1
    assert (end > start), 'Invalid arguments!'
    convKit: ConvKit = ConvKit()
    i: int = int(start)
    outputList: list[str] = []
    while (i <= end - 2):
        outputList.append(convKit.dna2ProDict[dnaSeq[i: i+3]])
        i += 3
    output: str = ''.join(outputList)
    return(output)


# RNA translation
def rna2Pro(rnaSeq: str, start: int = 0, end: int = -1) -> str:
    if (end == -1):
        end = len(rnaSeq) - 1
    assert (end > start), 'Invalid arguments!'
    convKit: ConvKit = ConvKit()
    i: int = int(start)
    outputList: list[str] = []
    while (i <= end - 2):
        outputList.append(convKit.rna2ProDict[rnaSeq[i: i+3]])
        i += 3
    output: str = ''.join(outputList)
    return(output)
    