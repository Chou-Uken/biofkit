'''
This script was created on Dec 16th, 2023 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.
Sorry for my poor English.
'''

import os
import typing

class ProteinKit():

    # Dictionary used to transfer abbreviation with the first letter capitalized to shorter abbreviation.
    aaDictTHREE2One: dict[str, str] = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',    \
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',    \
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',    \
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',    \
        'SEC': 'U', 'PYL': 'O'                                          # Rare amino acid
    }

    aaDictThree2One: dict[str, str] = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',    \
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',    \
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',    \
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',    \
        'Sec': 'U', 'Pyl': 'O'                                          # Rare amino acid
    }

    # Dictionary used to transfer abbreviation to shorter abbreviation.
    aaDictthree2One: dict[str, str] = {
        'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'cys': 'C',    \
        'gln': 'Q', 'glu': 'E', 'gly': 'G', 'his': 'H', 'ile': 'I',    \
        'leu': 'L', 'lys': 'K', 'met': 'M', 'phe': 'F', 'pro': 'P',    \
        'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V',    \
        'sec': 'U', 'pyl': 'O'                                          # Rare amino acid
    }

    # Dictionary used to transfer shorter abbreviation to abbreviation with the first letter capitalized.
    aaDictOne2Three: dict[str, str] = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',    \
        'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',    \
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',    \
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',    \
        'U': 'Sec', 'O': 'Pyl'                                          # Rare amino acid
    }

    # Dictionary used to transfer shorter abbreviation to abbreviation.
    aaDictOne2three: dict[str, str] = {
        'A': 'ala', 'R': 'arg', 'N': 'asn', 'D': 'asp', 'C': 'cys',    \
        'Q': 'gln', 'E': 'glu', 'G': 'gly', 'H': 'his', 'I': 'ile',    \
        'L': 'leu', 'K': 'lys', 'M': 'met', 'F': 'phe', 'P': 'pro',    \
        'S': 'ser', 'T': 'thr', 'W': 'trp', 'Y': 'tyr', 'V': 'val',    \
        'U': 'sec', 'O': 'pyl'                                          # Rare amino acid
    }


    def __init__(self) -> None:
        pass



# transfer protein structure file (pdb) into sequence string.
def pdb2Seq(pdbFilePath: str, fasta: bool = False) -> dict[str, str]:
    with open(file=pdbFilePath) as pdbFile:
        thisChainId: str = 'defined'
        line: str = pdbFile.readline()
        chainSeq: str = ''
        resSeq: int = 0
        resName: str = ''
        output: dict[str, str] = {}

        if (line.startswith('ATOM')):
            thisChainId = line[21]
            resSeq = int(line[22:26].strip())
            resName = line[17:20]
            chainSeq += ProteinKit.aaDictTHREE2One[resName]

        while (line):
            line = pdbFile.readline()
            if (line.startswith('ATOM')):
                if (resSeq == 0):                                       # chain begin
                    thisChainId = line[21]
                    resSeq = int(line[22:26].strip())
                    resName = line[17:20]
                    chainSeq += ProteinKit.aaDictTHREE2One[resName]
                elif (int(line[22:26].strip()) == resSeq):
                    continue
                else:
                    assert (int(line[22:26].strip()) == resSeq + 1), 'wrong index'
                    resSeq = int(line[22:26].strip())
                    resName = line[17:20]
                    chainSeq += ProteinKit.aaDictTHREE2One[resName]
            elif (line.startswith('TER')):                              # recognizing termination with the line 'TER'
                output[thisChainId] = chainSeq
                resSeq = 0
                chainSeq = ''
    if (fasta):
        fileName: str = pdbFilePath.split(os.sep)[-1].rstrip('.pdb')
        with open(file=fileName+'.fasta', mode='w') as fastaFile:
            for key in output.keys():
                fastaFile.write('>'+fileName+'_chain_'+key+'\n')
                fastaFile.write(output[key]+'\n')
    return (output)


