'''
This script was created on Dec 16th, 2023 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.
Sorry for my poor English.
'''

import os

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

    # PDB file infos import.
    # Column meanings:
    # Serial: Atom serial number.
    # Atom: Atom name.
    # ResName: Residue name.
    # ResSeq: Residue sequence number.
    # ChainId: Chain identifier.
    # X: Orthogonal coordinates for X in angstroms (A).
    # Y: Orthogonal coordinates for Y in angstroms (A).
    # Z: Orthogonal coordinates for Z in angstroms (A).
    pdbInfoColumns: [str] = ['Serial', 'Atom', 'ResName', 'ResSeq', 'ChainId', 'X', 'Y', 'Z']


    def __init__(self) -> None:
        pass



# transfer protein structure file (pdb) into sequence string.
def pdb2Seq(pdbFilePath: str, fastaFilePath: str = None, fastaLineLen: int = 80) -> dict[str, str]:
    proteinKit: ProteinKit = ProteinKit()
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
            if (resName == 'UNK'):
                chainSeq += '-'
            else:
                chainSeq += proteinKit.aaDictTHREE2One[resName]

        while (line):
            line = pdbFile.readline()
            if (line.startswith('ATOM')):
                if (resSeq == 0):                                       # chain begin
                    thisChainId = line[21]
                    resSeq = int(line[22:26].strip())
                    resName = line[17:20]
                    if (resName == 'UNK'):
                        chainSeq += '-'
                    else:
                        chainSeq += proteinKit.aaDictTHREE2One[resName]
                elif (int(line[22:26].strip()) == resSeq):
                    continue
                else:
                    if (int(line[22:26].strip()) == resSeq + 1):
                        resSeq = int(line[22:26].strip())
                        resName = line[17:20]
                        if (resName == 'UNK'):
                            chainSeq += '-'
                        else:
                            chainSeq += proteinKit.aaDictTHREE2One[resName]
                    else:
                        gap: int = int(line[22:26].strip()) - resSeq - 1
                        gapSeq: str = '-' * gap
                        chainSeq += gapSeq
                        resSeq = int(line[22:26].strip())
                        resName = line[17:20]
                        if (resName == 'UNK'):
                            chainSeq += '-'
                        else:
                            chainSeq += proteinKit.aaDictTHREE2One[resName]
            elif (line.startswith('TER')):                              # recognizing termination with the line 'TER'
                output[thisChainId] = chainSeq
                resSeq = 0
                chainSeq = ''
    # output the fasta files.
    if (fastaFilePath is not None):
        fileName: str = pdbFilePath.split(os.sep)[-1].rstrip('.pdb')
        with open(file=os.path.join(fastaFilePath, fileName+'.fasta'), mode='w') as fastaFile:
            for key in output.keys():
                fastaFile.write('>'+fileName+'_chain_'+key+'\n')
                thisLine: list[str] = [output[key][i:i+fastaLineLen] for i in range(0, len(output[key]), fastaLineLen)]
                for i in thisLine:
                    fastaFile.write(i + '\n')
            print(fileName+' converted!\n')
    return (output)



# load the information of all amino-acid-residue atoms into a list. output[idx] shows the information of an atom.
def pdb2List(pdbFilePath: str, csvPath: str = None, colName: bool = False) -> list[list[int, str, str, int, str, float, float, float]]:
    output: list[list] = []
    # pdbInfoColumns: [str] = ['Serial', 'Atom', 'ResName', 'ResSeq', 'ChainId', 'X', 'Y', 'Z']
    if (colName):
        proteinKit: ProteinKit = ProteinKit()
        output.append(proteinKit.pdbInfoColumns)
    with open(file=pdbFilePath, mode='r') as pdbFile:
        line: str = pdbFile.readline()
        if (line.startswith('ATOM')):
            output.append([int(line[6:11].strip()), str(line[12:16].strip()), \
                            str(line[17:20].strip()), int(line[22:26].strip()), \
                            str(line[21]), float(line[30:38].strip()), \
                            float(line[38:46].strip()), float(line[46:54])])
        while (line):
            line = pdbFile.readline()
            if (line.startswith('ATOM')):
                output.append([int(line[6:11].strip()), str(line[12:16].strip()), \
                            str(line[17:20].strip()), int(line[22:26].strip()), \
                            str(line[21]), float(line[30:38].strip()), \
                            float(line[38:46].strip()), float(line[46:54])])
    # output the csv file.
    if (csvPath is not None):
        for atomListIdx in range(len(output)):
            output[atomListIdx] = list(map(str, output[atomListIdx]))
        with open(file=csvPath, mode='w') as csvFile:
            csvFile.writelines([','.join(atomInfoList)+'\n' for atomInfoList in output])
    return (output)



# load the information of all amimo-acid-residue atoms into a dictionary, which can be converted into a dataframe with famous `pandas`.
def pdb2Dict(pdbFilePath: str) -> dict[str, list]:
    # pdbInfoColumns: [str] = ['Serial', 'Atom', 'ResName', 'ResSeq', 'ChainId', 'X', 'Y', 'Z']
    proteinKit: ProteinKit = ProteinKit()
    keyNameList: list[str] = proteinKit.pdbInfoColumns
    output: dict[str, list] = {keyName: [] for keyName in keyNameList}
    with open(file=pdbFilePath, mode='r') as pdbFile:
        line: str = pdbFile.readline()
        if (line.startswith('ATOM')):
            output[keyNameList[0]].append(int(line[6:11].strip()))
            output[keyNameList[1]].append(str(line[12:16].strip()))
            output[keyNameList[2]].append(str(line[17:20].strip()))
            output[keyNameList[3]].append(int(line[22:26].strip()))
            output[keyNameList[4]].append(str(line[21]))
            output[keyNameList[5]].append(float(line[30:38].strip()))
            output[keyNameList[6]].append(float(line[38:46].strip()))
            output[keyNameList[7]].append(float(line[46:54].strip()))
        while (line):
            line = pdbFile.readline()
            if (line.startswith('ATOM')):
                output[keyNameList[0]].append(int(line[6:11].strip()))
                output[keyNameList[1]].append(str(line[12:16].strip()))
                output[keyNameList[2]].append(str(line[17:20].strip()))
                output[keyNameList[3]].append(int(line[22:26].strip()))
                output[keyNameList[4]].append(str(line[21]))
                output[keyNameList[5]].append(float(line[30:38].strip()))
                output[keyNameList[6]].append(float(line[38:46].strip()))
                output[keyNameList[7]].append(float(line[46:54].strip()))
    return (output)

