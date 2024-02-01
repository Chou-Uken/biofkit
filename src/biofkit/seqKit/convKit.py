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
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', \
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
        outputList.append(convKit.dna2ProDict[dnaSeq[i: i+3].upper()])
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
        outputList.append(convKit.rna2ProDict[rnaSeq[i: i+3].upper()])
        i += 3
    output: str = ''.join(outputList)
    return(output)
    

# Alignment
def pairwiseAlign(seqA: str, seqB: str, method: str = 'NW', consoleWidth = 50):
    assert (len(seqA) * len(seqB) != 0), 'empty string'
    scoreMatrix: list[list[int]] = []
    for i in range(len(seqB)+1):
        scoreMatrix.append([])
        for j in range(len(seqA)+1):
            scoreMatrix[i].append(0)
    # Needleman-Wunsch (NW)
    # Matrix initiation
    for rowIdx in range(len(seqB)+1):
        scoreMatrix[rowIdx][0] = -rowIdx
    for columnIdx in range(len(seqA)+1):
        scoreMatrix[0][columnIdx] = -columnIdx
    # Matrix scan
    for rowIdx in range(1, len(seqB)+1):
        for columnIdx in range(1, len(seqA)+1):
            if (seqB[rowIdx-1] == seqA[columnIdx-1]):
                scoreMatrix[rowIdx][columnIdx] = max(scoreMatrix[rowIdx-1][columnIdx]-1, \
                                                     scoreMatrix[rowIdx-1][columnIdx]-1, \
                                                        scoreMatrix[rowIdx-1][columnIdx-1]+1)
                
            else:
                scoreMatrix[rowIdx][columnIdx] = max(scoreMatrix[rowIdx-1][columnIdx]-1, \
                                                        scoreMatrix[rowIdx-1][columnIdx]-1, \
                                                        scoreMatrix[rowIdx-1][columnIdx-1]-1)
                                                    
    finalScore = scoreMatrix[-1][-1]
    # backforward
    steps: list[str] = []
    rowIdx = len(seqB)
    columnIdx = len(seqA)
    while ((rowIdx != 0) or (columnIdx != 0)):

        if ((rowIdx != 0) and (columnIdx != 0)):
            bestPathScore = max(scoreMatrix[rowIdx-1][columnIdx], \
                                scoreMatrix[rowIdx][columnIdx-1], \
                                    scoreMatrix[rowIdx-1][columnIdx-1])
            # up-leftward
            if (bestPathScore == scoreMatrix[rowIdx-1][columnIdx-1]):
                steps.insert(0, 'c')
                rowIdx -= 1
                columnIdx -= 1
            # upward
            elif (bestPathScore == scoreMatrix[rowIdx-1][columnIdx]):
                steps.insert(0, 'd')
                rowIdx -= 1
            # leftward
            else:
                steps.insert(0, 'r')
                columnIdx -= 1
            
        elif ((rowIdx == 0) and (columnIdx != 0)):
            steps.insert(0, 'r')
            columnIdx -= 1
            
        else:
            steps.insert(0, 'd')
            rowIdx -= 1

    seqAIdx: int = 0
    seqBIdx: int = 0
    seqAOut: str = ''
    seqBOut: str = ''
    for step in steps:
        if (step == 'r'):
            seqAOut += seqA[seqAIdx]
            seqBOut += '-'
            seqAIdx += 1
        elif (step == 'd'):
            seqAOut += '-'
            seqBOut += seqB[seqBIdx]
            seqBIdx += 1
        else:
            seqAOut += seqA[seqAIdx]
            seqBOut += seqB[seqBIdx]
            seqBIdx += 1
            seqAIdx += 1


    seqAlign: str = ''
    for idx in range(len(seqAOut)):
        if (seqAOut[idx] == seqBOut[idx]):
            seqAlign += '|'
        else:
            seqAlign += ' '
    
    
    seqAOut = [seqAOut[i:i+consoleWidth] for i in range(0, len(seqAOut), consoleWidth)]
    seqAlign = [seqAlign[i:i+consoleWidth] for i in range(0, len(seqAlign), consoleWidth)]
    seqBOut = [seqBOut[i:i+consoleWidth] for i in range(0, len(seqBOut), consoleWidth)]
    for idx in range(len(seqAOut)):
        print(seqAOut[idx])
        print(seqAlign[idx])
        print(seqBOut[idx])
        print()

        
if __name__ == '__main__':
    pairwiseAlign('TGTTTGAACGTTTACAGACTAAACTTCACCTGAAATCCTCCCAGCAGAGAGCAAAGGTGGTGCCTCCCTCCCTACAAAACCCCCGTCTGTCTGCAGATTAACCTTTCTCTGGACGGACGGACGGCAGGTGAAGGACGGAG', 'GGCACCATGGCAACCGCTGCAGATCAGAACGTGGAGTTTGTTAGAACCGGCTACGGGAAGAACTCGGTGAAGGTTCTGTTCATCCGGAGGCAGAGGAACCACCACGAGATCATCGAGCTGAAGGCCGACGTGGAGCTGAC')
