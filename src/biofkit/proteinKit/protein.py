'''
This script was created on Jan 3rd, 2024 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.
Sorry for my poor English.
'''

import os
import pdbKit

# This funtion 'proDict2ProList' is used to tranfer pdbDict into pdbList.
def proDict2ProList(rawDict: dict) -> list:
    output: list = [[] for i in range(len(rawDict['Serial']))]
    for atomIdx in range(len(rawDict['Serial'])):
        output[atomIdx].append(rawDict['Serial'][atomIdx])
        output[atomIdx].append(rawDict['Atom'][atomIdx])
        output[atomIdx].append(rawDict['ResName'][atomIdx])
        output[atomIdx].append(rawDict['ResSeq'][atomIdx])
        output[atomIdx].append(rawDict['ChainId'][atomIdx])
        output[atomIdx].append(rawDict['X'][atomIdx])
        output[atomIdx].append(rawDict['Y'][atomIdx])
        output[atomIdx].append(rawDict['Z'][atomIdx])
    return (output)


# This function 'proList2ProDIct' is used to tranfer pdbList into pdbDict.
def proList2ProDict(rawList: list) -> dict:
    output: dict = {'Serial': [], 'Atom': [], 'ResName': [], 'ResSeq': [], \
                    'ChainId': [], 'X': [], 'Y': [], 'Z': []}
    for atomIdx in range(len(rawList)):
        output['Serial'].append(rawList[atomIdx][0])
        output['Atom'].append(rawList[atomIdx][1])
        output['ResName'].append(rawList[atomIdx][2])
        output['ResSeq'].append(rawList[atomIdx][3])
        output['ChainId'].append(rawList[atomIdx][4])
        output['X'].append(rawList[atomIdx][5])
        output['Y'].append(rawList[atomIdx][6])
        output['Z'].append(rawList[atomIdx][7])
    return (output)

# This new type 'protein' is for operation pdb file.
class Protein:
    def __init__(self, pdbFilePath: str, chainId: list[str] = ['all']):
        self.name: str = pdbFilePath.split(os.sep)[-1].rstrip('.pdb')
        self.infoList: list = pdbKit.pdb2List(pdbFilePath)
        self.infoDict: dict = pdbKit.pdb2Dict(pdbFilePath)
        if (chainId != ['all']):
            for atom in self.info:
                if (atom[4] not in chainId):
                    self.info.remove(atom)

    def __str__(self) -> str:
        return ('{name}: containing Chain {chainList}'.format(name = self.name, chainList = list(set(self.infoDict['ChainId']))))
    
    def __repr__(self) -> str:
        return ('{name}: containing Chain {chainList}'.format(name = self.name, chainList = list(set(self.infoDict['ChainId']))))
    
