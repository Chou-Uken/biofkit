'''
This script was created on Jan 3rd, 2024 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.
Sorry for my poor English.
'''

import os
import pdbKit

# This funtion 'proDict2ProList' is used to tranfer pdbDict into pdbList.
def proDict2ProList(rawDict: dict) -> list:
    output: list = [[] for count in range(len(rawDict['Serial']))]
    for atomIdx in range(len(rawDict['Serial'])):
        output[atomIdx].append(rawDict['Serial'][atomIdx])
        output[atomIdx].append(rawDict['Atom'][atomIdx])
        output[atomIdx].append(rawDict['ResName'][atomIdx])
        output[atomIdx].append(rawDict['ResSeq'][atomIdx])
        output[atomIdx].append(rawDict['ChainId'][atomIdx])
        output[atomIdx].append(rawDict['X'][atomIdx])
        output[atomIdx].append(rawDict['Y'][atomIdx])
        output[atomIdx].append(rawDict['Z'][atomIdx])


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


# Define proList is valid or not?
def proListIsValid(proList: list) -> bool:
    # form and type test
    try:
        for atomInfo in proList:
            # form test:
            if (len(atomInfo != 8)):
                raise (IndexError)
            # type test:
            if ((type(atomInfo[0]) != int) or (type(atomInfo[1]) != str) or \
                (type(atomInfo[2]) != str) or (type(atomInfo[3]) != int) or \
                (type(atomInfo[4]) != str) or (type(atomInfo[5]) != float) or \
                (type(atomInfo[6]) != float) or (type(atomInfo[7]) != float):
                raise (TypeError)
    except (IndexError):
        print('Uncomplete information for some atom!')
        raise
    except (TypeError):
        print('Atom property type error!')
        raise
    except (Exception):
        raise
    else:
        return (True)
    
        
# Define proDict is valid or not?
def proDictIsValid(proDict: dict) -> bool:
    # form and type test
    try:
        # form test
        if (len(set([len(proDict['Serial']), len(proDict['Atom']), \
                    len(proDict['ResName']), len(proDict['ResSeq']), \
                    len(proDict['ChainId']), len(proDict['X']), \
                    len(proDict['Y']), len(proDict['Z'])])) != 1):
            raise (IndexError)
        # type test
        if (any([type(i) != int for i in proDict['Serial']]) or \
            any([type(i) != str for i in proDict['Atom']]) or \
            any([type(i) != str for i in proDict['ResName']]) or \
            any([type(i) != int for i in proDict['ResSeq']]) or \
            any([type(i) != str for i in proDict['ChainId']]) or \
            any([type(i) != float for i in proDict['X']]) or \
            any([type(i) != float for i in proDict['Y']]) or \
            any([type(i) != float for i in proDict['Z'])):
            raise (TypeError)
    except (IndexError):
        print('Uncomplete information for some atom!')
        raise
    except (TypeError):
        print('Atom property type error!')
        raise
    else:
        return (True)


# This new type 'protein' is for operation pdb file.
class OldProtein:
    def __init__(self, pdbFile: str = '', infoDict: dict = '', infoList: dict = '', chainId: list[str] = ['all']):
        # PDB file as INPUT
        if (pdbFile):
            self.name: str = pdbFilePath.split(os.sep)[-1].rstrip('.pdb')
            self.infoList: list = pdbKit.pdb2List(pdbFilePath)
            _ = proListIsValid(self.infoList)
            self.infoDict: dict = pdbKit.pdb2List(pdbFilePath=pdbFile)
        # proDict as INPUT
        elif (infoDict):
            self.name: str = 'unnamed'
            _ = proDictIsValid(infoDict)
            self.infoDict: dict = infoDict
            self.infoList: list = pdbKit.proDict2ProList(infoDict)
        # proList as INPUT
        elif (infoList):
            self.name: str = 'unnamed'
            _ = proListIsValid(infoList)
            self.infoList: list = infoList
            self.infoDict: dict = pdbKit.proList2ProDict(infoList)
        # Remove unused atom
        if (chainId != ['all']):
            for atom in self.infoList:
                if (atom[4] not in chainId):
                    self.infoList.remove(atom)

    def __str__(self) -> str:
        return ('{name}: containing Chain {chainList}'.format(name = self.name, chainList = list(set(self.infoDict['ChainId']))))
    
    def __repr__(self) -> str:
        return ('{name}: containing Chain {chainList}'.format(name = self.name, chainList = list(set(self.infoDict['ChainId']))))

    def __sizeof__(self) -> int:
        return (len(set(self.infoDict['ChainId'])))

    def __len__(self) -> int:
        return (len(set(self.infoDict['ChainId'])))



