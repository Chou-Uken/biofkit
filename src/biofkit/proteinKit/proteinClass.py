# This script was created on Feb 22nd, 2024 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
# Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.

class Atom(object):
    def __init__(self, serial: int, atom: str, x: float, y: float, z: float):
        self.Serial: int = int(serial)
        self.Atom: str = str(atom)
        self.X: float = float(x)
        self.Y: float = float(y)
        self.Z: float = float(z)
    
    def __str__(self) -> str:
        return ('Atom: {name} at location of ({x}, {y}, {z})'.format(name=self.Atom, x=self.X, y=self.Y, z=self.Z))
    
    def __repr__(self) -> str:
        return ('Atom: {name} at location of ({x}, {y}, {z}))'.format(name=self.Atom, x=self.X, y=self.Y, z=self.Z))


class Residue(object):
    def __init__(self, atomList: list[Atom], resSeq: int = 0, resName = '*'):
        self.ResSeq: int = int(resSeq)
        self.ResName: str = str(resName)
        # type test
        try:
            for atom in atomList:
                if (type(atom) != Atom):
                    raise (TypeError)
        except (TypeError):
            print('TypeError! atomList must only contain elements of Atom type.')
            raise
        except (Exception):
            raise
        else:
            self.AtomSet: list[Atom] = atomList

    def __str__(self) -> str:
        return (self.ResName + str(self.ResSeq))

    def __repr__(self) -> str:
        return (self.ResName + str(self.ResSeq))

    def getName(self) -> str:
        return (self.ResName)
        
        
class Peptide(object):
    def __init__(self, resList: list[Residue], chainId: str = 'A'):
        self.ChainId: str = chainId
        # type test
        try:
            for res in resList:
                if (type(res) != Residue):
                    raise (TypeError)
        except (TypeError):
            print('TypeError! resList must only contain elements of Residue type.')
            raise
        except (Exception):
            raise
        else:
            self.ResSet: list[Residue] = resList
    
    def __str__(self) -> str:
        return ('Chain {name}: {seq}'.format(name=self.ChainId, seq=''.join([i.getName() for i in self.ResSet])))

    def __repr__(self) -> str:
        return ('Chain {name}: {seq}'.format(name=self.ChainId, seq=''.join([i.getName() for i in self.ResSet])))

    def __len__(self) -> int:
        return (len(self.ResSet))

    def getChainId(self) -> str:
        return (self.ChainId)
    
    def cut(self, end: int, beginning: int = 0, newChainId: str = 'X'):
        if ((beginning == 0) and (end == len(self.ResSet))):
            print('An identical peptide is created!')
        return (Peptide(resList=self.ResSet[beginning: end], chainId=newChainId))
        


class Protein(object):
    def __init__(self, pepList: list[Peptide], proteinName: str = 'Unnamed'):
        self.name = proteinName
        # type test
        try:
            for pep in pepList:
                if (type(pep) != Peptide):
                    raise (TypeError)
        except (TypeError):
            print('TypeError! PepList must only contain elements of Peptide type.')
            raise
        except (Exception):
            raise
        else:
            self.PepSet: list[Peptide] = pepList
            
    def __str__(self) -> str:
        return (self.name + ': ' + str(self.PepSet))

    def __repr__(self) -> str:
        return (self.name + ': ' + str(self.PepSet))

    def pick(self, chainId: str) -> list[Peptide]:
        output: list[Peptide] = []
        for peptide in self.PepSet:
            if (peptide.getChainId() == chainId):
                output.append(peptide)
        return (output)

