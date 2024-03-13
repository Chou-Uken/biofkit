# This script was created on Feb 22nd, 2024 by Zhang Yujian as a doctoral candidate in Institute of Zoology, CAS.
# Thanks for using. Please report bugs (if any) at zhangyujian23@mails.ucas.ac.cn.

class Atom(object):
    """A class containing number, type and coordinate of an atom.

    Attributes:
        Serial (int): Id of atom.
        Atom (str): Type of atom.
        X (float): Coordinate of x.
        Y (float): Coordinate of y.
        Z (float): Coordinate of z.
    """

    def __init__(self, serial: int, atom: str, x: float, y: float, z: float):
        """Inits Atom with its id, type and coordinate.
        
        Args:
            Serial (int): Id of atom.
            Atom (str): Type of atom.
            X (float): Coordinate of x.
            Y (float): Coordinate of y.
            Z (float): Coordinate of z.
        """

        self.Serial: int = int(serial)
        self.Atom: str = str(atom)
        self.X: float = float(x)
        self.Y: float = float(y)
        self.Z: float = float(z)
    
    def __str__(self) -> str:
        return ('Atom: {name} at location of ({x}, {y}, {z})'.format(name=self.Atom, x=self.X, y=self.Y, z=self.Z))

    def __repr__(self) -> str:
        return ('Atom: {name} at location of ({x}, {y}, {z}))'.format(name=self.Atom, x=self.X, y=self.Y, z=self.Z))

    def getCoord(self) -> dict[str, float]:
        return ({'x': self.X, 'y': self.Y, 'z': self.Z})


class Residue(object):
    """A class containing number, type and its Atoms of a residue in a peptide.

    Attributes:
        ResSeq (int): Id of residue.
        ResName (str): Type of residue.
        AtomSet (list[Atom]): Set of Atom belonging to this residue.
    """

    def __init__(self, atomList: list[Atom], resSeq: int = 0, resName = '*'):
        """Inits residue with its id, type and atoms.

        Args:
            ResSeq (int): Id of residue.
            ResName (str): Type of residue.
            AtomSet (list[Atom]): Set of Atoms belonging to this residue.
        """

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

    def __getitem__(self, idx: int) -> Atom:
        return (self.AtomSet[idx])

    def getName(self) -> str:
        """To get type of this residue.
        
        Returns:
            str: Type of this residue.
        """

        return (self.ResName)
        
        
class Peptide(object):
    """a class of peptide with its id and set of residues belonging to it.

    Attributes:
        ChainId (str): Id of this chain.
        ResSet (list[Residue]): Set of Residues belonging to this chain.
    """

    def __init__(self, resList: list[Residue], chainId: str = 'A'):
        """Inits Peptide with its id and residues belonging to it.
        
        Args:
            ChainId (str): Id of this chain.
            ResSet (list[Residue]): Set of Residues belonging to this chain.
        """

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

    def __getitem__(self, idx: int) -> Residue:
        return (self.ResSet[idx])

    def getChainId(self) -> str:
        """To get id of this chain.
        
        Returns:
            str: Id of this chain.
        """

        return (self.ChainId)
    
    def cut(self, end: int, beginning: int = 0, newChainId: str = 'X'):
        """To cut a peptide at specific site.
        
        Args:
            end (int): Id of the residue at the end of output.
            beginning (int): Id of the residue at the beginning of output.
            newChainId (str): Id of the Peptide you created.

        Returns:
            Peptide: A new peptide.
        """

        if ((beginning == 0) and (end == len(self.ResSet))):
            print('An identical peptide is created!')
        return (Peptide(resList=self.ResSet[beginning: end], chainId=newChainId))
        


class Protein(object):
    """a class containing protein with its peptides.

    Attributes:
        ame (str): Name of this protein.
        pepSet (list[Peptide]): The list of Peptides belonging to the Protein.
    """

    def __init__(self, pepList: list[Peptide], proteinName: str = 'Unnamed'):
        """Inits Protein with its Peptides and its name.

        Args:
            pepList (list[Peptide]): a list comprising of Peptides.
            name (str): Name of this protein.
        """

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
            self.pepSet: list[Peptide] = pepList
            
    def __str__(self) -> str:
        return (self.name + ': ' + str(self.pepSet))

    def __repr__(self) -> str:
        return (self.name + ': ' + str(self.pepSet))

    def __getitem__(self, idx: int) -> Peptide:
        return (self.pepSet[idx])

    def pick(self, chainId: str) -> list[Peptide]:
        """pick specific peptide from a Protein.
        
        Args:
            chainId (str): Id of selected peptide.

        Returns:
            list[Peptide]: a list of Peptides with given Id.
        """

        output: list[Peptide] = []
        for peptide in self.pepSet:
            if (peptide.getChainId() == chainId):
                output.append(peptide)
        return (output)

