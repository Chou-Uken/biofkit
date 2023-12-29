# biofkit
Biofkit is now containing only one module with one package which can deal with pdb file (extract sequence, atom information and so on) very easily. 

## How to Install
### Conda
```console
conda install -c chou_uken biofkit
```
### Mamba
``` 
mamba install -c chou_uken biofkit
```

### Pip
```console
pip install biofkit
```

## How to Use
### ProteinKit
```python
from proteinKit import pdbKit

# if argument fasta is True, a fasta file will be created in the same path as the pdb file.
pdbKit.pdb2Seq(pdbFilePath: str, fasta: bool = False, fastaLineLen: int = 80) -> dict[str, str]
    # pdfFilePath: the file path of the pdb you want to transfer.
    # fasta: Whether to transfer into fasta file. If false, pdb will only be transferred to a dictionary.{chainId: Seq}
    # fastaLineLen: How many residues are contained in a single line of the fasta, only work when `fasta` is true.

# load the information of all amino-acid-residue atoms into a list. output[idx] shows the information of an atom.
pdb2List(pdbFilePath: str, csvPath: str = None, colName: bool = False) -> list[list[int, str, str, int, str, float, float, float]]:
    # pdbFilePath: the file path of the pdb you want to load.
    # csvPath: if given, then write a csv file for your atom information.
    # colName: a prepared colname list for the output. If colName is true, the output[0] will be such a list shown below. May help when column names are needed.
    pdbInfoColumns: [str] = ['Serial', 'Atom', 'ResName', 'ResSeq', 'ChainId', 'X', 'Y', 'Z']

# load the information of all amimo-acid-residue atoms into a dictionary, which can be converted into a dataframe with famous `pandas`.
pdb2Dict(pdbFilePath: str) -> dict[str, list]:
    # pdbFilePath: the file path of the pdb you want to load.

```
