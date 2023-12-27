# biofkit
Biofkit is now containing only one module with one package which can transfer pdb into fasta (or a sequence dictionary) very easily. 

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

# load the information of all amino-acid-residue atoms into a list which can be converted to a dataframe with famous `pandas`.
pdbKit.pdb2dfList(pdbFilePath: str, colName: bool = True) -> list[list[int, str, str, int, str, float, float, float]]:
    # pdbFilePath: the file path of the pdb you want to load.
    # colName: a prepared colname list for the output. If colName is true, pdb can be directly transfer to a pandas dataframe.
```