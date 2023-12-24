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

## How to Use
### ProteinKit
```python
from proteinKit import pdbKit

# if argument fasta is True, a fasta file will be created in the same path as the pdb file.
pdbKit.pdb2Seq(pdbFilePath: str, fasta: bool = False) -> dict[str, str]
```