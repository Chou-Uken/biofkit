---
sort: 2
---

# Module - oldProtein

这个`module`仅在alpha开发时期维护，在初期作为权宜之计创建。

## 导入
```python3
from biofkit.proteinKit import oldProtein
```

## PDB文件读取
### Function - pdb2List
加载PDB文件并返回一个列表，输出列表由若干个表示单个原子信息的列表构成，子列表的元素顺序为`[原子序号(int), 原子类型(str), 残基名(str), 残基标号(int), 链名(str), X坐标(float), Y坐标(float), Z坐标(float)]`。
```python3
def pdb2List(pdbFilePath: str, csvPath: str = '', colName: bool = False) -> list[list[int, str, str, int, str, float, float, float]]:
  """读取PDB文件并转换为列表

  参数：
    pdbFilePath (str)：PDB文件的路径
    csvPath (str)：如果需要输出一个csv文件，则填入欲写入文件的路径，不存在则创建
    colName (bool)：列表（和）csv文件是否需要列名

  输出：
    list：蛋白质结构（列表形式）
  """

例：
proList = pdb2List(pdbFilePath='~/example.pdb')
```

### Function - pdb2Dict
加载PDB文件并输出一个字典，键为原子的属性，值为所有原子该属性组成的列表。`{'Serial: [], 'Atom': [], 'ResName': [], 'ResSeq': [], 'ChainId': [], 'X': [], 'Y': [], 'Z': []}`
```python3
def pdb2Dict(pdbFilePath: str) -> dict[str, list]:
  """读取PDB文件并输出一个字典
  
  参数：
    pdbFilePath (str)：PDB文件的路径

  输出：
    dict：蛋白质结构（字典形式）
  """

例：
proDict = pdb2Dict(pdbFilePath='~/example.pdb')
```

### Function - proDict2ProList
把一个包含蛋白结构的字典转换为列表，函数包含了查错机制。
```python3
def proDict2ProList(rawDict: dict) -> list:
  """转换蛋白列表为蛋白字典

  参数：
    rawDict (dict)：蛋白列表

  输出：
    dict：蛋白字典
  """
```

### Function - proList2ProDict
{% include list.liquid %}
