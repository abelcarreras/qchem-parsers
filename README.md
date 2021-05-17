QChem-parsers
=============
This is a repository of python parsers for Q-Chem (https://www.q-chem.com)

How to use
----------
These parsers are simple python functions with a single
positional argument. This argument  corresponds to a
Q-Chem output  in plain text (string) format. 
The return of the parser is a python dictionary in which 
entry corresponds to a parsed property.

```python
from qcparsers.parsers import parser_optimization

with open('output_file.out') as f:
    qc_output = f.read()

molecule = parser_optimization(qc_output)
```

Version system
--------------
As optional feature the parsers can include a docstring with
the list of compatible Q-Chem versions. 

```python
def parser_basic(output):
    """
    This showcases the format of  Q-Chem version parser compatibility.
    Just create a docstring with the following line:

    compatibility: 5.1, 5.2+

    more text can be added to the docstring 

    """

```
This information can be extracted using provided helper functions

```python
from qcparsers.tools.version import get_version_output
from qcparsers.tools.version import get_compatibility_list_from_parser
from qcparsers.parsers import parser_basic


with open('output_file.out') as f:
    output = f.read()

version = get_version_output(output)
compatibility_list = get_compatibility_list_from_parser(parser_basic)

if version in compatibility_list:
    print('Parser compatible')
else:
    print('Parser not compatible')
```

How to contribute
-----------------
Collaboration is welcomed. 
Contribution is specially needed in this points:
- Expanding: Adding new parsers and tools into repository
- Fixing: remove bugs in existing implementation
- Testing: Adding tests and finding limitations in existing parsers

To contribute to this repository use the following steps:
- Fork this repository using GitHub web interface
- Clone your fork to your local machine.
- Create a new branch using git setting a name related to your contribution
- Checkout this branch and make your changes/contributions to the code
- Push new commits to Github
- On Github create a pull request merging your new branch into master of the original project
