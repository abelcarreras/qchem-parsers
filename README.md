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
