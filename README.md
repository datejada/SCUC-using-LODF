# Security Constrained Unit Commitment Using Line Outage Distribution Factors

- [Description](#description)
- [Citation](#citation)
- [Contact](#contact)
- [Files](#files)
- [How to Run the Model](#how-to-run-the-model)
- [Input Data](#input-data)
- [Output Results](#output-results)
- [License](#license)

## Description
Security-constrained unit commitment (SCUC) problem is one of the necessary tools for system operators to make operational planning and real-time operation. The internalization of transmission network and security constraints (e.g., N-1 criterion) could lead to different decisions in the generation dispatch. However, the computational burden of this problem is challenging mainly due to its inherent large problem size. This model proposes an N-1 security constrained formulation based on the line outage distribution factors (LODF) instead of the one based on injection sensitivity factors (ISF). Additionally, an iterative methodology for filtering the active N-1 congestion constraints is used. The computational efficiency of the proposed model is shown by solving the SCUC of the IEEE 118 bus system.

The model is documented in:
[D. A. Tejada-Arango, P. Sánchez-Martın and A. Ramos, "Security Constrained Unit Commitment Using Line Outage Distribution Factors," in IEEE Transactions on Power Systems, vol. 33, no. 1, pp. 329-337, Jan. 2018.
doi: 10.1109/TPWRS.2017.2686701](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7886335&isnumber=8231802)

The entire model is coded in [GAMS](https://www.gams.com/) and solved with [GUROBI](http://www.gurobi.com/).

## Citation
Whenever you use this code, please refer to:

[D. A. Tejada-Arango, P. Sánchez-Martın and A. Ramos, "Security Constrained Unit Commitment Using Line Outage Distribution Factors," in IEEE Transactions on Power Systems, vol. 33, no. 1, pp. 329-337, Jan. 2018.
doi: 10.1109/TPWRS.2017.2686701](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7886335&isnumber=8231802)

## Contact
Diego Alejandro Tejada-Arango, Universidad Pontificia Comillas, dtejada@comillas.edu

## Files
* StarGenLite_SCUC.xlsm **Microsoft Excel interface** with IEEE 118 bus system information
* StarGenLite_SCUC.gms  has the GAMS code
* StarGenLite_SCUC.gpr  has the GAMS project

## How to Run the Model
There are two possible ways to run the model.

1. _From the Microsoft Excel interface_: Click on the Run button. A command window will appear showing the problem being solved. After the execution click on the Load results button and the output results will be loaded in the corresponding worksheets of the Excel interface.
2. _GAMSIDE Application_: Open GAMS, create a new project in the directory where you have put the SCUC files and then open **StarGenLite_SCUC.gms**. The user parameter user1 has to be defined in the options window with the name of the Excel interface, e.g., **user1="StarGenLite_SCUC.xlsm"**. In both modes a Microsoft Excel workbook named **tmp.xlsx** that is read the Excel interface clicking the Load results button. Besides the Microsoft Excel interface, an additional file is written after the execution. It is the GAMS listing, **StarGenLite_SCUC.lst**. This file may display execution errors and always have to be checked to look for them.

## Input Data
Input data are located in named ranges contained in green worksheets of the Microsoft Excel interface. The user can introduce or delete input data but the named ranges have to be kept. To avoid errors, it is advisable to insert/delete rows or columns between the first and the last ones, i.e., in any intermediate row or column. In order to ease the introduction of data the interface follows this coloring code:
* blue is reserved for comments,
* black for headings of rows or columns and
* green for data introduced by the user.
The light green cells can be modified by the user. It is suggested that the other ones are not modified.
Some cells have conditional format to help in detecting inconsistencies in input data.
The writing criteria follow the GAMS syntaxes. For example, the * in the first column indicates a comment.

## Output Results
The output results are in salmon worksheets of the Microsoft Excel interface.

## License
[License file](../blob/master/LICENSE)
