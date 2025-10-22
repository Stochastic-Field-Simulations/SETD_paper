# SETD paper

This is the code to generate data and plots of Model A, Model B and the KPZ equation in _Efficient Pseudospectral Algorithms for Statistical Field Theories_ by Martin Kjøllesdal Johnsrud and Navdeep Rana.
A preprint of the paper, where we derive the SETD method discuss implementing it for stochastic field simulations is available at [LINK].
If you use these methods or code from this repository, please cite:

`
Citation style
`

To clone repository, including the necessary `SFS` submodule, run

`git clone --recurse-submodules https://github.com/Stochastic-Field-Simulations/SETD_paper`

To activate Julia environment do run `project=SFS`, or select `SFS/` as current environment in VS Code.
All data generated to make the plots are available in `data/SETD_paper/`, or can be re-generated in the vaius `.jl` fiels.
The plots featured in the paper, plus a few more, are generated in the `.ipynb` python notebooks.
These files utilizes files from the `SFS/` module, including implementation of the numeric algorithms.
