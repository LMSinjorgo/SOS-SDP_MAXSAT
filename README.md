# On solving MAX-SAT using Sum of Squares
## By L. Sinjorgo and R. Sotirov (2024)

MATLAB implementation of the Peaceman-Rachford splitting method (PRSM) for computing bounds for the maximum-satisfiability problem (MAX-SAT).

Includes functions for parsing MAX-SAT instances, as given by CNF files, into the format required by the PRSM. More details and numerical results are presented in the paper, available [here](https://pubsonline.informs.org/doi/full/10.1287/ijoc.2023.0036).

Also includes the MATLAB implementation of the LOBPCG algorithm (``lobpcg.m``), as taken from <https://github.com/lobpcg/blopex>.


If you use this software, please cite it according to BibTeX file below.
````
@article{sinjorgo2024solving,
  title={On Solving MAX-SAT Using Sum of Squares},
  author={Sinjorgo, Lennart and Sotirov, Renata},
  journal={INFORMS Journal on Computing},
  volume={36},
  number={2},
  pages={417--433},
  year={2024},
  publisher={INFORMS}
}
````
