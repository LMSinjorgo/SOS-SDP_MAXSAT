# On solving the MAX-SAT using sum of squares
## By L. Sinjorgo and R. Sotirov (February, 2023)

MATLAB implementation of the Peaceman-Rachford splitting method (PRSM) for computing bounds for the maximum-satisfiability problem (MAX-SAT).

Includes functions for parsing MAX-SAT instances, as given by CNF files, into the format required by the PRSM. More details and numerical results are presented in the paper, available [here](https://arxiv.org/abs/2302.06931).

Also includes the MATLAB implementation of the LOBPCG algorithm (``lobpcg.m``), as taken from <https://github.com/lobpcg/blopex>.


If you use this software, please cite it according to BibTeX file below.
>@article{Sinjorgo2023_SOSandMAXSAT,  
>  title={On solving the MAX-SAT using sum of squares},  
>  author={Sinjorgo, Lennart and Sotirov, Renata},  
>  journal={arXiv preprint arXiv:2302.06931},  
>  year={2023}  
>}
