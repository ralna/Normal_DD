# Normal_DD
Preconditioner for Normal Equation Matrix
## __IMPORTANT:__
This code is a proof of concept to demonstrate the effectiveness of the preconditioner proposed in [1]. 
Note that this code may not be efficient as it is purely sequential. Performance can be achieved by using a parallel implementation of the method, see for example the implementation available in [PCHPDDM](https://petsc.org/main/docs/manualpages/PC/PCHPDDM/) in [PETSc](https://petsc.org/release/).
Potentially, the code can be used to test the preconditioner in a prototype code.

## Introduction
This is a proof of concept Matlab code for the two-level additive Schwarz preconditioner for the normal equaiton proposed in [1]

The code requires a Matlab interface for metis 5. Download [metis 5](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz) and [metismex.c](https://github.com/dgleich/metismex) then compile it on your platform.

## Running test examples
### Normal equation
The file example_normal.m is a script that compares the two-level additive Schwarz preconditioner against the one-level additive Schwarz preconditioner.

The parameters are the following:
- N: the number of subdomains
- k: Number of vectors contributed by each subdomain to the coarse space.
If kappa is provided and positive, k might be increased sufficiently to achieve 
a guaranteed upper bound on the condition number. (The best k is problem-dependent. Recommended value is 50, if kappa is negative and 50 is not enough then use double k and so on).
- kappa: desired bound on the condition number of the preconditioned matrix.
If kappa is negative, it is ignored. (Recommended value 100 but higher values may be desirable as they lead to faster setup and preconditioner application with certain increase in the iteration count)
- verbosity: control the information output
Note: if kappa and k are provided k will be ignored and the code will compute all necessary vectors to satisfy the bound on the condition number

### Weighted normal equation
The preconditioner can also be used to solve linear systems arising from the interior point method for solving Linear Programming problems.
The linear system matrix is usually represented in the augmented form

    A_aug = [T    B  ]
            [B^T  -sI]

where T is a diagonal matrix whose diagonal entries are strictly positive, B^T is the transpose of B (the constraint matrix), s is a regularization parameter, and I is the identity matrix.
Another representation is the reduced matrix

    A_red = B^T T^{-1} B + sI

The file example_lp.m is a script that compares the two-level additive Schwarz preconditioner against the one-level additive Schwarz preconditioner when solving the reduced system

The parameters are the following:
- N: the number of subdomains
- shift: the regularization parameter
- T: the weights (diagonal entries of the (1,1)-block in the augmented system)
- k: Number of vectors contributed by each subdomain to the coarse space.
If kappa is provided and positive, k might be increased sufficiently to achieve 
a guaranteed upper bound on the condition number. (The best k is problem-dependent. Recommended value is 50, if kappa is negative and 50 is not enough then use double k and so on).
- kappa: desired bound on the condition number of the preconditioned matrix.
If kappa is negative, it is ignored. (Recommended value 100 but higher values may be desirable as they lead to faster setup and preconditioner application with certain increase in the iteration count)
- verbosity: control the information output
Note: if kappa and k are provided k will be ignored and the code will compute all necessary vectors to satisfy the bound on the condition number

      [1] A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations. 
          SIAM Journal on Scientific Computing. 2022.
          Hussam Al Daas and Pierre Jolivet and Jennifer Scott
    
 To cite this work:
 
    @article{doi:10.1137/21M1434891,
    author = {Al Daas, Hussam and Jolivet, Pierre and Scott, Jennifer A.},
    title = {A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations},
    journal = {SIAM Journal on Scientific Computing},
    volume = {44},
    number = {3},
    pages = {A1047--A1068},
    year = {2022},
    doi = {10.1137/21M1434891},
    URL = {https://doi.org/10.1137/21M1434891},
    eprint = {https://doi.org/10.1137/21M1434891}
    }

## Acknowledgement
This work was supported by EPSRC grant EP/W009676/1.