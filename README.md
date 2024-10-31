# Installation

Clone repository and run from root directory:
```
mkdir build && cd build 
cmake ..
make
```

## For developers

Additionally you can download python package manager and create venv
```
curl -sSL https://pdm-project.org/install-pdm.py | python3 -
pdm install
```

# Overview:

This program solves Dirichlet problem
```math
k(x, y) = 
\begin{cases} 
-\Delta{u} = f(x,y), & (x,y) \in D \\ 
u(x, y) = 0, & (x,y) \in \gamma
\end{cases} 
```


for a D, that is right trapezoid with next nodes:
$A(0,0)$, $B(3,0)$, $C(2,3)$, $D(0,3)$, - $\gamma$ is its boarder and $f(x, y) = 1 \quad \forall (x,y) \in D$.

The goal is to build accurate and parallel solution, using **OpenMP** and 
**MPI** libraries. Moreover launch it on Moscow State University 
[Polus](http://hpc.cmc.msu.ru/polus) computing complex.


## Project description

1. [Mathematical model of solution](docs/math_model.md).
2. [Difference Scheme for Solving the Problem](docs/diff_schema.md)
3. [Method for Solving Systems of Linear Algebraic Equations](docs/math_model.md)
4. [OpenMP implementation](docs/omp.md)
5. [Solution plots](docs/sol_images.md)