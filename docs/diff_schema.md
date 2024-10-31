# Difference Scheme for Solving the Problem

The problem (3) is proposed to be solved numerically using the finite difference
method. In the notation of adjoining $\Pi$, we define a uniform rectangular grid
$\bar{\omega_h} = \omega_1 \times \omega_2$, where

```math
\bar{\omega_1} = \{x_i = A_1 + ih_1, \quad i = \overline{0, M}\}, \quad 
\bar{\omega_2} = \{y_j = A_2 + jh_2, \quad j = \overline{0, N}\}.
```

Here, $h_1 = (B_1 - A_1)/M$, $h_2 = (B_2 - A_2)/N$. The set of interior nodes of
the grid \bar{$\omega_h$} is denoted by $\omega_h$, which includes the sets of
grid points of the rectangle that do not lie on the boundary $Г$.

Let’s consider the linear space $H$ of functions defined on the grid $w$.
We denote by $w_{ij}$ the value of the grid function $w$ at the node of the grid
$(x_i, y_j) \in \omega_h$. We will assume that the inner product and Euclidean
norm in the space $H$ are defined as follows:

```math
(u, v) = \sum_{i=1}^{M-1} \sum_{j=1}^{N-1} h_1h_2u_{ij}v_{ij},
\quad \|u\|_E = \sqrt{(u, u)}.
```

In the finite difference method, the differential problem of mathematical
physics is replaced by the finite difference operation of a problem of the form

```math
Aw = B,
```

where $A : H \to H$ is an operator acting in the space of grid functions, and
$B \subset H$ is the known right-hand side. The task (9) is called a difference
scheme. The solution to this problem is considered a numerical solution to the 
original differential problem.

When constructing a difference scheme, one should approximate
(replace approximately) the equations of the boundary problem with their
difference analogs—grid equations that relate the values of the sought function.


The differential equation of the problem (3) is approximated by a difference
equation at all interior points of the grid:

```math
-\frac{1}{h_1} \left(
    a_{i+1,j} \frac{w_{i+1,j} - w_{ij}}{h_1} - 
    a_{ij} \frac{w_{i,j} - w_{i-1,j}}{h_1} 
    \right) -  \\
-\frac{1}{h_2} \left(
    b_{i,j+1} \frac{w_{i,j+1} - w_{ij}}{h_2} -
    b_{ij} \frac{w_{ij} - w_{i,j-1}}{h_2}
    \right) = \\
= F_{ij},
```

where 

```math
a_{ij} = \frac{1}{h_2} \int_{y_{j-\frac{1}{2}}}^{y_{j+\frac{1}{2}}} 
k(x_{i-\frac{1}{2}}, t) dt,
```

```math
b_{ij} = \frac{1}{h_1} \int_{x_{i-\frac{1}{2}}}^{x_{i+\frac{1}{2}}} 
k(t, y_{j-\frac{1}{2}}) dt,
```

for all $i = \overline{1, M}$, $j = \overline{1, N}$. Here, the half-nodes are
defined as

```math
x_{i\pm\frac{1}{2}} = x_i \pm 0.5h_1, \quad 
y_{j\pm\frac{1}{2}} = y_j \pm 0.5h_2.
```