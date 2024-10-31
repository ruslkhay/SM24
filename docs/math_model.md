# Mathematical model of solution

Dirichlet problem in $D$ area with boarder $\gamma$ :
```math
k(x, y) = 
\begin{cases} 
-\Delta{u} = f(x,y), & (x,y) \in D \\ 
u(x, y) = 0, & (x,y) \in \gamma
\end{cases} 
\tag{1}
```

Let $\Pi$ be a rectangle with vertices 
$A_1(0,0),\space B_1(3,0),\space A_2(0,3),\space B_2(0,3)$.
Let the area $D$ belong to the rectangle $\Pi$. We denote the closure of area
$D$ and the rectangle $\Pi$ by $D$ and $\Pi$ respectively, and the boundary of
the rectangle by $\Gamma$. The difference of the sets 
$\hat{D} = \Pi \setminus D$ is called the fictitious area. 
Let's choose and fix a small $\varepsilon > 0$. In the rectangle $\Pi$ we 
consider the Dirichlet problem

```math
-\frac{\partial}{\partial x} \left(
    k(x, y) \frac{\partial v}{\partial x} 
    \right) - 
    \frac{\partial}{\partial y} \left( 
        k(x, y) \frac{\partial v}{\partial y} 
        \right) = F(x, y), 
\tag{2}
```
with the boundary condition 

```math
v(x, y) = 0, \quad (x, y) \in \Gamma.
```

with a piecewise constant coefficient
```math
k(x, y) = 
\begin{cases} 
1, & (x, y) \in D, \\ 
\frac{1}{\varepsilon}, & (x, y) \in \hat{D} 
\end{cases} 
\tag{3}
```
and the right-hand side

```math
F(x, y) = 
\begin{cases} 
f(x, y), & (x, y) \in D, \\ 
0, & (x, y) \in \hat{D} 
\end{cases} 
\tag{4}
```

It is required to find a continuous function $v(x, y)$ in $\Pi$ that 
satisfies the differential equation of problem (3) everywhere in $\Pi$ 
$\setminus \gamma$, is equal to zero on the boundary $\Gamma$ of the rectangle, 
and such that the flow vector

```math
W(x, y) = -k(x, y) \left(
    \frac{\partial v}{\partial x}, \frac{\partial v}{\partial y}
    \right)
```

has a continuous normal component on the common part of the curvilinear boundary
of the domain $D$ and the rectangle $\Pi$. The latter means that in each point
$(x_0, y_0) \in \gamma \cap \Pi$ the following equality must hold:

```math
\lim_{D \ni (x,y) \to (x_0,y_0)} (W(x, y), n(x_0, y_0)) = \lim_{\hat{D} \ni (x,y) \to (x_0,y_0)}(W(x, y), n(x_0, y_0)),
\tag{5}
```

where $n(x, y)$ is the unit normal vector to the boundary $\gamma$ at the point
$(x, y)$, defined everywhere or almost everywhere on the curve.

It is known that the function $v(x, y)$ uniformly approximates the solution
$u(x, y)$ of problems (1) and (2) in the domain $D$, namely,

```math
\max_{P \in D} |v(x, y) - u(x, y)| < C_\epsilon, \; C_\epsilon > 0.
```

