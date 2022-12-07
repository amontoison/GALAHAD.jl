# Introduction

## Purpose

This package uses a primal-dual interior-point method
to **find a well-centered interior point** $x$ for a set of
general linear constraints
$\mbox{(1)} \;\; c_i^l\leqa_i^Tx\leq c_i^u, \;\;\; i = 1, \ldots , m,$
\n
 (1)c_i^l \[<=] a_i^Tx \[<=] c_i^u, i = 1, ... , m,
\n
and the simple bound constraints
$\mbox{(2)} \;\; x_j^l\leqx_j \leq x_j^u, \;\;\; j = 1, \ldots , n,$
\n
 (2) x_j^l \[<=] x_j \[<=] x_j^u, j = 1, ... , n,
\n
where the vectors
$a_{i}$, $c^l$, $c^u$, $x^l$ and $x^u$ are given.
More specifically, if possible, the package finds a solution to the
system ofprimal optimality equations
$\mbox{(3)} \;\; A x = c,$
\n
(3) A x = c,
\n
dual optimality equations
$\mbox{(4) $\hspace{3mm} g = A^T y + z, \;\; y = y^l + y^u, \;\mbox{and} \; z = z^l + z^u,$}$
\n
(4) g = A^T y + z, y = y^l + y^u and z = z^l + z^u,
\n
and perturbed complementary slackness equations
$\mbox{(5)} \;\;
( c_i - c^l_i ) y^l_i = (\mu_c^l)_i \;\mbox{and}\;
( c_i - c_i^u ) y^u_i = (\mu_c^u)_i, \;\;\;
 i = 1, \ldots , m,$
\n
(c_i - c^l_i) y^l_i = (mu_c^l)_i and
(c_i - c_i^u) y^u_i = (mu_c^u)_i, i = 1,...,m,
\n
and
$\mbox{(6)} \;\;
((x_j - x^l_j ) z_j^l = (\mu_x^l)_j\;\mbox{and}\;
( x_j - x^u_j ) z_j^u = (\mu_x^u)_j, \;\;\;
 j = 1, \ldots , n,$
\n
(x_j - c^l_j) z^l_j = (mu_x^l)_j and
(x_j - x_j^u) z^u_j = (mu_x^u)_i, j = 1,...,n,
\n
for which
\[
\mbox{(7)} \;\; c^l \leq c \leq c^u, \;\; x^l \leq x \leq x^u, \;\;
y^l \geq 0 , \;\;y^u \leq 0 , \;\; z^l \geq 0 \;\; \mbox{and} \;\; z^u \leq 0
\]
$$
\mbox{(7)} \;\; c^l \leq c \leq c^u, \;\; x^l \leq x \leq x^u, \;\;
y^l \geq 0 , \;\;y^u \leq 0 , \;\; z^l \geq 0 \;\; \mbox{and} \;\; z^u \leq 0
$$
\n
(7) c^l \[<=] c \[<=] c^u, x^l \[<=] x \[<=] x^u,
y^l \[>=] 0, y^u \[<=] 0, z^l \[>=] 0 and z^u \[<=] 0
\n
Here $A$ is the matrix whose rows are the $a_i^T$, $i = 1,
\ldots , m$, $\mu_c^l$, $\mu_c^u$, $\mu_x^l$ and
$\mu_x^u$ are vectors of strictly positive {\em targets}, $g$
is another given vector, and $(y^l, y^u)$ and $(z^l,
z^u)$ are dual variables for the linear constraints and simple
bounds respectively; $c$ gives the constraint value $A x$.
Since (5)-(7) normally imply that
\[
\mbox{(8)} \;\; c^l < c < c^u, \;\; x^l < x < x^u, \;\;
y^l > 0 , \;\;y^u < 0 , \;\; z^l > 0 \;\; \mbox{and} \;\; z^u < 0
\]
$$
\mbox{(8)} \;\; c^l < c < c^u, \;\; x^l < x < x^u, \;\;
y^l > 0 , \;\;y^u < 0 , \;\; z^l > 0 \;\; \mbox{and} \;\; z^u < 0
$$
\n
(8) c^l < c < c^u, x^l <; x < x^u,
y^l > 0, y^u < 0, z^l > 0 and z^u < 0
\n
such a primal-dual point $(x, c, y^l, y^u, z^l, z^l)$
may be used, for example, as a feasible starting point for primal-dual
interior-point methods for solving the linear programming problem
of minimizing $g^T x$ subject to (1) and (2).

Full advantage is taken of any zero coefficients in the vectors
$a_{i}$.Any of the constraint bounds $c_{i}^l$,
$c_{i}^u$, $x_{j}^l$ and $x_{j}^u$ may be infinite.
The package identifies infeasible problems, and problems for which
there is no strict interior, that is one or more of (8)
only holds as an equality for all feasible points.

## Authors

C. Cartis and N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

July 2006, C interface January 2022.

## Terminology

## Method

The algorithm is iterative, and at each major iteration attempts
to find a solution to the perturbed system (3), (4),
$\mbox{(9)}\;\;
( c_i - c^l_i + (\theta_c^l)_i )
( y^l_i + (\theta_y^l)_i )
= (\mu_c^l)_i \;\mbox{and}\;
( c_i - c_i^u - (\theta_c^u)_i )
( y^u_i - (\theta_y^u)_i )
= (\mu_c^u)_i, \;\;\;
 i = 1, \ldots , m,$
\n
 ( c_i - c^l_i + (theta_c^l)_i ) ( y^l_i + (theta_y^l)_i )
(9) = (mu_c^l)_i and
 ( c_i - c_i^u - (theta_c^u)_i )( y^u_i - (theta_y^u)_i )
= (mu_c^u)_i,i = 1,...,m
\n
$\mbox{(10)}\;\;
( x_j - x^l_j + (\theta_x^l)_j )
( z^l_j + (\theta_z^l)_j )
= (\mu_x^l)_j \;\mbox{and}\;
( x_j - x_j^u - (\theta_x^u)_j )
( z^u_j - (\theta_z^u)_j )
= (\mu_x^u)_j, \;\;\;
 j = 1, \ldots , n,$
\n
 ( x_j - x^l_j + (\theta_x^l)_j )( z^l_j + (\theta_z^l)_j )
(10) = (\mu_x^l)_j and
 ( x_j - x_j^u - (\theta_x^u)_j ) ( z^u_j - (\theta_z^u)_j )
 = (\mu_x^u)_j, j = 1,...,n,
\n
and
$\mbox{(11)}\;\;
c^l - \theta_c^l < c < c^u + \theta_c^u, \;\;
x^l - \theta_x^l < x < x^u + \theta_x^u, \;\;
 y^l > - \theta_y^l , \;\;
 y^u < \theta_y^u , \;\;
 z^l > - \theta_z^l \;\; \mbox{and} \;\;
 z^u < \theta_z^u ,$
\n
c^l - theta_c^l < c < c^u + theta_c^u,
x^l - theta_x^l < x < x^u + theta_x^u,
y^l > - theta_y^l, y^u < theta_y^u,
z^l > - theta_z^l and z^u < theta_z^,
\n
where the vectors of perturbations
$\theta^l_c$,$\theta^u_c$,$\theta^l_x$,$\theta^u_x$,
$\theta^l_x$,$\theta^u_x$,$\theta^l_y$,$\theta^u_y$,
$\theta^l_z$ and$\theta^u_z$,
are non-negative. Rather than solve (3)-(4) and (9)-(11) exactly,
we instead seek a feasible point for the easier relaxation(3)-(4) and
$\mbox{(12)}\;\;
\begin{array}{rcccll}
\gamma (\mu_c^l)_i & \leq &
( c_i - c^l_i + (\theta_c^l)_i ) ( y^l_i + (\theta_y^l)_i )
& \leq & (\mu_c^l)_i / \gamma & \mbox{and}\; \\
\gamma (\mu_c^u)_i & \leq &
( c_i - c_i^u - (\theta_c^u)_i ) ( y^u_i - (\theta_y^u)_i )
& \leq & (\mu_c^u)_i, /\gamma &
i = 1, \ldots , m, \;\mbox{and}\; \\
\gamma (\mu_x^l)_j & \leq &
( x_j - x^l_j + (\theta_x^l)_j ) ( z^l_j + (\theta_z^l)_j )
& \leq & (\mu_x^l)_j /\gamma & \mbox{and}\; \\
\gamma (\mu_x^u)_j& \leq &
( x_j - x_j^u - (\theta_x^u)_j )
( z^u_j - (\theta_z^u)_j )
& \leq & (\mu_x^u)_j /\gamma , &j = 1, \ldots , n,
\end{array}$
\n
 gamma (mu_c^l)_i
\[<=] ( c_i - c^l_i + (theta_c^l)_i ) ( y^l_i + (theta_y^l)_i )
\[<=](mu_c^l)_i / gamma and
 gamma (mu_c^u)_i
\[<=] ( c_i - c_i^u - (theta_c^u)_i ) ( y^u_i - (theta_y^u)_i )
 (12) \[<=](mu_c^u)_i, /gamma i = 1,...,m, and
 gamma (mu_x^l)_j
\[<=] ( x_j - x^l_j + (theta_x^l)_j ) ( z^l_j + (theta_z^l)_j )
\[<=](mu_x^l)_j /gamma and
 gamma (mu_x^u)_j
\[<=] ( x_j - x_j^u - (theta_x^u)_j ) ( z^u_j - (theta_z^u)_j )
\[<=](mu_x^u)_j /gamma , j = 1,...,n,
\n
for some $\gamma \in (0,1]$ which is allowed to be smaller than one
if there is a nonzero perturbation.

Given any solution to (3)-(4) and (12) satisfying (11),
the perturbations are reduced (sometimes to zero) so as to ensure that the
current solution is feasible for the next perturbed problem. Specifically,
the perturbation $(\theta^l_c)_i$ for the constraint $c_i \geq c^l_i$
is set to zero if $c_i$ is larger than some given parameter
$\epsilon > 0$.
If not, but $c_i$ is strictly positive, the perturbation will be
reduced by a multiplier $\rho \in (0,1)$. Otherwise, the new perturbation
will be set to $\xi (\theta^l_c)_i + ( 1 - \xi ) ( c_i^l - c_i )$
for some factor $\xi \in (0,1)$. Identical rules are used to reduce the
remaining primal and dual perturbations.
The targets $\mu_c^l$, $\mu_c^u$, $\mu_x^l$ and $\mu_x^u$
will also be increased by the factor $\beta \geq 1$ for those
(primal and/or dual) variables with strictly
positive perturbations so as to try to accelerate the convergence.

Ultimately the intention is to drive all the perturbations to zero.
It can be shown that if the original problem (3)-(6) and (8)
has a solution, the perturbations will be zero after a finite number of major
iterations. Equally, if there is no interior solution(8)
the sets of (primal and dual) variables that are necessarily at (one of) their
bounds for all feasible points---we refer to these as {\em implicit}
equalities---will be identified, as will the possibility that there is
no point (interior or otherwise) in the primal and/or dual feasible regions.

Each major iteration requires the solution $u = (x,c,z^l,z^u,y^l,y^u)$
of the nonlinear system (3), (4) and (9)-(11)
for fixed perturbations, using a minor iteration. The minor iteration
uses a stabilized (predictor-corrector) Newton method, in which the arc
$u(\alpha) = u + \alpha \dot{u} + \alpha^2 \ddot{u}, \alpha \in [0,1],$
$$u(\alpha) = u + \alpha \acute{u} + \alpha^2 \ddot{u}, \alpha \in [0,1],$$
u(alpha) = u + alpha u' + alpha^2 u”, alpha in [0,1], \
involving the standard Newton step
$\dot{u}$
&uacute;
u' \
for the equations (3), (4), (9) and (10), optionally augmented by a corrector
$\ddot{u}$
&uuml;
u” \
account for the nonlinearity in (9) and (10), is truncated so as to
ensure that
$(c_i(\alpha) - c^l_i + (\theta_c^l)_i)(y^l_i(\alpha) + (\theta_y^l)_i)
\geq \tau (\mu_c^l)_i \;\mbox{and}\;
(c_i(\alpha) - c_i^u - (\theta_c^u)_i)(y^u_i(\alpha) - (\theta_y^u)_i)
\geq \tau (\mu_c^u)_i, \;\;\; i = 1, \ldots , m,$
\n
(c_i(alpha) - c^l_i + (theta_c^l)_i)(y^l_i(alpha) + (theta_z^l)_i)
\[>=] tau (mu_c^l)_i and
(c_i(alpha) - c_i^u - (theta_c^u)_i ) (y^u_i(alpha) - (theta_z^u)_i)
\[>=] tau (mu_c^u)_i, i = 1,...,m
\n
and
$(x_j(\alpha) - x^l_j + (\theta_x^l)_j)(z^l_j(\alpha) + (\theta_z^l)_j)
\geq \tau (\mu_x^l)_j \;\mbox{and}\;
(x_j(\alpha) - x_j^u - (\theta_x^u)_j ) (z^u_j(\alpha) - (\theta_z^u)_j)
\geq \tau (\mu_x^u)_j, \;\;\; j = 1, \ldots , n,$
\n
(x_j(alpha) - x^l_j + (theta_x^l)_j)(z^l_j(alpha) + (theta_z^l)_j)
\[>=] tau (mu_x^l)_j and
(x_j(alpha) - x_j^u - (theta_x^u)_j ) (z^u_j(alpha) - (theta_z^u)_j)
\[>=] tau (mu_x^u)_j, j = 1,...,n
\n
for some $\tau \in (0,1)$, always holds, and also so that the norm
of the residuals to (3), (4), (9) and (10)
is reduced as much as possible.
The Newton and corrector systems are solved using a factorization of
the Jacobian of its defining functions (the so-called “augmented system”
approach) or of a reduced system in which some of the trivial equations are
eliminated (the “Schur-complement” approach).
The factors are obtained using the GALAHAD package SBLS.

In order to make the solution as efficient as possible, the
variables and constraints are reordered internally
by the GALAHAD package QPP prior to solution.
In particular, fixed variables, and
free (unbounded on both sides) constraints are temporarily removed.
In addition, an attempt to identify and remove linearly dependent
equality constraints may be made by factorizing
\[
\mat{cc}{\alpha I & A^T_E \\ A_E & 0},
\]
$$
 \left( \begin{array}{cc} \alpha I & A^T_E \\ A_E & 0 \end{array}, \right)
$$
\n
( alpha I A_E^T ),
(A_E0 )
\n
where $A_E$denotes the gradients of the equality constraints and
$\alpha > 0$ is a given scaling factor,
using the GALAHAD package SBLS, and examining small pivot blocks.

## Reference

The basic algorithm, its convergence analysis and results of
numerical experiments are given in

C. Cartis and N. I. M. Gould (2006).
Finding a point n the relative interior of a polyhedron.
Technical Report TR-2006-016, Rutherford Appleton Laboratory.

## Call order

To solve a given problem, functions from the wcp package must be called
in the following order:

- wcp\_initialize - provide default control parameters and
set up initial data structures
- wcp\_read\_specfile (optional) - override control values
by reading replacement values from a file
- wcp\_import - set up problem data structures and fixed
values
- wcp\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- wcp_find_wcp - find a well-centered point
- wcp\_information (optional) - recover information about
the solution and solution process
- wcp\_terminate - deallocate data structures

##  Unsymmetric matrix storage formats

The unsymmetric $m$ by $n$ constraint matrix $A$ may be presented
and stored in a variety of convenient input formats.

Both C-style (0 based)and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; the discussion below presumes C style, but add 1 to
indices for the corresponding fortran version.

Wrappers will automatically convert between 0-based (C) and 1-based
(fortran) array indexing, so may be used transparently from C. This
conversion involves both time and memory overheads that may be avoided
by supplying data that is already stored using 1-based indexing.

### Dense storage format
The matrix $A$ is stored as a compactdense matrix by rows, that is,
the values of the entries of each row in turn are
stored in order within an appropriate real one-dimensional array.
In this case, component $n \ast i + j$of the storage array A_val
will hold the value $A_{ij}$ for $0 \leq i \leq m-1$,
$0 \leq j \leq n-1$.

###  Sparse co-ordinate storage format

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $A$,
its row index i, column index j
and value $A_{ij}$,
$0 \leq i \leq m-1$,$0 \leq j \leq n-1$,are stored as
the $l$-th components of the integer arrays A_row and
A_col and real array A_val, respectively, while the number of nonzeros
is recorded as A_ne = $ne$.

###  Sparse row-wise storage format

Again only the nonzero entries are stored, but this time
they are ordered so that those in row i appear directly before those
in row i+1. For the i-th row of $A$ the i-th component of the
integer array A_ptr holds the position of the first entry in this row,
while A_ptr(m) holds the total number of entries plus one.
The column indices j, $0 \leq j \leq n-1$, and values
$A_{ij}$ of thenonzero entries in the i-th row are stored in components
l = A_ptr(i), $\ldots$, A_ptr(i+1)-1,$0 \leq i \leq m-1$,
of the integer array A_col, and real array A_val, respectively.
For sparse matrices, this scheme almost always requires less storage than
its predecessor.
