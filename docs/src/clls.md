# Introduction

## Purpose

This package uses a primal-dual interior-point crossover method
to solve the **constrained linear least-squares problem**
$\mbox{minimize}\;\; f(x) = \frac{1}{2} \| A x - b \|^2
+ \frac{1}{2} \sigma \| x \|^2 $
\n
minimize f(x) := 1/2 ||Ax-b||^2 + 1/2 sigma ||x||^2
\n
subject to the general linear constraints
$c_i^l\leql_i^Tx\leq c_i^u, \;\;\; i = 1, \ldots , m,$
\n
 c_i^l \[<=] l_i^Tx \[<=] c_i^u, i = 1, ... , m,
\n
and the simple bound constraints
$x_j^l\leqx_j \leq x_j^u, \;\;\; j = 1, \ldots , n,$
\n
 x_j^l \[<=] x_j \[<=] x_j^u, j = 1, ... , n,
\n
where the $m$ by $n$ matrix
$A$, the vectors $l_i$, $c^l$, $c^u$, $x^l$,
$x^u$ and the scalar $\sigma$ are given.
Any of the constraint bounds $c_i^l$, $c_i^u$,
$x_j^l$ and $x_j^u$ may be infinite.
Full advantage is taken of any zero coefficients in the matrix $A$
or the matrix $A$ of vectors $l_i$.

## Authors

N. I. M. Gould and D. P. Robinson, STFC-Rutherford Appleton Laboratory,
England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

July 2022.

## Terminology

The required solution $x$ necessarily satisfies
the primal optimality conditions
$\mbox{(1a) $\hspace{66mm} A x = c\hspace{66mm}$}$
\n
(1a) L x = c
\n
and
$\mbox{(1b) $\hspace{52mm} c^l \leq c \leq c^u, \;\; x^l \leq x \leq x^u,\hspace{52mm}$} $
\n
(1b) c^l \[<=] c \[<=] c^u, x^l \[<=] x \[<=] x^u,
\n
the dual optimality conditions
$\mbox{(2a) $\hspace{3mm} A^T ( Ax-b ) = L^T y + z$
\n
(2a) A^T ( A x - b ) = A^T y + z
\n
where
$\mbox{(2b) $\hspace{24mm} y = y^l + y^u, \;\; z = z^l + z^u, \,\,
 y^l \geq 0 , \;\;y^u \leq 0 , \;\;
 z^l \geq 0 \;\; \mbox{and} \;\; z^u \leq 0,\hspace{24mm}$} $
\n
 (2b) y = y^l + y^u, z = z^l + z^u, y^l \[>=] 0, y^u \[<=] 0,
z^l \[>=] 0 and z^u \[<=] 0,
\n
and the complementary slackness conditions
$\mbox{(3) $\hspace{12mm}
( A x - c^l )^T y^l = 0,\;\;( A x - c^u )^T y^u = 0,\;\;
(x -x^l )^T z^l = 0 \;\;\mbox{and} \;\; (x -x^u )^T z^u = 0,\hspace{12mm} $}$
\n
(3) (A x - c^l)^T y^l = 0, (A x - c^u)^T y^u = 0,
(x -x^l)^T z^l = 0 and (x -x^u)^T z^u = 0,
\n
where the diagonal matrix $W^2$ has diagonal entries $w_j^2$,
$j = 1, \ldots , n$, where the vectors $y$ and $z$ are
known as the Lagrange multipliers for
the general linear constraints, and the dual variables for the bounds,
respectively, and where the vector inequalities hold component-wise.

## Method

Primal-dual interior point methods iterate towards a point
that satisfies these conditions by ultimately aiming to satisfy
(1a), (2a) and (3), while ensuring that (1b) and (2b) are
satisfied as strict inequalities at each stage.Appropriate norms of the
amounts bywhich (1a), (2a) and (3) fail to be satisfied are known as the
primal and dual infeasibility, and the violation of complementary slackness,
respectively. The fact that (1b) and (2b) are satisfied as strict
inequalities gives such methods their other title, namely
interior-point methods.

The method aims at each stage to reduce the
overall violation of (1a), (2a) and (3),
rather than reducing each of the terms individually. Given an estimate
$v = (x, c, y, y^l, y^u, z, z^l, z^u)$
of the primal-dual variables, a correction
$\Delta v = \Delta (x, c, y, y^l, y^u z, z^l, z^u)$
is obtained by solving a suitable linear system of Newton equations for the
nonlinear systems (1a), (2a) and a parameterized “residual
trajectory” perturbation of (3); residual trajectories
proposed by Zhang (1994) and Zhao and Sun (1999) are possibilities.
An improved estimate $v + \alpha \Delta v$
is then used, where the step-size $\alpha$
is chosen as close to 1.0 as possible while ensuring both that
(1b) and (2b) continue to hold and that the individual components
which make up the complementary slackness
(3) do not deviate too significantly
from their average value. The parameter that controls the perturbation
of (3) is ultimately driven to zero.

The Newton equations are solved by applying the
GALAHAD matrix factorization package SBLS, but there are options
to factorize the matrix as a whole (the so-called "augmented system"
approach), to perform a block elimination first (the "Schur-complement"
approach), or to let the method itself decide which of the two
previous options is more appropriate.
The "Schur-complement" approach is usually to be preferred when all the
weights are nonzero or when every variable is bounded (at least one side),
but may be inefficient if any of the columns of $A$ is too dense.

Optionally, the problem may be pre-processed temporarily to eliminate
dependent constraints using the GALAHAD package FDC. This may
improve the performance of the subsequent iteration.

## Reference

The basic algorithm is a generalisation of those of

Y. Zhang (1994),
 On the convergence of a class of infeasible interior-point methods for the
 horizontal linear complementarity problem,
 SIAM J. Optimization 4(1) 208-227,

and

G. Zhao and J. Sun (1999).
On the rate of local convergence of high-order infeasible path-following
algorithms for the $P_\ast$ linear complementarity problems,
Computational Optimization and Applications 14(1) 293-307,

with many enhancements described by

N. I. M. Gould, D. Orban and D. P. Robinson (2013).
Trajectory-following methods for large-scaledegenerate convex quadratic
programming,
Mathematical Programming Computation 5(2) 113-142.

## Call order

To solve a given problem, functions from the clls package must be called
in the following order:

- clls\_initialize - provide default control parameters and
set up initial data structures
- clls\_read\_specfile (optional) - override control values
by reading replacement values from a file
- clls\_import - set up problem data structures and fixed
values
- clls\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- solve the problem by calling one of
 - clls\_solve_qp - solve the quadratic program
 - clls\_solve_sldqp - solve the shifted least-distance problem
- clls\_information (optional) - recover information about
the solution and solution process
- clls\_terminate - deallocate data structures

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

##  Symmetric matrix storage formats

Likewise, the symmetric $n$ by $n$ objective Hessian matrix
$H$ may be presented
and stored in a variety of formats. But crucially symmetry is exploited
by only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal).

### Dense storage format

The matrix $H$ is stored as a compactdense matrix by rows, that is,
the values of the entries of each row in turn are
stored in order within an appropriate real one-dimensional array.
Since $H$ is symmetric, only the lower triangular part (that is the part
$h_{ij}$ for $0 \leq j \leq i \leq n-1$) need be held.
In this case the lower triangle should be stored by rows, that is
component $i \ast i / 2 + j$of the storage array H_val
will hold the value $h_{ij}$ (and, by symmetry, $h_{ji}$)
for $0 \leq j \leq i \leq n-1$.

###  Sparse co-ordinate storage format

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $H$,
its row index i, column index j
and value $h_{ij}$, $0 \leq j \leq i \leq n-1$,are stored as
the $l$-th components of the integer arrays H_row and
H_col and real array H_val, respectively, while the number of nonzeros
is recorded as H_ne = $ne$.
Note that only the entries in the lower triangle should be stored.

###  Sparse row-wise storage format

Again only the nonzero entries are stored, but this time
they are ordered so that those in row i appear directly before those
in row i+1. For the i-th row of $H$ the i-th component of the
integer array H_ptr holds the position of the first entry in this row,
while H_ptr(n) holds the total number of entries plus one.
The column indices j, $0 \leq j \leq i$, and values
$h_{ij}$ of theentries in the i-th row are stored in components
l = H_ptr(i), $\ldots$, H_ptr(i+1)-1 of the
integer array H_col, and real array H_val, respectively.
Note that as before only the entries in the lower triangle should be stored.
For sparse matrices, this scheme almost always requires less storage than
its predecessor.

### symmetric\_matrix_diagonal Diagonal storage format

If $H$ is diagonal (i.e., $H_{ij} = 0$ for all
$0 \leq i \neq j \leq n-1$) only the diagonals entries
$H_{ii}$, $0 \leq i \leq n-1$ need
be stored, and the first n components of the array H_val may be
used for the purpose.

### symmetric\_matrix_scaled_identity Multiples of the identity storage format

If $H$ is a multiple of the identity matrix, (i.e., $H = \alpha I$
where $I$ is the n by n identity matrix and $\alpha$ is a scalar),
it suffices to store $\alpha$ as the first component of H_val.

### symmetric\_matrix_identity The identity matrix format

If $H$ is the identity matrix, no values need be stored.

### symmetric\_matrix_zero The zero matrix format

The same is true if $H$ is the zero matrix.

