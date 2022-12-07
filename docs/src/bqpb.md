# Introduction

## Purpose

This package uses a primal-dual interior-point method
to solve the **convex quadratic programming problem**
$\mbox{minimize}\;\; q(x) = \frac{1}{2} x^T H x + g^T x + f $
\[
minimize q(x) := 1/2 x^T H x + g^T x + f
\]
or the **shifted least-distance problem**
$\mbox{minimize}\;\; \frac{1}{2} \sum_{j=1}^n w_j^2 ( x_j - x_j^0 )^2
 + g^T x + f $
\n
 minimize 1/2 \sum_{j=1}^n w_j^2 ( x_j - x_j^0 )^2+ g^T x + f
\n
subject to the simple bound constraints
$x_j^l\leqx_j \leq x_j^u, \;\;\; j = 1, \ldots , n,$
\n
 x_j^l \[<=] x_j \[<=] x_j^u, j = 1, ... , n,
\n
where the $n$ by $n$ symmetric, positive-semi-definite matrix
$H$, the vectors $g$, $w$, $x^{0}$,
$x^l$,$x^u$ and the scalar $f$ are given. Any of
the constraint bounds $x_j^l$ and $x_j^u$ may be infinite.
Full advantage is taken of any zero coefficients in the matrix $H$.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

July 2021, C interface December 2021.

## Method

The required solution $x$ necessarily satisfies
the primal optimality conditions
$\mbox{(1) $\hspace{52mm} x^l \leq x \leq x^u,\hspace{52mm}$} $
\n
(1) x^l \[<=] x \[<=] x^u,
\n
the dual optimality conditions
$\mbox{(2a) $\hspace{3mm} H x + g = z \;\; (\mbox{or}
\;\;W^{2} (x -x^0) + g = z \;\; \mbox{for the shifted-least-distance type objective})$}$
\n
(2a) H x + g = z
 (or W^2 (x -x^0) + g = z for the shifted-least-distance type objective)
\n
where
$\mbox{(2b) $\hspace{24mm} z = z^l + z^u, \,\,
 z^l \geq 0 \;\; \mbox{and} \;\; z^u \leq 0,\hspace{24mm}$} $
\n
 (2b) z = z^l + z^u, z^l \[>=] 0 and z^u \[<=] 0,
\n
and the complementary slackness conditions
$\mbox{(3) $\hspace{12mm}
(x -x^l )^{T} z^l = 0 \;\;\mbox{and} \;\; (x -x^u )^{T} z^u = 0,\hspace{12mm} $}$
\n
(3) (x -x^l)^T z^l = 0 and (x -x^u)^T z^u = 0,
\n
where the diagonal matrix $W^2$ has diagonal entries $w_j^2$,
$j = 1, \ldots , n$, where the vector $z$ is known as
the dual variables for the bounds,
respectively, and where the vector inequalities hold component-wise.

Primal-dual interior point methods iterate towards a point
that satisfies these conditions by ultimately aiming to satisfy
(2a) and (3), while ensuring that (1) and (2b) are
satisfied as strict inequalities at each stage.Appropriate norms of the
amounts bywhich (2a) and (3) fail to be satisfied are known as the
primal and dual infeasibility, and the violation of complementary slackness,
respectively. The fact that (1) and (2b) are satisfied as strict
inequalities gives such methods their other title, namely
interior-point methods.

The method aims at each stage to reduce the
overall violation of (2a) and (3),
rather than reducing each of the terms individually. Given an estimate
$v = (x, c, z, z^l, z^u)$
of the primal-dual variables, a correction
$\Delta v = \Delta (x, c, z, z^l, z^u)$
is obtained by solving a suitable linear system of Newton equations for the
nonlinear systems (2a) and a parameterized “residual
trajectory” perturbation of (3); residual trajectories
proposed by Zhang (1994) and Zhao and Sun (1999) are possibilities.
An improved estimate $v + \alpha \Delta v$
is then used, where the step-size $\alpha$
is chosen as close to 1.0 as possible while ensuring both that
(1) and (2b) continue to hold and that the individual components
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

The package is actually just a front-end to the more-sophisticated
GALAHAD package CQP that saves users from setting unnecessary arguments.

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

To solve a given problem, functions from the bqpb package must be called
in the following order:

- bqpb\_initialize - provide default control parameters and set up initial data structures
- bqpb\_read\_specfile (optional) - override control values by reading replacement values from a file
- bqpb\_import - set up problem data structures and fixed values
- bqpb\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- solve the problem by calling one of
- bqpb\_solve_qp - solve the bound-constrained
quadratic program
- bqpb\_solve_sldqp - solve the bound-constrained
 shifted least-distance problem
- bqpb\_information (optional) - recover information about the solution and solution process
- bqpb\_terminate - deallocate data structures

##  Symmetric matrix storage formats

The symmetric $n$ by $n$ objective Hessian matrix $H$ may be
presented and stored in a variety of convenient formats. But crucially
symmetry is exploitedby only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal).

Both C-style (0 based)and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; the discussion below presumes C style, but add 1 to
indices for the corresponding fortran version.

Wrappers will automatically convert between 0-based (C) and 1-based
(fortran) array indexing, so may be used transparently from C. This
conversion involves both time and memory overheads that may be avoided
by supplying data that is already stored using 1-based indexing.

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


