# Introduction

## Purpose

This package uses a primal-dual interior-point trust-region method
to solve the **quadratic programming problem**
$\mbox{minimize}\;\; q(x) = \frac{1}{2} x^T H x + g^T x + f $
\n
minimize q(x) := 1/2 x^T H x + g^T x + f
\n
subject to the general linear constraints
$c_i^l\leqa_i^Tx\leq c_i^u, \;\;\; i = 1, \ldots , m,$
\n
 c_i^l \[<=] a_i^Tx \[<=] c_i^u, i = 1, ... , m,
\n
and the simple bound constraints
$x_j^l\leqx_j \leq x_j^u, \;\;\; j = 1, \ldots , n,$
\n
 x_j^l \[<=] x_j \[<=] x_j^u, j = 1, ... , n,
\n
where the $n$ by $n$ symmetric matrix $H$,
the vectors $g$, $a_i$, $c^l$, $c^u$, $x^l$,
$x^u$ and the scalar $f$ are given.
Any of the constraint bounds $c_i^l$, $c_i^u$,
$x_j^l$ and $x_j^u$ may be infinite.
Full advantage is taken of any zero coefficients in the matrix $H$
or the matrix $A$ of vectors $a_i$.

If the matrix $H$ is positive semi-definite, a global
solution is found. However, if $H$ is indefinite,
the procedure may find a (weak second-order) critical point
that is not the global solution to the given problem.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England, and
Philippe L. Toint, University of Namur, Belgium.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

December 1999, C interface January 2022.

## Terminology

The required solution $x$ necessarily satisfies
the primal optimality conditions
$\mbox{(1a) $\hspace{66mm} A x = c\hspace{66mm}$}$
\n
(1a) A x = c
\n
and
$\mbox{(1b) $\hspace{52mm} c^l \leq c \leq c^u, \;\; x^l \leq x \leq x^u,\hspace{52mm}$} $
\n
(1b) c^l \[<=] c \[<=] c^u, x^l \[<=] x \[<=] x^u,
\n
the dual optimality conditions
$\mbox{(2a) $\hspace{58mm} H x + g = A^T y + z\hspace{58mm}$}$
\n
(2a) H x + g = A^T y + z
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
where the vectors $y$ and $z$ are known as the Lagrange multipliers
for2 the general linear constraints, and the dual variables for the bounds,
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

The problem is solved in two phases. The goal of the first "initial
feasible point" phase is to find a strictly interior point which is
primal feasible, that is that {1a} is satisfied. The GALAHAD package
LSQP is used for this purpose, and offers the options of either
accepting the first strictly feasible point found, or preferably of
aiming for the so-called "analytic center" of the feasible region.
Having found such a suitable initial feasible point, the second
"optimality" phase ensures that \req{4.1a} remains satisfied while
iterating to satisfy dual feasibility (2a) and complementary
slackness (3).The optimality phase proceeds by approximately
minimizing a sequence of barrier functions
\[\frac{1}{2} x^T H x + g^T x + f -
 \mu \left[ \sum_{i=1}^{m} \log ( c_{i}-c_{i}^{l} )
 + \sum_{i=1}^{m} \log ( c_{i}^{u}-c_{i} )
 + \sum_{j=1}^{n} \log ( x_{j}-x_{j}^{l} )
 + \sum_{j=1}^{n} \log ( x_{j}^{u}-x_{j} ) \right],\]
$$\frac{1}{2} x^T H x + g^T x + f -
 \mu \left[ \sum_{i=1}^{m} \log ( c_{i}-c_{i}^{l} )
 + \sum_{i=1}^{m} \log ( c_{i}^{u}-c_{i} )
 + \sum_{j=1}^{n} \log ( x_{j}-x_{j}^{l} )
 + \sum_{j=1}^{n} \log ( x_{j}^{u}-x_{j} ) \right],$$
\n
1/2 x^T H x + g^T x + f -
 mu [ sum_{i=1}^m log (c_i-c_i^l)+ sum_{i=1}^m log (c_i^u-c_i ) +
sum_{j=1}^n log (x_j-x_j^l ) + sum_{j=1}^n log (x_j^u-x_j ) ]
\n
for an approriate sequence of positive barrier parameters $\mu$
converging to zero
while ensuring that (1a) remain satisfied and that
$x$ and $c$ are strictly interior points for (1b).
Note that terms in the above sumations corresponding to infinite bounds are
ignored, and that equality constraints are treated specially.

Each of the barrier subproblems is solved using a trust-region method.
Such a method generates a trial correction step $\Delta (x, c)$
to the current iterate $(x, c)$ by replacing the nonlinear
barrier function locally by a suitable quadratic model, and
approximately minimizing this model in the intersection of \req{4.1a}
and a trust region $\|\Delta (x, c)\| \leq \Delta$ for some
appropriate strictly positive trust-region radius $\Delta$ and norm
$\| \cdot \|$.The step is accepted/rejected and the radius adjusted
on the basis of how accurately the model reproduces the value of
barrier function at the trial step. If the step proves to be
unacceptable, a linesearch is performed along the step to obtain an
acceptable new iterate. In practice, the natural primal "Newton" model
of the barrier function is frequently less successful than an
alternative primal-dual model, and consequently the primal-dual model
is usually to be preferred.

Once a barrier subproblem has been solved, extrapolation based on
values and derivatives encountered on the central path is optionally
used to determine a good starting point for the next subproblem.
Traditional Taylor-series extrapolation has been superceded by more
accurate Puiseux-series methods as these are particularly suited to
deal with degeneracy.

The trust-region subproblem is approximately solved using the combined
conjugate-gradient/Lanczos method implemented in the GALAHAD package
GLTR.Such a method requires a suitable preconditioner, and in our
case, the only flexibility we have is in approximating the model of
the Hessian. Although using a fixed form of preconditioning is
sometimes effective, we have provided the option of an automatic
choice, that aims to balance the cost of applying the preconditioner
against the needs for an accurate solution of the trust-region
subproblem.The preconditioner is applied using the GALAHAD matrix
factorization package SBLS, but options at this stage are to factorize
the preconditioner as a whole (the so-called "augmented system"
approach), or to perform a block elimination first (the
"Schur-complement" approach). The latter is usually to be prefered
when a (non-singular) diagonal preconditioner is used, but may be
inefficient if any of the columns of $A$ is too dense.

In order to make the solution as efficient as possible, the variables
and constraints are reordered internally by the GALAHAD package QPP
prior to solution.In particular, fixed variables, and free
(unbounded on both sides) constraints are temporarily removed.

## Reference

The basic algorithm is a generalisation of those of

Y. Zhang (1994),
 On the convergence of a class of infeasible interior-point methods for the
 horizontal linear complementarity problem,
 SIAM J. Optimization 4(1) 208-227,

with a number of enhancements described by

A. R. Conn, N. I. M. Gould, D. Orban and Ph. L. Toint (1999).
A primal-dual trust-region algorithm for minimizing a non-convex
function subject to general inequality and linear equality constraints.
Mathematical Programming **87** 215-249.

## Call order

To solve a given problem, functions from the qpb package must be called
in the following order:

- qpb\_initialize - provide default control parameters and
set up initial data structures
- qpb\_read\_specfile (optional) - override control values
by reading replacement values from a file
- qpb\_import - set up problem data structures and fixed
values
- qpb\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- qpb\_solve_qp - solve the quadratic program
- qpb\_information (optional) - recover information about
the solution and solution process
- qpb\_terminate - deallocate data structures

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

