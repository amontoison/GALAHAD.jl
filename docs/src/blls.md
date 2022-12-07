# Introduction

## Purpose

This package uses a preconditioned, projected-gradient method to solve the
 **bound-constrained regularized linear least-squares problem**
$\mbox{minimize}\;\; q(x) = q(x) = \frac{1}{2} \| A x - b\|_2^2 + \frac{1}{2} \sigma \|x\|^2$
\n
minimize q(x) := 1/2 || A x - b ||^2 + sigma ||x||^2
\n
subject to the simple bound constraints
$x_j^l\leqx_j \leq x_j^u, \;\;\; j = 1, \ldots , n,$
\n
 x_j^l \[<=] x_j \[<=] x_j^u, j = 1, ... , n,
\n
where the $m$ by $n$ real matrix $A$, the vectors
$b$, $x^{l}$, $x^{u}$ and the non-negative weight
$\sigma$ are given. Any of the constraint bounds $x_j^l$ and
$x_j^u$ may be infinite.Full advantage is taken of any zero
coefficients of the Jacobian matrix $A$ of the **residuals**
$c(x) = A x - b$;the matrix need not be provided as there are options
to obtain matrix-vector products involving $A$ and its transpose either
by reverse communication or from a user-provided subroutine.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

October 2019, C interface March 2022.

## Terminology

The required solution $x$ necessarily satisfies
the primal optimality conditions
$x^l \leq x \leq x^u,$
\n
 x^l \[<=] x \[<=] x^u,
\n
the dual optimality conditions
$(A^T A + \sigma I ) x = A^T b + z$
\n
 ( A^T A + sigma I ) x = A^T b + z
\n
where
$ z = z^l + z^u, \,\,
 z^l \geq 0 \;\; \mbox{and} \;\; z^u \leq 0,$
\n
 z = z^l + z^u, z^l \[>=] 0 and z^u \[<=] 0,
\n
and the complementary slackness conditions
$(x -x^l )^{T} z^l = 0 \;\;\mbox{and} \;\; (x -x^u )^{T} z^u = 0,\hspace{12mm} $
\n
(x -x^l)^T z^l = 0 and (x -x^u)^T z^u = 0,
\n
where the vector $z$ is known asthe dual variables for the bounds,
respectively, and where the vector inequalities hold component-wise.

## Method

The method is iterative. Each iteration proceeds in two stages.
Firstly, a search direction $s$ from the current estimate of the solution
$x$ is computed. This may be in a scaled steepest-descent direction, or,
if the working set of variables on bounds has not changed dramatically,
in a direction that provides an approximate minimizer of the objective
over a subspace comprising the currently free-variables. The latter is
computed either using an appropriate sparse factorization by the
GALAHAD package SBLS, or by theconjugate-gradient least-squares (CGLS)
method; tt may be necessary to regularize the subproblem very slightly to
avoid a ill-posedness. Thereafter, a piecewise linesearch (arc search) is
carried out along the arc $x(\alpha) = P( x + \alpha s)$ for
$\alpha > 0$, where the projection operator is defined component-wise
at any feasible point $v$ to be
$P_j(v) = \min( \max( x_j, x_j^l), x_j^u);$
thus this arc bends the search direction into the feasible region.
The arc search is performed either exactly, by passing through a set
of increasing breakpoints at which it changes direction, or inexactly,
by evaluating a sequence of different $\alpha$on the arc.
All computation is designed to exploit sparsity in $A$.

## Reference

Full details are provided in

N. I. M. Gould (2022).
Numerical methods for solving bound-constrained linear least squares problems.
In preparation.

## Call order

To solve a given problem, functions from the blls package must be called
in the following order:

- blls\_initialize - provide default control parameters and
set up initial data structures
- blls\_read\_specfile (optional) - override control values
by reading replacement values from a file
- set up problem data structures and fixed values by caling one of
 - blls\_import - in the case that $A$ is explicitly
available
 - blls\_import\_without_a - in the case that only the
effect of applying $A$ and its transpose to a vector is possible
- blls\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- solve the problem by calling one of
 - blls\_solve_given_a - solve the problem using values
of $A$
 - blls\_solve\_reverse_a_prod - solve the problem by returning
 to the caller for products of $A$ and its transpose with specified
vectors
- blls\_information (optional) - recover information about
the solution and solution process
- blls\_terminate - deallocate data structures

##  Unsymmetric matrix storage formats

The unsymmetric $m$ by $n$ matrix $A$ may be presented
and stored in a variety of convenient input formats.

Both C-style (0 based)and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; the discussion below presumes C style, but add 1 to
indices for the corresponding fortran version.

Wrappers will automatically convert between 0-based (C) and 1-based
(fortran) array indexing, so may be used transparently from C. This
conversion involves both time and memory overheads that may be avoided
by supplying data that is already stored using 1-based indexing.

### unsymmetric\_matrix_dense_row Dense row storage format

The matrix $A$ is stored as a compactdense matrix by rows, that is,
the values of the entries of each row in turn are
stored in order within an appropriate real one-dimensional array.
In this case, component $n \ast i + j$of the storage array A_val
will hold the value $A_{ij}$ for $0 \leq i \leq m-1$,
$0 \leq j \leq n-1$.

### unsymmetric\_matrix_dense_column Dense column storage format

The matrix $A$ is stored as a compactdense matrix by columns, that is,
the values of the entries of each column in turn are
stored in order within an appropriate real one-dimensional array.
In this case, component $m \ast j + i$of the storage array A_val
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
its predecessors.

### unsymmetric\_matrix_column_wise Sparse column-wise storage format

Again only the nonzero entries are stored, but this time
they are ordered so that those in column j appear directly before those
in column j+1. For the j-th column of $A$ the j-th component of the
integer array A_ptr holds the position of the first entry in this column,
while A_ptr(n) holds the total number of entries plus one.
The row indices i, $0 \leq i \leq m-1$, and values $A_{ij}$
of thenonzero entries in the j-th column are stored in components
l = A_ptr(j), $\ldots$, A_ptr(j+1)-1,$0 \leq j \leq n-1$,
of the integer array A_row, and real array A_val, respectively.
Once again, for sparse matrices, this scheme almost always requires less
storage than the dense of coordinate formats.

