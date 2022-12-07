# Introduction

## Purpose

This package uses an iterative method to solve the
**equality-constrained quadratic programming problem**
$\mbox{minimize}\;\; q(x) = \frac{1}{2} x^T H x + g^T x + f $
\[
minimize q(x) := 1/2 x^T H x + g^T x + f
\]
subject to the linear constraints
$(1) \;\; A x + c = 0,$
where the $n$ by $n$ symmetric matrix $H$,
the $m$ by $n$ matrix $A$, the vectors $g$ and $c$
Full advantage is taken of any zero coefficients in the matrices $H$
and $A$.

The package may alternatively be used to minimize the (shifted) squared-
least-distance objective
$\frac{1}{2} \sum_{j=1}^n w_j^2 ( x_j - x_j^0 )^2 + g^T x + f,$
\n
 minimize 1/2 \sum_{j=1}^n w_j^2 ( x_j - x_j^0 )^2+ g^T x + f,
\n
subject to the linear constraint (1), for given vectors $w$ and
$x^0$.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

March 2006, C interface January 2021.

## Terminology

The required solution $x$ necessarily satisfies
the primal optimality conditions
$(2) \;\; A x + c = 0$
\n
(2) A x + c = 0
\n
and the dual optimality conditions
\[H x + g - A^T y = 0 \;\; (\mbox{or} \;\; W^{2} (x -x^0) + g - A^T y = 0 \;\;\mbox{for the shifted-least-distance type objective})\]
$$H x + g - A^T y = 0 \;\; (\mbox{or} \;\;W^{2} (x -x^0) + g - A^T y = 0 \;\; \mbox{for the shifted-least-distance type objective})$$
\n
(3) H x + g - A^T y = 0
 (or W^2 (x -x^0) + g - A^T y = 0
for the shifted-least-distance type objective)
\n
where the diagonal matrix $W^2$ has diagonal entries $w_j^2$,
$j = 1, \ldots , n$, and where the vector $y$ is
known as the Lagrange multipliers for the linear constraints.

 ## Method

A solution to the problem is found in two phases.
In the first, a point $x_F$ satisfying (2) is found.
In the second, the required solution $x = x_F + s$
is determined by finding $s$ to minimize
$q(s) = \frac{1}{2} s^T H s + g_F^T s + f_F$
subject to the homogeneous constraints $A s = zero$,
where $g_F = H x_F + g$ and
$f_F = \frac{1}{2} x_F^T H x_F + g^T x_F + f$.
The required constrained minimizer of $q(s)$ is obtained
by implictly applying the preconditioned conjugate-gradient method
in the null space of $A$. Any preconditioner of the form
$ K_G = \mat{cc}{ G & A^T \\ A& 0 }$
\n
K_G = ( GA^T )
( A 0)
\n
is suitable, and the GALAHAD package SBLS provides a number of
possibilities. In order to ensure that the minimizer obtained is
finite, an additional, precautionary trust-region constraint $\|s\|
\leq \Delta$ for some suitable positive radius $\Delta$ is
imposed, and the GALAHAD package GLTR is used to solve this
additionally-constrained problem.

## Reference

The preconditioning aspcets are described in detail in

H. S. Dollar, N. I. M. Gould and A. J. Wathen.
“On implicit-factorization constraint preconditioners”.
InLarge Scale Nonlinear Optimization (G. Di Pillo and M. Roma, eds.)
Springer Series on Nonconvex Optimization and Its Applications, Vol. 83,
Springer Verlag (2006) 61-82

and

H. S. Dollar, N. I. M. Gould, W. H. A. Schilders and A. J. Wathen
“On iterative methods and implicit-factorization preconditioners for
regularized saddle-point systems”.
SIAM Journal on Matrix Analysis and Applications, **28(1)** (2006)
170-189,

while the constrained conjugate-gradient method is discussed in

N. I. M. Gould, S. Lucidi, M. Roma and Ph. L. Toint,
Solving the trust-region subproblem using the Lanczos method.
SIAM Journal on Optimization **9:2** (1999), 504-525.

## Call order

To solve a given problem, functions from the eqp package must be called
in the following order:

- eqp\_initialize - provide default control parameters and set up initial data structures
- eqp\_read\_specfile (optional) - override control values by reading replacement values from a file
- eqp\_import - set up problem data structures and fixed values
- eqp\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- solve the problem by calling one of
- eqp\_solve_qp - solve the quadratic program
- eqp\_solve_sldqp - solve the shifted least-distance problem
- eqp_resolve_qp (optional) - resolve the problem with the
same Hessian and Jacobian, but different $g$, $f$ and/or $c$
- eqp\_information (optional) - recover information about the solution and solution process
- eqp\_terminate - deallocate data structures

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
