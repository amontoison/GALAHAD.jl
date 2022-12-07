# Introduction

## Purpose

Given real $n$ by $n$ symmetric matrices $H$ and $M$
(with $M$ diagonally dominant), another real $m$ by $n$
matrix $A$, a real $n$ vector $c$ and scalars
$\Delta>0$ and $f$, this package finds a
**global minimizer of the quadratic objective function
 $\frac{1}{2} x^T Hx + c^T x + f$,
where the vector $x$ is required to satisfy
the constraint $\|x\|_M \leq \Delta$ and possibly $A x =0$,**
and where the $M$-norm of $x$ is $\|x\|_M = \sqrt{x^T M x}$.
This problem commonly occurs as atrust-region subproblem in nonlinear
optimization calculations.The package may also be used to solve the
related problem in which $x$ is
instead required to satisfy the **equality constraint
$\|x\|_M = \Delta$**.The matrix $M$ need not be provided in
the commonly-occurring $\ell_2$-trust-region case for which
$M = I$, the $n$ by $n$ identity matrix.

Factorization of matrices of the form $H + \lambda M$---or
$\mbox{(1)}\;\;\; \mat{cc}{ H + \lambda M & A^T \\ A & 0}$
\n
(1) ( H + lambda M A^T )
(A0)
\n
in cases where $A x = 0$ is imposed---for a succession of
scalars $\lambda$ will be required, so this package is most suited for
the case where such a factorization may be found efficiently. If this is
not the case, the GALAHAD package GLTR may be preferred.

## Authors

N. I. M. Gould and H. S. Thorne, STFC-Rutherford Appleton Laboratory, England,
and D. P. Robinson, Oxford University, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

October 2008, C interface December 2021.

## Method

The method is iterative, and proceeds in two phases.Firstly,
lower and upper bounds, $\lambda_L$ and
$\lambda_U$, on $\lambda_*$ are computed using
Gershgorin's theorems and other eigenvalue bounds. The first phase of
the computation proceeds by progressively shrinking the bound interval
$[\lambda_L,\lambda_U]$
until a value $\lambda$ for which $\|x(\lambda)\|_M \geq \Delta$
is found.Here $x(\lambda)$ and its companion $y(\lambda)$ are
defined to be a solution of
$\mbox{(2)}\;\;\; (H + \lambda M)x(\lambda)
 + A^T y(\lambda) = - c \;\mbox{and}\; A x(\lambda) = 0.$
\n
 (2)(H + lambda M)x(lambda) + A^T y(lambda) = - c and
A x(lambda) = 0;
\n
along the way the possibility that $H$ might be positive definite on
the null-space of $A$ and $\|x(0)\|_M \leq \Delta$ is examined, and
if this transpires the process is terminated with $x_* = x(0)$.
Once the terminating $\lambda$ from the first phase has
been discovered, the second phase consists of applying Newton or
higher-order iterations to the nonlinear “secular” equation
$\|x(\lambda)\|_M = \Delta$ with the knowledge that such
iterations are both globally and ultimately rapidly convergent. It is
possible in the “hard” case that the interval in the first-phase will
shrink to the single point $\lambda_*$, and precautions are taken, using
inverse iteration with Rayleigh-quotient acceleration to ensure that
this too happens rapidly.

The dominant cost is the requirement that we solve a sequence of linear
systems (2). In the absence of linear constraints, an efficient
sparse Cholesky factorization with precautions to detect indefinite
$H + \lambda M$ is used. If $A x = 0$ is required, a sparse
symmetric, indefinite factorization of (1) is used rather than a
Cholesky factorization.

## Reference

The method is described in detail in

H. S. Dollar, N. I. M. Gould and D. P. Robinson.
On solving trust-region and other regularised subproblems in optimization.
Mathematical Programming Computation **2(1)** (2010) 21--57.

## Call order

To solve a given problem, functions from the trs package must be called
in the following order:

- trs\_initialize - provide default control parameters and
set up initial data structures
- trs\_read\_specfile (optional) - override control values
by reading replacement values from a file
- trs\_import - set up problem data structures and fixed
values
- trs\_import_m - (optional) set up problem data structures
and fixed values for the scaling matrix $M$, if any
- trs\_import_a - (optional) set up problem data structures
and fixed values for the constraint matrix $A$, if any
- trs\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- trs\_solve_problem - solve the trust-region problem
- trs\_information (optional) - recover information about
the solution and solution process
- trs\_terminate - deallocate data structures

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
$H$ and scaling matrix $M$ may be presented
and stored in a variety of formats. But crucially symmetry is exploited
by only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal).
In what follows, we refer to $H$ but this applies equally to $M$.

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
