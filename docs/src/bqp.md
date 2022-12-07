# Introduction

## Purpose

This package usesa preconditioned, projected-gradient method
to solve the **convex bound-constrained quadratic programming problem**
$\mbox{minimize}\;\; q(x) = \frac{1}{2} x^T H x + g^T x + f $
\[
minimize q(x) := 1/2 x^T H x + g^T x + f
\]
subject to the simple bound constraints
$x_j^l\leqx_j \leq x_j^u, \;\;\; j = 1, \ldots , n,$
\n
 x_j^l \[<=] x_j \[<=] x_j^u, j = 1, ... , n,
\n
where the $n$ by $n$ symmetric postive semi-definite matrix $H$,
the vectors $g$, $x^l$, $x^u$ and the scalar $f$ are given.
Any of the constraint bounds
$x_j^l$ and $x_j^u$ may be infinite.
Full advantage is taken of any zero coefficients in the matrix $H$;
the matrix need not be provided as there are options to obtain matrix-vector
products involving $H$ by reverse communication.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

November 2009, C interface February 2022.

## Terminology

The required solution $x$ necessarily satisfies
the primal optimality conditions
$x^l \leq x \leq x^u,$
\n
 x^l \[<=] x \[<=] x^u,
\n
the dual optimality conditions
$H x + g = z$
\n
 H x + g = z
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
Firstly, the so-called generalized Cauchy point for the quadratic
objective is found.(The purpose of this point is to ensure that the
algorithm converges and that the set of bounds which are satisfied as
equations at the solution is rapidly identified.)Thereafter an
improvement to the objective is sought using either a
direct-matrix or truncated conjugate-gradient algorithm.

## Reference

This is a specialised version of the method presented in

A. R. Conn, N. I. M. Gould and Ph. L. Toint (1988).
Global convergence of a class of trust region algorithms
for optimization with simple bounds.
SIAM Journal on Numerical Analysis **25** 433-460,

## Call order

To solve a given problem, functions from the bqp package must be called
in the following order:

- bqp\_initialize - provide default control parameters and set up initial data structures
- bqp\_read\_specfile (optional) - override control values by reading replacement values from a file
- set up problem data structures and fixed values by caling one of
- bqp\_import - in the case that $H$ is explicitly
available
- bqp\_import\_without_h - in the case that only the
effect of applying $H$ to a vector is possible
- bqp\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- solve the problem by calling one of
- bqp\_solve_given_h - solve the problem using values
of $H$
- bqp\_solve\_reverse_h_prod - solve the problem by returning
 to the caller for products of $H$ with specified vectors
- bqp\_information (optional) - recover information about the solution and solution process
- bqp\_terminate - deallocate data structures

##  Symmetric matrix storage formats

If it is explicitly available, the symmetric $n$ by $n$
objective Hessian matrix $H$ may be presented
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
