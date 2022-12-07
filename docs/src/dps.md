# Introduction

## Purpose

Given a real $n$ by $n$ symmetric matrix $H$, this package
**construct a
symmetric, positive definite matrix $M$ so that $H$
is diagonal in the norm $\|v\|_{M} = \sqrt{v^T M v}$
induced by $M$**. Subsequently the package can be use to
**solve the trust-region subproblem**
$\mbox{(1)}\;\; \mbox{minimize}\;\; q(x) = \frac{1}{2} x^T H x + c^T x
+ f \;\; \mbox{subject to}\;\; \|x\||_{M} \leq \Delta$
or the **regularized quadratic problem**
$\mbox{(2)}\;\;\mbox{minimize}\;\; q(x) + \frac{1}{p} \sigma \|x\||_{M}^p\hspace{50mm} \mbox{$$}$
for a real $n$ vector $c$ and scalars $f$,
$\Delta>0$, $\sigma>0$ and $p \geq 2$.

A factorization of the matrix $H$ will be required, so this package is
most suited for the case where such a factorization, either dense or sparse,
may be found efficiently.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

August 2011, C interface December 2021.

## Terminology

## Method

The required solution $x_*$ necessarily satisfies the optimality
condition $H x_* + \lambda_* M x_* + c = 0$,
where $\lambda_* \geq 0$ is a Lagrange
multiplier that corresponds to the constraint
$\|x\|_{M}\leq\Delta$ in the trust-region case (1),
and is given by $\lambda_* = \sigma \|x_*\|^{p-2}$
for the regularization problem (2).
In addition $H + \lambda_* M$ will be positive semi-definite; in
most instances it will actually be positive definite, but in special
“hard” cases singularity is a possibility.

The matrix $H$ is decomposed as
$H = P L D L^T P^T$
by calling the GALAHAD package SLS.
Here $P$ is a permutation matrix,
$L$ is unit lower triangular and $D$ is block diagonal, with
blocks of dimension at most two. The spectral decomposition of each diagonal
block of $D$ is computed, and each eigenvalue $\theta$ is replaced by
$\max ( | \theta | , \theta_{\min} ) $,
where $\theta_{\min}$ is a positive user-supplied value. The resulting
block diagonal matrix is $B$, from which we define the
**modified-absolute-value**
$M = P L B L^T P^T;$
an alternative due to Goldfarb uses instead the simpler
$M = P L L^T P^T.$

Given the factors of $H$ (and $M$), the required solution is
found by making the change of variables $y = B^{1/2} L^T P^T x$
(or $y = L^T P^T x$ in the Goldfarb case)
which results in “diagonal” trust-region and regularization subproblems,
whose solution may be easily obtained suing a Newton or higher-order iteration
of a resulting “secular” equation.If subsequent problems, for which
$H$ and $c$ are unchanged, are to be attempted, the existing
factorization and solution may easily be exploited.

The dominant cost is that for the factorization of the symmetric, but
potentially indefinite, matrix $H$ using the GALAHAD package SLS.

## Reference

The method is described in detail for the trust-region case in

N. I. M. Gould and J. Nocedal (1998).
The modified absolute-value factorization for trust-region minimization.
In “High Performance Algorithms and Software in Nonlinear Optimization”
(R. De Leone, A. Murli, P. M. Pardalos and G. Toraldo, eds.),
Kluwer Academic Publishers, pp. 225-241,

while the adaptation for the regularization case is obvious. The method used
to solve the diagonal trust-region and regularization subproblems are as
given by

H. S. Dollar, N. I. M. Gould and D. P. Robinson (2010).
On solving trust-region and other regularised subproblems in optimization.
Mathematical Programming Computation **2(1)** 21-57

with simplifications due to the diagonal Hessian.

## Call order

To solve a given problem, functions from the dps package must be called
in the following order:

- dps\_initialize - provide default control parameters and
set up initial data structures
- dps\_read\_specfile (optional) - override control values
by reading replacement values from a file
- dps\_import - import control and matrix data structures
- dps\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- one of
- dps\_solve_tr_problem - solve the trust-region problem (1)
- dps\_solve_rq_problem - solve the regularized-quadratic
problem (2)
- optionally one of
- dps_resolve_tr_problem - resolve the trust-region problem
(1) when the non-matrix data has changed
- dps_resolve_rq_problem - resolve the regularized-quadratic
problem (2) when the non-matrix data has changed
- dps\_information (optional) - recover information about
the solution and solution process
- dps\_terminate - deallocate data structures

##  Symmetric matrix storage formats

The symmetric $n$ by $n$ coefficient matrix $H$ may be presented
and stored in a variety of convenient input formats.Crucially symmetry
is exploitedby only storing values from the lower triangular part
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
$H_{ij}$ for $0 \leq j \leq i \leq n-1$) need be held.
In this case the lower triangle should be stored by rows, that is
component $i \ast i / 2 + j$of the storage array val
will hold the value $H_{ij}$ (and, by symmetry, $H_{ji}$)
for $0 \leq j \leq i \leq n-1$.

###  Sparse co-ordinate storage format

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $H$,
its row index i, column index j
and value $H_{ij}$, $0 \leq j \leq i \leq n-1$,are stored as
the $l$-th components of the integer arrays row and
col and real array val, respectively, while the number of nonzeros
is recorded as ne = $ne$.
Note that only the entries in the lower triangle should be stored.

###  Sparse row-wise storage format

Again only the nonzero entries are stored, but this time
they are ordered so that those in row i appear directly before those
in row i+1. For the i-th row of $H$ the i-th component of the
integer array ptr holds the position of the first entry in this row,
while ptr(n) holds the total number of entries plus one.
The column indices j, $0 \leq j \leq i$, and values
$H_{ij}$ of theentries in the i-th row are stored in components
l = ptr(i), $\ldots$, ptr(i+1)-1 of the
integer array col, and real array val, respectively.
Note that as before only the entries in the lower triangle should be stored.
For sparse matrices, this scheme almost always requires less storage than
its predecessor.
