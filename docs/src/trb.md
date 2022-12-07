# Introduction

## Purpose

The trb package uses a **trust-region method to find a (local)
 minimizer of a differentiable objective function
$\mathbf{f(x)}$ of many variables $\mathbf{x}$,
where the variables satisfy the simple bounds
$\mathbf{x^l \leq x \leq x^u}$.**
The method offers the choice of
direct and iterative solution of the key subproblems, and
is most suitable for large problems. First derivatives are required,
and if second derivatives can be calculated, they will be exploited---if
the product of second derivatives with a vector may be found but
not the derivatives themselves, that may also be exploited.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

July 2021, C interface August 2021.

## Terminology

The gradient $\nabla_x f(x)$ of $f(x)$ is the vector whose
$i$-th component is $\partial f(x)/\partial x_i$.
The Hessian $\nabla_{xx} f(x)$ of $f(x)$ is the symmetric matrix
whose $i,j$-th entry is $\partial^2 f(x)/\partial x_i \partial x_j$.
The Hessian is sparse if a significant and useful proportion of the
entries are universally zero.

## Method

A trust-region method is used. In this, an improvement to a current estimate
of the required minimizer, $x_k$ is sought by computing a step $s_k$.
The step is chosen to approximately minimize a model $m_k(s)$ of
 $f(x_k + s)$ within the intersection of the boundconstraints
$x^l \leq x \leq x^u$ and a trust region $\|s_k\| \leq \Delta_k$
for some specified positive "radius" $\Delta_k$. The quality of the
resulting step $s_k$ is assessed by computing the "ratio"
$(f(x_k) - f(x_k + s_k))/ (m_k(0) - m_k(s_k))$. The step is deemed
to have succeeded if the ratio exceeds a given $\eta_s > 0$,
and in this case $x_{k+1} = x_k + s_k$. Otherwise
$x_{k+1} = x_k$, and the radius is reduced by powers of a given
reduction factor until it is smaller than $\|s_k\|$. If the ratio is
larger than$\eta_v \geq \eta_d$, the radius will be increased so that
it exceeds $\|s_k\|$ by a given increase factor. The method will
terminate as soon as $\|\nabla_x f(x_k)\|$ is smaller than a
specified value.

Either linear or quadratic models $m_k(s)$ may be used. The former will
be taken as the first two terms $f(x_k) + s^T \nabla_x f(x_k)$
of a Taylor series about $x_k$, while the latter uses an
approximation to the first three terms
$f(x_k) + s^T \nabla_x f(x_k) + \frac{1}{2} s^T B_k s$,
for which $B_k$ is a symmetric approximation to the Hessian
$\nabla_{xx}f(x_k)$; possible approximations include the true Hessian,
limited-memory secant and sparsity approximations and a scaled identity
matrix. Normally a two-norm trust region will be used, but this may change
if preconditioning is employed.

The model minimization is carried out in two stages.
Firstly, the so-called generalized Cauchy point for the quadratic
subproblem is found---the purpose of this point is to ensure that the
algorithm converges and that the set of bounds which are satisfied as
equations at the solution is rapidly identified.Thereafter an
improvement to the quadratic model on the face of variables predicted
to be active by the Cauchy point is sought using either a
direct approach involving factorization or an
iterative (conjugate-gradient/Lanczos) approach based on approximations
to the required solution from a so-called Krlov subspace. The direct
approach is based on the knowledge that the required solution
satisfies the linear system of equations
$(B_k + \lambda_k I) s_k= - \nabla_x f(x_k)$, involving a scalar
Lagrange multiplier $\lambda_k$, on the space of inactive variables.
This multiplier is found by uni-variate root finding, using a safeguarded
Newton-like process, by the GALAHAD package TRS. The iterative approach
uses GALAHAD package GLTR, and is best accelerated by preconditioning
with good approximations to $B_k$ using GALAHAD's PSLS. The
iterative approach has the advantage that only matrix-vector products
$B_k v$ are required, and thus $B_k$ is not required explicitly.
However when factorizations of $B_k$ are possible, the direct approach
is often more efficient.

The iteration is terminated as soon as the Euclidean norm of the
projected gradient,
$\|\min(\max( x_k - \nabla_x f(x_k), x^l), x^u) -x_k\|_2,$
is sufficiently small. At such a point, $\nabla_x f(x_k) = z_k$,
where the $i$-th dual variable $z_i$ is non-negative if
$x_i$ is on its lower bound $x^l_i$, non-positive if $x_i$
is on its upper bound $x^u_i$, and zero if $x_i$ lies strictly
between its bounds.

## References

The generic bound-constrained trust-region method is described in detail in

A. R. Conn, N. I. M. Gould and Ph. L. Toint,
"Trust-region methods",
SIAM/MPS Series on Optimization (2000).

# Call order

To solve a given problem, functions from the trb package must be called
in the following order:

- trb\_initialize - provide default control parameters and
set up initial data structures
- trb\_read\_specfile (optional) - override control values
by reading replacement values from a file
- trb\_import - set up problem data structures and fixed
values
- trb\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- solve the problem by calling one of
 - trb\_solve\_with\_mat - solve using function calls to
 evaluate function, gradient and Hessian values
 - trb\_solve\_without\_mat - solve using function calls to
 evaluate function and gradient values and Hessian-vector products
 - trb\_solve\_reverse\_with\_mat - solve returning to the
 calling program to obtain function, gradient and Hessian values, or
 - trb\_solve\_reverse\_without\_mat - solve returning to the
 calling prorgram to obtain function and gradient values and
 Hessian-vector products
- trb\_information (optional) - recover information about
the solution and solution process
- trb\_terminate - deallocate data structures

#  Symmetric matrix storage formats

The symmetric $n$ by $n$ matrix $H = \nabla_{xx}f$ may be
presented and stored in a variety of formats. But crucially symmetry
is exploited by only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal).

Both C-style (0 based)and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; the discussion below presumes C style, but add 1 to
indices for the corresponding fortran version.

Wrappers will automatically convert between 0-based (C) and 1-based
(fortran) array indexing, so may be used transparently from C. This
conversion involves both time and memory overheads that may be avoided
by supplying data that is already stored using 1-based indexing.

## Dense storage format

The matrix $H$ is stored as a compactdense matrix by rows, that is,
the values of the entries of each row in turn are
stored in order within an appropriate real one-dimensional array.
Since $H$ is symmetric, only the lower triangular part (that is the part
$H_{ij}$ for $0 \leq j \leq i \leq n-1$) need be held.
In this case the lower triangle should be stored by rows, that is
component $i \ast i / 2 + j$of the storage array H_val
will hold the value $H_{ij}$ (and, by symmetry, $H_{ji}$)
for $0 \leq j \leq i \leq n-1$.

##  Sparse co-ordinate storage format

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $H$,
its row index i, column index j
and value $H_{ij}$, $0 \leq j \leq i \leq n-1$,are stored as
the $l$-th components of the integer arrays H_row and
H_col and real array H_val, respectively, while the number of nonzeros
is recorded as H_ne = $ne$.
Note that only the entries in the lower triangle should be stored.

##  Sparse row-wise storage format

Again only the nonzero entries are stored, but this time
they are ordered so that those in row i appear directly before those
in row i+1. For the i-th row of $H$ the i-th component of the
integer array H_ptr holds the position of the first entry in this row,
while H_ptr(n) holds the total number of entries plus one.
The column indices j, $0 \leq j \leq i$, and values
$H_{ij}$ of theentries in the i-th row are stored in components
l = H_ptr(i), $\ldots$, H_ptr(i+1)-1 of the
integer array H_col, and real array H_val, respectively.
Note that as before only the entries in the lower triangle should be stored.
For sparse matrices, this scheme almost always requires less storage than
its predecessor.
