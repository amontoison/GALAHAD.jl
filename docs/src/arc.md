# Introduction

## Purpose

The arc package uses a **regularization method to find a (local)
unconstrained minimizer of a differentiable objective function
$\mathbf{f(x)}$ of many variables $\mathbf{x}$.**
The method offers the choice of
direct and iterative solution of the key regularization subproblems, and
is most suitable for large problems. First derivatives are required,
and if second derivatives can be calculated, they will be exploited---if
the product of second derivatives with a vector may be found, but
not the derivatives themselves, that may also be exploited.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England,
and M. Porcelli, University of Bologna, Italy.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

May 2011, C interface August 2021.

## Terminology

The gradient $\nabla_x f(x)$ of $f(x)$ is the vector whose
$i$-th component is $\partial f(x)/\partial x_i$.
The Hessian $\nabla_{xx} f(x)$ of $f(x)$ is the symmetric matrix
whose $i,j$-th entry is $\partial^2 f(x)/\partial x_i \partial x_j$.
The Hessian is sparse if a significant and useful proportion of the
entries are universally zero.

## Method

An adaptive cubic regularization method is used.
In this, an improvement to a current
estimate of the required minimizer, $x_k$ is sought by computing a
step $s_k$. The step is chosen to approximately minimize a model
$m_k(s)$ of $f(x_k + s)$ that includes a weighted term
$\sigma_k \|s_k\|^3$ for some specified positive weight
$\sigma_k$. The quality of the resulting step $s_k$ is
assessed by computing the "ratio"
$(f(x_k) - f(x_k + s_k))/ (m_k(0) - m_k(s_k))$.
The step is deemed to have succeeded if the ratio exceeds a given
$\eta_s > 0$, and in this case $x_{k+1} = x_k + s_k$.
Otherwise $x_{k+1} = x_k$, and the weight is increased by
powers of a given increase factor up to a given limit. If the ratio is
larger than $\eta_v \geq \eta_d$, the weight will be decreased by
powers of a given decrease factor again up to a given limit. The method
will terminate as soon as $\|\nabla_x f(x_k)\|$ is smaller than
a specified value.

Either linear or quadratic models $m_k(s)$ may be used.
The former will be taken as the first two terms
$f(x_k) + s^T \nabla_x f(x_k)$
of a Taylor series about $x_k$, while the latter uses an
approximation to the first three terms
$f(x_k) + s^T \nabla_x f(x_k) + \frac{1}{2} s^T B_k s$,
for which $B_k$ is a symmetric approximation to the Hessian
$\nabla_{xx}f(x_k)$; possible approximations include the
true Hessian, limited-memory secant and sparsity approximations and
a scaled identity matrix. Normally a two-norm regularization will be used,
but this may change if preconditioning is employed.

An approximate minimizer of the cubic model
is found using either a direct approach involving factorization or an
iterative (conjugate-gradient/Lanczos) approach based on approximations
to the required solution from a so-called Krlov subspace. The direct
approach is based on the knowledge that the required solution
satisfies the linear system of equations
$(B_k + \lambda_k I) s_k = - \nabla_x f(x_k)$
involving a scalar Lagrange multiplier $\lambda_k$.
This multiplier is found by uni-variate root finding, using a safeguarded
Newton-like process, by the GALAHAD packages RQS or DPS
(depending on the norm chosen). The iterative approach
uses the GALAHAD packag GLRT, and is best accelerated by preconditioning
with good approximations to $B_k$ using GALAHAD's PSLS.
The iterative approach has the advantage that only matrix-vector products
$B_k v$ are required, and thus $B_k$ is not required
explicitly. However when factorizations of $B_k$ are possible,
the direct approach is often more efficient.

## References

The generic adaptive cubic regularization method is described in detail in

C. Cartis, N. I. M. Gould and Ph. L. Toint,
“Adaptive cubic regularisation methods for unconstrained optimization.
Part I: motivation, convergence and numerical results”
Mathematical Programming 127(2) (2011) 245-295,

and uses “tricks” as suggested in

N. I. M. Gould, M. Porcelli and Ph. L. Toint,
“Updating the regularization parameter in the adaptive cubic regularization
algorithm”.
Computational Optimization and Applications 53(1) (2012) 1-22.

# Call order

To solve a given problem, functions from the arc package must be called
in the following order:

- arc\_initialize - provide default control parameters and
set up initial data structures
- arc\_read\_specfile (optional) - override control values
by reading replacement values from a file
- arc\_import - set up problem data structures and fixed
values
- arc\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- solve the problem by calling one of
 - arc\_solve\_with\_mat - solve using function calls to
 evaluate function, gradient and Hessian values
 - arc\_solve\_without\_mat - solve using function calls to
 evaluate function and gradient values and Hessian-vector products
 - arc\_solve\_reverse\_with\_mat - solve returning to the
 calling program to obtain function, gradient and Hessian values, or
 - arc\_solve\_reverse\_without\_mat - solve returning to the
 calling prorgram to obtain function and gradient values and
 Hessian-vector products
- arc\_information (optional) - recover information about
the solution and solution process
- arc\_terminate - deallocate data structures

#  Symmetric matrix storage formats

The symmetric $n$ by $n$ matrix $H = \nabla_{xx}f$ may be
presented and stored in a variety of formats. But crucially symmetry is
exploited by only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal).

Both C-style (0 based) and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; the discussion below presumes C style, but add 1 to
indices for the corresponding fortran version.

Wrappers will automatically convert between 0-based (C) and 1-based
(fortran) array indexing, so may be used transparently from C. This
conversion involves both time and memory overheads that may be avoided
by supplying data that is already stored using 1-based indexing.

## Dense storage format

The matrix $H$ is stored as a compact dense matrix by rows, that is,
the values of the entries of each row in turn are
stored in order within an appropriate real one-dimensional array.
Since $H$ is symmetric, only the lower triangular part (that is the part
$H_{ij}$ for $0 \leq j \leq i \leq n-1$) need be held.
In this case the lower triangle should be stored by rows, that is
component $i \ast i / 2 + j$ of the storage array H_val
will hold the value $H_{ij}$ (and, by symmetry, $H_{ji}$)
for $0 \leq j \leq i \leq n-1$.

##  Sparse co-ordinate storage format

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $H$,
its row index i, column index j
and value $H_{ij}$, $0 \leq j \leq i \leq n-1$, are stored as
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
$H_{ij}$ of the entries in the i-th row are stored in components
l = H_ptr(i), $\ldots$, H_ptr(i+1)-1 of the
integer array H_col, and real array H_val, respectively.
Note that as before only the entries in the lower triangle should be stored.
For sparse matrices, this scheme almost always requires less storage than
its predecessor.
