# Introduction

## Purpose

This package uses a **regularization method to find a (local)
unconstrained minimizer of a differentiable weighted sum-of-squares objective
function
$\mathbf{f(x) :=
 \frac{1}{2} \sum_{i=1}^m w_i^{} c_i^2(x) \equiv \frac{1}{2} \|c(x)\|^2_W}$
\n
f(x):= 1/2 sum_{i=1}^m w_i c_i^2(x) = 1/2 ||c(x)||^2_W
\n
of many variables $\mathbf{x}$ involving positive weights
$\mathbf{w_i}$, $\mathbf{i=1,\ldots,m}$.**
The method offers the choice of direct and iterative solution of the key
regularization subproblems, and is most suitable for large problems.
First derivatives of the <i>residual function</i>
$c(x)$ are required, and if second derivatives of the
$c_i(x)$ can be calculated, they may be exploited---if suitable products
of the first or second derivatives with a vector may be found but not the
derivatives themselves, that can also be used to advantage.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

October 2016, C interface August 2021.

## Terminology

The gradient $\nabla_x f(x)$ of a function $f(x)$ is the vector
whose $i$-th component is $\partial f(x)/\partial x_i$.
The Hessian $\nabla_{xx} f(x)$ of $f(x)$ is the symmetric matrix
whose $i,j$-th entry is $\partial^2 f(x)/\partial x_i \partial x_j$.
The Hessian is sparse if a significant and useful proportion of the
entries are universally zero.

The algorithm used by the package is iterative. From the current best estimate
of the minimizer $x_k$, a trial improved point $x_k + s_k$ is sought.
The correction $s_k$ is chosen to improve a model $m_k(s)$ of
%the stabilised objective function $f_{\rho,p}(x_k+s)$ built around
the objective function $f(x_k+s)$ built around
$x_k$. The model is the sum of two basic components,
a suitable approximation $t_k(s)$ of $f(x_k+s)$,
%another approximation of $(\rho/r) \|x_k+s\|_r^r$
(if $\rho > 0$),
and a regularization term $(\sigma_k/p) \|s\|_{S_k}^p$
involving a weight $\sigma_k$, power $p$ and
a norm $\|s\|_{S_k} := \sqrt{s^T S_k s}$ for a given positive
definite scaling matrix $S_k$ that is included to prevent large
corrections. The weight$\sigma_k$ is adjusted as the algorithm
progresses toensure convergence.

The model $t_k(s)$ is a truncated Taylor-series approximation, and this
relies on being able to compute or estimate derivatives of $c(x)$.
Various models are provided, and each has different derivative requirements.
We denote the $m$ by $n$ <i>residual Jacobian</i>
$J(x) \equiv \nabla_x c(x)$ as the matrixwhose $i,j$-th component
$J(x)_{i,j} := \partial c_i(x) / \partial x_j \;\;
\mbox{for $i=1,\ldots,m$ and $j=1,\ldots,n$.}$
\n
J(x)_{i,j} := partial c_i(x) / \partial x_j
\n
\manonly for i=1,...,m and j=1,...,n.For a given $m$-vector $y$, the
<i>weighted-residual Hessian</i> is the sum
$H(x,y) := \sum_{\ell=1}^m y_\ell H_\ell(x),
\;\; \mbox{where}\;\;
H_\ell(x)_{i,j} := \partial^2 c_\ell(x) / \partial x_i \partial x_j
\;\; \mbox{for $i,j=1,\ldots,n$}$
\n
H(x,y) := sum_{\ell=1}^m y_\ell H_\ell(x), where
\n
H_l(x)_{i,j} := partial^2 c_l(x) / partial x_i partial x_j
\n
for i,j=1,...,nis the Hessian of $c_\ell(x)$.
Finally, for a given vector $v$, we define
the <i>residual-Hessians-vector product matrix</i>
$P(x,v) := (H_1(x) v, \ldots, H_m(x) v).$
\n
P(x,v) := (H_1(x) v, ..., H_m(x) v).
\n
The models $t_k(s)$ provided are,
-# the first-order Taylor approximation $f(x_k) + g(x_k)^T s$,
where
$g(x) = J^T(x) W c(x)$,
-# a barely second-order approximation
$f(x_k) + g(x_k)^T s + \frac{1}{2} s^T W s$,
-# the Gauss-Newton approximation
$\frac{1}{2} \| c(x_k) + J(x_k) s\|^2_W$,
-# the Newton (second-order Taylor) approximation
$f(x_k) + g(x_k)^T s + \frac{1}{2} s^T
[ J^T(x_k) W J(x_k) + H(x_k,W c(x_k))] s$, and
-# the tensor Gauss-Newton approximation
$\frac{1}{2} \| c(x_k) + J(x_k) s +
 \frac{1}{2} s^T \cdot P(x_k,s) \|^2_W$,
where the $i$-th component of $s^T \cdot P(x_k,s)$ is
shorthand for the scalar $s^T H_i(x_k) s$,
where $W$ is the diagonal matrix of weights
$w_i$, $i = 1, \ldots m$.

Access to a particular model requires that the user is either able to
provide the derivatives needed (“<i>matrix available</i>”)
or that the products
of these derivatives (and their transposes) with specified vectors are
possible (“<i>matrix free</i>”).

## Method

An adaptive regularization method is used.
In this, an improvement to a current
estimate of the required minimizer, $x_k$ is sought by computing a
step $s_k$. The step is chosen to approximately minimize a model
$t_k(s)$ of $f_{\rho,r}(x_k+s)$
that includes a weighted regularization term
$(\sigma_k/p) \|s\|_{S_k}^p$
for some specified positive weight $\sigma_k$. The quality of the
resulting step $s_k$ is assessed by computing the "ratio"
%$(f_{\rho,p}(x_k) - f_{\rho,p}(x_k+s_k))/(t_k(0)-t_k(s_k))$.
$(f(x_k) - f(x_k + s_k))/(t_k(0) - t_k(s_k))$.
The step is deemed to have succeeded if the ratio exceeds a given
$\eta_s > 0$,
and in this case $x_{k+1} = x_k + s_k$. Otherwise
$x_{k+1} = x_k$, and the weight is increased by powers of a given
increase factor up to a given limit. If the ratio is larger than
$\eta_v \geq \eta_d$, the weight will be decreased by powers of a given
decrease factor again up to a given limit. The method will terminate
as soon as $f(x_k)$ or
$\|\nabla_x f(x_k)\|$ is smaller than a specified value.

A choice of linear, quadratic or quartic models $t_k(s)$ is available
(see the \ref section), and normally a two-norm
regularization willbe used, but this may change if preconditioning
is employed.

If linear or quadratic models are employed, an appropriate,
approximate model minimizer is found using either a direct approach
involving factorization of a shift of the model Hessian $B_k$ or an
iterative (conjugate-gradient/Lanczos) approach based on approximations
to the required solution from a so-called Krlov subspace. The direct
approach is based on the knowledge that the required solution
satisfies the linear system of equations $(B_k + \lambda_k I) s_k
= - \nabla_x f(x_k)$ involving a scalar Lagrange multiplier $\lambda_k$.
This multiplier is found by uni-variate root finding, using a safeguarded
Newton-like process, by the GALAHAD packages RQS. The iterative approach
uses the GALAHAD packag GLRT, and is best accelerated by preconditioning
with good approximations to the Hessian of the model using GALAHAD's PSLS. The
iterative approach has the advantage that only Hessian matrix-vector products
are required, and thus the Hessian $B_k$ is not required explicitly.
However when factorizations of the Hessian are possible, the direct approach
is often more efficient.

When a quartic model is used, the model is itself of least-squares form,
and the package calls itself recursively to approximately minimize its
model. The quartic model often gives a better approximation, but at the
cost of more involved derivative requirements.

## Reference

The generic adaptive cubic regularization method is described in detail in

C. Cartis,N. I. M. Gould and Ph. L. Toint,
“Adaptive cubic regularisation methods for unconstrained optimization.
Part I: motivation, convergence and numerical results”,
Mathematical Programming 127(2) (2011) 245-295,

and uses “tricks” as suggested in

N. I. M. Gould, M. Porcelli and Ph. L. Toint,
“Updating the regularization parameter in the adaptive cubic regularization
algorithm”.
Computational Optimization and Applications 53(1) (2012) 1-22.

The specific methods employed here are discussed in

N. I. M. Gould, J. A. Scott and T. Rees,
“Convergence and evaluation-complexity analysis of a regularized
tensor-Newton method for solving nonlinear least-squares problems”.
Computational Optimization and Applications 73(1) (2019) 1–35.

## Call order

To solve a given problem, functions from the nls package must be called
in the following order:

- nls\_initialize - provide default control parameters and
set up initial data structures
- nls\_read\_specfile (optional) - override control values
by reading replacement values from a file
- nls\_import - set up problem data structures and fixed
values
- nls\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- solve the problem by calling one of
 - nls\_solve\_with\_mat - solve using function calls to
 evaluate function, gradient and Hessian values
 - nls\_solve\_without\_mat - solve using function calls to
 evaluate function and gradient values and Hessian-vector products
 - nls\_solve\_reverse\_with\_mat - solve returning to the
 calling program to obtain function, gradient and Hessian values, or
 - nls\_solve\_reverse\_without\_mat - solve returning to the
 calling prorgram to obtain function and gradient values and
 Hessian-vector products
- nls\_information (optional) - recover information about
the solution and solution process
- nls\_terminate - deallocate data structures

##  Unsymmetric matrix storage formats

The unsymmetric $m$ by $n$ Jacobian matrix
$J \equiv \nabla_x c(x)$ and the residual-Hessians-vector
product matrix $P(x,v)$ may be presented
and stored in a variety of convenient input formats. Let
$A$ be $J$ or $P$ as appropriate.

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

### unsymmetric\_matrix_dense_cols Dense by columns storage format

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
its predecessor.

### unsymmetric\_matrix_column_wise Sparse column-wise storage format

Once again only the nonzero entries are stored, but this time
they are ordered so that those in column j appear directly before those
in column j+1. For the j-th column of $A$ the j-th component of the
integer array A_ptr holds the position of the first entry in this column,
while A_ptr(n) holds the total number of entries plus one.
The row indices i, $0 \leq i \leq m-1$, and values $A_{ij}$
of thenonzero entries in the j-th columnsare stored in components
l = A_ptr(j), $\ldots$, A_ptr(j+1)-1, $0 \leq j \leq n-1$,
of the integer array A_row, and real array A_val, respectively.
As before, for sparse matrices, this scheme almost always requires less
storage than the co-ordinate format.

##  Symmetric matrix storage formats

Likewise, the symmetric $n$ by $n$ weighted-residual
Hessian matrix $H = H(x,y)$ may be presented
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

