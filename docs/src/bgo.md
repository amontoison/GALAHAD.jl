# Introduction

## Purpose

The bgo package uses a **multi-start trust-region method to find an
approximation to the global minimizer of a differentiable objective
function $f(x)$ of $n$ variables $x$, subject to simple bounds
$x^l \leq x \leq x^u$ on the variables.**
Here, any of the components of the vectors of bounds $x^l$ and $x^u$
may be infinite. The method offers the choice of direct
and iterative solution of the key trust-region subproblems, and
is suitable for large problems. First derivatives are required,
and if second derivatives can be calculated, they will be exploited---if
the product of second derivatives with a vector may be found but
not the derivatives themselves, that may also be exploited.

The package offers both random multi-start and local-minimize-and probe
methods to try to locate the global minimizer. There are no theoretical
guarantees unless the sampling is huge, and realistically the success of
the methods decreases as the dimension and nonconvexity increase.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

July 2016, C interface August 2021.

## Terminology

The gradient $\nabla_x f(x)$ of $f(x)$ is the vector whose
$i$-th component is $\partial f(x)/\partial x_i$.
The Hessian $\nabla_{xx} f(x)$ of $f(x)$ is the symmetric matrix
whose $i,j$-th entry is $\partial^2 f(x)/\partial x_i \partial x_j$.
The Hessian is sparse if a significant and useful proportion of the
entries are universally zero.

## Method

A choice of two methods is available.
In the first, local-minimization-and-probe, approach, local minimization
and univariate global minimization are intermixed. Given a current
champion $x^S_k$, a local minimizer $x_k$ of $f(x)$ within the
feasible box $x^l \leq x \leq x^u$ is found using the GALAHAD package trb.
Thereafter $m$ random directions $p$ are generated, and univariate
local minimizer of $f(x_k + \alpha p)$ as a function of the scalar
$\alpha$ along each $p$ within the interval $[\alpha^L,\alpha^u]$,
where $\alpha^L$ and $\alpha^u$ are the smallest and largest
$\alpha$ for which $x^l \leq x_k + \alpha p \leq x^u$,
is performed using the GALAHAD package ugo. The point $x_k + \alpha p$
that gives the smallest value of $f$ is then selected as the new champion
$x^S_{k+1}$.

The random directions $p$ are chosen in one of three ways. The simplest is
to select the components as
(ignore next phrase - doxygen bug!)
$p_i = \mbox{pseudo random} \in
\left\{
\begin{array}{rl}
\mbox{[-1,1]} & \mbox{if} \;\; x^l_i < x_{k,i} < x^u_i \\
\mbox{[0,1]} & \mbox{if} \;\; x_{k,i} = x^l_i \\
\mbox{[-1,0]} & \mbox{if} \;\; x_{k,i} = x^u_i
\end{array}
\right.
$
\n
 ( [-1,1] if x^l_i < x_{k,i} < x^u_i
p_i = pseudo random in ( [0,1] if x_{k,i} = x^l_i
 ( [-1,0] if x_{k,i} = x^u_i
\n
for each $1 \leq i \leq n$. An alternative is to
pick $p$ by partitioning each dimension of the feasible “hypercube” box
into $m$ equal segments, and then selecting sub-boxes
randomly within this hypercube using GALAHAD's Latin hypercube sampling
package, lhs.
Each components of $p$ is then selected in its sub-box, either uniformly
or pseudo randomly.

The other, random-multi-start, method provided selects
$m$ starting points
at random, either componentwise pseudo randomly in the feasible box, or by
 partitioning each component into $m$ equal segments, assigning each to
a sub-box using Latin hypercube sampling, and finally choosing the
values either uniformly or pseudo randomly. Local minimizers within the
feasible box are then computed by the GALAHAD package trb, and
the best is assigned as the current champion. This process is then
repeated until evaluation limits are achieved.

If $n=1$, the GALAHAD package UGO is called directly.

We reiterate that there are no theoretical guarantees unless the sampling
is huge, and realistically the success of the methods decreases as the
dimension and nonconvexity increase. Thus the methods used should best
be viewed as heuristics.

## References

The generic bound-constrained trust-region method is described in detail in

A. R. Conn, N. I. M. Gould and Ph. L. Toint (2000),
Trust-region methods.
SIAM/MPS Series on Optimization,

the univariate global minimization method employed is an extension of that
due to

D. Lera and Ya. D. Sergeyev (2013),
“Acceleration of univariate global optimization algorithms working with
Lipschitz functions and Lipschitz first derivatives”
SIAM J. Optimization Vol. 23, No. 1, pp. 508–529,

while the Latin-hypercube sampling method employed is that of

B. Beachkofski and R. Grandhi (2002).
“Improved Distributed Hypercube Sampling”,
43rd AIAA structures, structural dynamics, and materials conference,
pp. 2002-1274.

# Call order

To solve a given problem, functions from the bgo package must be called
in the following order:

- bgo\_initialize - provide default control parameters and set up initial data structures
- bgo\_read\_specfile (optional) - override control values by reading replacement values from a file
- bgo\_import - set up problem data structures and fixed values
- bgo\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- solve the problem by calling one of
- bgo\_solve\_with\_mat - solve using function calls to evaluate function, gradient and Hessian values
- bgo\_solve\_without\_mat - solve using function calls to evaluate function and gradient values and Hessian-vector products
- bgo\_solve\_reverse\_with\_mat - solve returning to the calling program to obtain function, gradient and Hessian values, or
- bgo\_solve\_reverse\_without\_mat - solve returning to the calling prorgram to obtain function and gradient values and Hessian-vector products
- bgo\_information (optional) - recover information about the solution and solution process
- bgo\_terminate - deallocate data structures

#  Symmetric matrix storage formats

The symmetric $n$ by $n$ matrix $H = \nabla_{xx}f$ may be
presented and stored in a variety of formats. But crucially symmetry
is exploited by only storing values from the lower triangular part
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
