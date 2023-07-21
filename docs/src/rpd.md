# Introduction

## Purpose

 Read and write data for the linear program (LP)
$\mbox{minimize}\;\; g^T x + f
\;\mbox{subject to}\; c_l \leq A x \leq c_u
\;\mbox{and}\; x_l \leqx\leq x_u,
$
\n
minimize g^T x + f
 subject toc_l <= A x <= c_u
 x_l <=x<= x_u,
\n
 the linear program with quadratic constraints (QCP)
$\mbox{minimize}\;\; g^T x + f
\;\mbox{subject to}\; c_l \leq A x + \frac{1}{2} \mbox{vec}(x.H_c.x) \leq c_u
\;\mbox{and}\; x_l \leqx\leq x_u,
$
\n
minimize g^T x + f
 subject toc_l <= A x + 1/2 vec(x.H_c.x) <= c_u
 x_l <=x<= x_u,
\n
 the bound-constrained quadratic program (BQP)
$\mbox{minimize}\;\; \frac{1}{2} x^T H x + g^T x + f
\;\mbox{subject to}\; x_l \leqx\leq x_u,
$
\n
 minimize 1/2 x^T H x + g^T x + f
 subject to x_l <=x<= x_u,
\n
 the quadratic program (QP)
$\mbox{minimize}\;\; \frac{1}{2} x^T H x + g^T x + f
\;\mbox{subject to}\; c_l \leq A x \leq c_u
\;\mbox{and}\; x_l \leqx\leq x_u,
$
\n
 minimize1/2 x^T H x + g^T x + f
 subject toc_l <= A x <= c_u
 x_l <=x<= x_u,
\n
 or the quadratic program with quadratic constraints (QCQP)
$\mbox{minimize}\;\; \frac{1}{2} x^T H x + g^T x + f
\;\mbox{subject to}\; c_l \leq A x + \frac{1}{2} \mbox{vec}(x.H_c.x) \leq c_u
\;\mbox{and}\; x_l \leqx\leq x_u,
$
\n
minimize 1/2 x^T H x + g^T x + f
 subject toc_l <= A x + 1/2 vec(x.H_c.x) <= c_u
 x_l <=x<= x_u,
\n
 where vec$( x . H_c . x )$ is the vector whose
 $i$-th component is$x^T (H_c)_i x$ for the $i$-th
 constraint, from and to a QPLIB-format data file.
 Variables may be continuous, binary or integer.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

January 2006, C interface January 2022.

## Reference

The QPBLIB format is defined in

F. Furini, E. Traversi, P. Belotti, A. Frangioni, A. Gleixner, N. Gould,
L. Liberti, A. Lodi, R. Misener, H. Mittelmann, N. V. Sahinidis,
S. Vigerske and A. Wiegele(2019).
QPLIB: a library of quadratic programming instances,
Mathematical Programming Computation **11** 237–265.

## Call order

To decode a given QPLIB file, functions from the rpd package must be called
in the following order:

- rpd\_initialize - provide default control parameters and set up initial data structures
- rpd_get_stats - read a given QPLIB file into internal
 data structures, and report vital statistics
- (optionally, and in any order, where relevant)
- rpd_get_g - get the objective gradient term $g$
- rpd_get_f - get the objective constant term $f$
- rpd_get_xlu - get the variable bounds
 $x_l$ and $x_u$
- rpd_get_xlu - get the constraint bounds
$c_l$ and $c_u$
- rpd_get_h - get the objective Hessian term $H$
- rpd_get_a - get the constrain Jacobian term $A$
- rpd_get_h_c - get the constraint Hessian terms $H_c$
- rpd_get_x_type - determine the type of each variable
 $x$
- rpd_get_x - get initial value of the variable $x$
- rpd_get_y - get initial value of Lagrange multipliers
$y$
- rpd_get_z - get initial value of the dual variables
$z$
- rpd\_terminate - deallocate data structures

##  Sparse unsymmetric co-ordinate storage format

The unsymmetric $m$ by $n$ constraint matrix $A$ will be
output in sparse co-ordinate format.

Both C-style (0 based)and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; the discussion below presumes C style, but add 1 to
indices for the corresponding fortran version.

Wrappers will automatically convert between 0-based (C) and 1-based
(fortran) array indexing, so may be used transparently from C. This
conversion involves both time and memory overheads that may be avoided
by supplying data that is already stored using 1-based indexing.

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $A$,
its row index i, column index j
and value $A_{ij}$,
$0 \leq i \leq m-1$,$0 \leq j \leq n-1$,are stored as
the $l$-th components of the integer arrays A_row and
A_col and real array A_val, respectively, while the number of nonzeros
is recorded as A_ne = $ne$.

##  Sparse symmetric co-ordinate storage format

Likewise, the symmetric $n$ by $n$ objective Hessian matrix
$H$ will be returned in a sparse co-ordinate format. But crucially
symmetry is exploited by only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal).

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $H$,
its row index i, column index j
and value $h_{ij}$, $0 \leq j \leq i \leq n-1$,are stored as
the $l$-th components of the integer arrays H_row and
H_col and real array H_val, respectively, while the number of nonzeros
is recorded as H_ne = $ne$.
Note that only the entries in the lower triangle should be stored.

## joint_ Joint sparse symmetric co-ordinate storage format

The symmetric $n$ by $n$ constraint Hessian matrices
$ (H_c)_i$ are stored as a whole in a joint symmetric co-ordinate
storage format.
In addition to the row and column indices and values of each lower
triangular matrix, record is also kept of the particular constraint invlved.

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $H$,
its constraint index k, row index i, column index j
and value $(h_k)_{ij}$, $0 \leq j \leq i \leq n-1$,are stored as
the $l$-th components of the integer arrays H_c_ptr, H_c_row and
H_c_col and real array H_c_val, respectively, while the number of nonzeros
is recorded as H_c_ne = $ne$.
Note as before that only the entries in the lower triangles should be stored.

