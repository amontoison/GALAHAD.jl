# Introduction

## Purpose

Given a **block, symmetric matrix**
$K_H = \mat{cc}{ H & A^T \\ A & - C },$
\n
K_H = ( HA^T )
( A- C )
\n
this package constructs a variety of **preconditioners** of the form
$K_{G} = \mat{cc}{ G & A^T \\ A & - C }.$
\n
K_G = ( GA^T ).
( A- C )
\n
Here, the leading-block matrix $G$ is a suitably-chosen
approximation to $H$; it may either be prescribed **explicitly**, in
which case a symmetric indefinite factorization of $K_G$
will be formed using the GALAHAD symmetric matrix factorization package SLS,
or **implicitly**, by requiring certain sub-blocks of $G$
be zero. In the latter case, a factorization of $K_G$ will be
obtained implicitly (and more efficiently) without recourse to SLS.
In particular, for implicit preconditioners, a reordering
$K_G = P
\mat{ccc}{ G_{11}^{} & G_{21}^T & A_1^T \\G_{21}^{} & G_{22}^{} & A_2^T \\
A_{1}^{} & A_{2}^{} & - C} P^T
$
\n
( G_11G_21^TA_1^T )
K_G = P ( G_21 G_22 A_2^T ) P^T
(A_1 A_2 - C)
\n
involving a suitable permutation $P$ will be found, for some
invertible sub-block (“basis”) $A_1$ of the columns of $A$;
the selection and factorization of $A_1$ uses
the GALAHAD unsymmetric matrix factorization package ULS.
Once the preconditioner has been constructed,
solutions to the preconditioning system
$\mat{cc}{ G & A^T \\ A& - C } \vect{ x \\ y }
 = \vect{a \\ b}
$
\n
( GA^T ) ( x ) = ( a )
( A- C ) ( y ) ( b )
\n
may be obtained by the package.
Full advantage is taken of any zero coefficients in the matrices $H$,
$A$ and $C$.

## Authors

H. S. Dollar and N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

April 2006, C interface November 2021.

## Method

The method used depends on whether an explicit or implicit
factorization is required. In the explicit case, the
package is really little more than a wrapper for the GALAHAD
symmetric, indefinite linear solver SLS in
which the system matrix $K_G$ is assembled from its constituents
$A$, $C$ and whichever $G$ is requested by the user.
Implicit-factorization preconditioners are more involved,
and there is a large variety of different possibilities. The
essential ideas are described in detail in

H. S. Dollar, N. I. M. Gould and A. J. Wathen.
“On implicit-factorization constraint preconditioners”.
In Large Scale Nonlinear Optimization (G. Di Pillo and M. Roma, eds.)
Springer Series on Nonconvex Optimization and Its Applications, Vol. 83,
Springer Verlag (2006) 61--82

and

H. S. Dollar, N. I. M. Gould, W. H. A. Schilders and A. J. Wathen
“On iterative methods and implicit-factorization preconditioners for
regularized saddle-point systems”.
SIAM Journal on Matrix Analysis and Applications, 28(1) (2006) 170--189.

The range-space factorization is based upon the decomposition
$K_{G} = \mat{cc}{ G & 0 \\ A & I}
\mat{cc}{ G^{-1} & 0 \\ 0 & - S}\mat{cc}{ G & A^T \\ 0 & I},
$
\n
K_G = ( G0 ) ( G^{-1} 0 ) ( G A^T )
( AI ) ( 0 -S ) ( 0I)
\n
where the “Schur complement” $S = C + A G^{-1} A^T$.
Such a method requires that $S$ is easily invertible. This is often the
case when $G$ is a diagonal matrix, in which case $S$ is frequently
sparse, or when $m \ll n$ in which case $S$
is small and a dense Cholesky factorization may be used.

When $C = 0$, the null-space factorization is based upon the decomposition
$K_{G} = P\mat{ccc}{G_{11}^{} & 0 & I \\
G_{21}^{} & I & A_{2}^{T} A_{1}^{-T} \\A_{1}^{} & 0 & 0 }
\mat{ccc}{0 & 0 & I \\ \;\;\; 0 \;\; & \;\; R \;\; & 0 \\ I & 0 & - G_{11}^{}}
\mat{ccc}{G_{11}^{} & G_{21}^T & A_{1}^T \\0 & I & 0 \\
I & A_{1}^{-1} A_{2}^{} & 0} P^T,
$
\n
( G_110I) ( 00 I )
K_G = P ( G_21IA_2^T A_1^{-T} ) ( 0R 0 )
( A_1 00) ( I0 -G_11 )

( G_11 G_21^T A_1^T )
. (0I0) P^T,
(IA_1^{-1} A_20 )
\n
where the “reduced Hessian”
$R = ( - A_{2}^{T} A_1^{-T} \;\; I )
\mat{cc}{G_{11}^{} & G_{21}^{T} \\ G_{21}^{} & G_{22}^{}}
\vect{ - A_1^{-1} A_2^{} \\ I}
$
\n
 R = ( -A_2^T A_1^{-T}I )( G_11G_21^T ) ( -A_1^{-1} A_2 )
 ( G_21 G_22) ( I )
\n
and $P$ is a suitably-chosen permutation for which $A_1$ is
invertible. The method is most useful when $m \approx n$ as then the
dimension of $R$ is small and a dense Cholesky factorization may be used.

## Call order

To solve a given problem, functions from the sbls package must be called
in the following order:

- sbls\_initialize - provide default control parameters and set up initial data structures
- sbls\_read\_specfile (optional) - override control values by reading replacement values from a file
- sbls\_import - set up matrix data structures
- sbls\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- sbls_factorize\_matrix - form and factorize the block
matrix from its components
- sbls\_solve_system - solve the block linear system of
equations
- sbls\_information (optional) - recover information about the solution and solution process
- sbls\_terminate - deallocate data structures

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

Likewise, the symmetric $n$ by $n$ matrix $H$, as well as
the $m$ by $m$ matrix $C$,may be presented
and stored in a variety of formats. But crucially symmetry is exploited
by only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal). We focus
on $H$, but everything we say applies equally to $C$.

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

