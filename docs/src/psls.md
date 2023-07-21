# Introduction

## Purpose

Given an $n$ by $n$sparse symmetric matrix $A =a_{ij}$,
this package **builds a suitable symmetric, positive definite (or
diagonally dominant)-preconditioner $P$ of $A$ or a symmetric
sub-matrix thereof**. The matrix $A$ need not be definite. Facilities
are provided to apply the preconditioner to a given vector, and to
remove rows and columns (symmetrically) from the initial
preconditioner without a full re-factorization.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

April 2008, C interface January 2022.

## Method and references

The basic preconditioners are described in detail in Section 3.3.10 of

A. R. Conn, N. I. M. Gould and Ph. L. Toint (1992).
LANCELOT. A fortran package for large-scale nonlinear optimization
(release A). Springer Verlag Series in Computational Mathematics 17,
Berlin,

along with the more modern versions implements in ICFS due to

C.-J. Lin and J. J. More' (1999).
Incomplete Cholesky factorizations with limited memory.
SIAM Journal on Scientific Computing **21** 21-45,

and in HSL_MI28 described by

J. A. Scott and M. Tuma (2013). HSL MI28: an efficient and robust
limited-memory incomplete Cholesky factorization code.
ACM Transactions on Mathematical Software **40(4)** (2014), Article 24.

The factorization methods used by the GALAHAD package SLS in conjunction
with some preconditioners are described in the documentation to that
package. The key features of the external solvers supported by SLS are
given in the following table.

(ignore next paragraph - doxygen bug!)

<table>
<caption>External solver characteristics</caption>
<tr><th> solver <th> factorization <th> indefinite $A$
<th> out-of-core <th> parallelised
<tr><td> SILS/MA27 <td> multifrontal <td> yes <td> no <td> no
<tr><td> HSL_MA57 <td> multifrontal <td> yes <td> no <td> no
<tr><td> HSL_MA77 <td> multifrontal <td> yes <td> yes <td> OpenMP core
<tr><td> HSL_MA86 <td> left-looking <td> yes <td> no <td> OpenMP fully
<tr><td> HSL_MA87 <td> left-looking <td> no <td> no <td> OpenMP fully
<tr><td> HSL_MA97 <td> multifrontal <td> yes <td> no <td> OpenMP core
<tr><td> SSIDS <td> multifrontal <td> yes <td> no <td> CUDA core
<tr><td> PARDISO <td> left-right-looking <td> yes <td> no <td> OpenMP fully
<tr><td> MKL_PARDISO <td> left-right-looking <td> yes <td> optionally
 <td> OpenMP fully
<tr><td> WSMP <td> left-right-looking <td> yes <td> no <td> OpenMP fully
<tr><td> POTR <td> dense <td> no <td> no <td> with parallel LAPACK
<tr><td> SYTR <td> dense <td> yes <td> no <td> with parallel LAPACK
<tr><td> PBTR <td> dense band <td> no <td> no <td> with parallel LAPACK
</table>

External solver characteristics (ooc = out-of-core factorization)

 solver factorization indefinite Aoocparallelised
 SILS/MA27 multifrontalyes nono
 HSL_MA57multifrontalyes nono
 HSL_MA77multifrontalyesyesOpenMP core
 HSL_MA86left-lookingyes noOpenMP fully
 HSL_MA87left-looking no noOpenMP fully
 HSL_MA97multifrontalyes noOpenMP core
 SSIDS multifrontalyes noCUDA core
 PARDISO left-right-lookingyes noOpenMP fully
 MKL_PARDISO left-right-lookingyesoptionallyOpenMP fully
 WSMPleft-right-lookingyes noOpenMP fully
 POTRdenseno nowith parallel LAPACK
 SYTRdense yes nowith parallel LAPACK
 PBTRdense band no nowith parallel LAPACK

Note that ** the solvers themselves do not form part of this package and
must be obtained separately.**
Dummy instances are provided for solvers that are unavailable.

Orderings to reduce the bandwidth, as implemented in HSL's MC61, are due to

J. K. Reid and J. A. Scott (1999)
Ordering symmetric sparse matrices for small profile and wavefront
International Journal for Numerical Methods in Engineering
**45** 1737-1755.

If a subset of the rows and columns are specified, the remaining rows/columns
are removed before processing. Any subsequent removal of rows and columns
is achieved using the GALAHAD Schur-complement updating package SCU
unless a complete re-factorization is likely more efficient.

## Call order

To solve a given problem, functions from the psls package must be called
in the following order:

- psls\_initialize - provide default control parameters and set up initial data structures
- psls\_read\_specfile (optional) - override control values by reading replacement values from a file
- psls\_import - set up matrix data structures for $A$
prior to solution
- psls\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- one of
- psls_form_preconditioner - form and factorize a
preconditioner $P$ of the matrix $A$
- psls_form_subset_preconditioner - form and factorize a
preconditioner $P$ of a symmetric submatrix of the matrix $A$
- psls_update_preconditioner (optional) - update the
 preconditioner $P$ when rows (amd columns) are removed
- psls_apply_preconditioner - solve the linear system of
equations $Px=b$
- psls\_information (optional) - recover information about
the preconditioner and solution process
- psls\_terminate - deallocate data structures

 ##  Symmetric matrix storage formats

The symmetric $n$ by $n$ coefficient matrix $A$ may be presented
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

The matrix $A$ is stored as a compactdense matrix by rows, that is,
the values of the entries of each row in turn are
stored in order within an appropriate real one-dimensional array.
Since $A$ is symmetric, only the lower triangular part (that is the part
$A_{ij}$ for $0 \leq j \leq i \leq n-1$) need be held.
In this case the lower triangle should be stored by rows, that is
component $i \ast i / 2 + j$of the storage array val
will hold the value $A_{ij}$ (and, by symmetry, $A_{ji}$)
for $0 \leq j \leq i \leq n-1$.

###  Sparse co-ordinate storage format

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $A$,
its row index i, column index j
and value $A_{ij}$, $0 \leq j \leq i \leq n-1$,are stored as
the $l$-th components of the integer arrays row and
col and real array val, respectively, while the number of nonzeros
is recorded as ne = $ne$.
Note that only the entries in the lower triangle should be stored.

###  Sparse row-wise storage format

Again only the nonzero entries are stored, but this time
they are ordered so that those in row i appear directly before those
in row i+1. For the i-th row of $A$ the i-th component of the
integer array ptr holds the position of the first entry in this row,
while ptr(n) holds the total number of entries plus one.
The column indices j, $0 \leq j \leq i$, and values
$A_{ij}$ of theentries in the i-th row are stored in components
l = ptr(i), $\ldots$, ptr(i+1)-1 of the
integer array col, and real array val, respectively.
Note that as before only the entries in the lower triangle should be stored.
For sparse matrices, this scheme almost always requires less storage than
its predecessor.
