# Introduction

## Purpose

This package **solves sparse unsymmetric system of linear equations.**.
Given an $m$ by $n$ sparse matrix $A = a_{ij}$, the package
solves the system $A x = b$ (or optionally $A^T x = b$).
The matrix $A$ can be rectangular.

**N.B.** The package is simply a sophisticated interface to the
HSL package MA33, and requires that a user has obtained the latter.
** MA33 is not included in GALAHAD**
but is available without charge to recognised academics, see
https://www.hsl.rl.ac.uk/archive/specs/ma33.pdf .
The package offers additional features to MA33.
The storage required for the factorization is chosen
automatically and, if there is insufficient space for the factorization,
more space is allocated and the factorization is repeated.The package
also returns the number of entries in the factors and has facilities
for identifying the rows and columns that are treated specially
when the matrix is singular or rectangular.

Currently, only the control and inform parameters are exposed;
these are provided and used by other GALAHAD packages with C interfaces.
Extended functionality is available using the GALAHAD package uls.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

March 2006, C interface December 2021.

##  Unsymmetric matrix storage formats

The unsymmetric $m$ by $n$matrix $A$ may be presented
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

