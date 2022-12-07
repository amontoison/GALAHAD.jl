# Introduction

## Purpose

This package **solves sparse symmetric system of linear equations.**.
Given an $n$ by $n$ sparse matrix $A = {a_{ij}}$, and an
$n$ vector $b$, the package solves the system $A x = b$.
The matrix $A$ need not be definite. There is an option for iterative
refinement.

Currently, only the control and inform parameters are exposed;
these are provided and used by other GALAHAD packages with C interfaces.
Extended functionality is available using the GALAHAD package sls.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

April 2001, C interface December 2021.

## Method

The method used is a direct method based on a sparse variant
of Gaussian elimination and is discussed further by

I. S. Duff and J. K. Reid (1983),
ACM Trans. Math. Software **9** pp.302-325.

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
