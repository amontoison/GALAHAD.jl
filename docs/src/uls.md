# Introduction

## Purpose

This package
** solves dense or sparse unsymmetric systems of linear equations**
using variants of Gaussian elimination.
Given a sparse symmetric $m \times n$ matrix $A = a_{ij}$, and an
$m$-vector $b$, this subroutine solves the system $A x = b$. If
$b$ is an $n$-vector, the package may solve instead the system
$A^T x = b$. Both square ($m=n$) and
rectangular ($m \neq n$)matrices are handled; one of an infinite
class of solutions for consistent systems will be returned
whenever $A$ is not of full rank.

The method provides a common interface to a variety of well-known solvers
from HSL. Currently supported solvers include MA28/GLS and HSL\_MA48.
Note that ** the solvers themselves do not form part of this package
and must be obtained separately.**
Dummy instances are provided for solvers that are unavailable.
Also note that additional flexibility may be obtained by calling the
solvers directly rather that via this package.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

August 2009,C interface December 2021.

## Terminology

The solvers used each produce an $P_R L U P_C$ factorization
of $A$, where $L$ and $U$ are lower and upper triangular
matrices, and $P_R$ and $P_C$ are row and column permutation
matrices respectively.

## Method

Variants of sparse Gaussian elimination are used.

The solver GLS is available as part of GALAHAD and relies on
the HSL Archive packages MA33. To obtain HSL Archive packages, see

http://hsl.rl.ac.uk/archive/ .

The solver HSL\_MA48 is part of HSL 2007. To obtain HSL 2007 packages, see

http://hsl.rl.ac.uk/hsl2007/ .

## Reference

The methods used are described in the user-documentation for

HSL 2007, A collection of {F}ortran codes for large-scale scientific
computation (2007).\n
http://www.cse.clrc.ac.uk/nag/hsl

## Call order

To solve a given problem, functions from the uls package must be called
in the following order:

- uls\_initialize - provide default control parameters and set up initial data structures
- uls\_read\_specfile (optional) - override control values by reading replacement values from a file
- uls_factorize\_matrix - set up matrix data structures,
 analyse the structure to choose a suitable order for factorization,
 and then factorize the matrix $A$
- uls\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- uls\_solve_system - solve the linear system of
equations $Ax=b$ or $A^Tx=b$
- uls\_information (optional) - recover information about the solution and solution process
- uls\_terminate - deallocate data structures

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

