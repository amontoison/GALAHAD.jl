# Introduction

## Purpose

This package
** solves dense or sparse symmetric systems of linear equations**
using variants of Gaussian elimination.Given a sparse symmetric
$n \times n$ matrix $A$, and an $n$-vector $b$, this
subroutine solves the system $A x = b$.The matrix $A$ need not
be definite.

The package provides a common interface to a variety of well-known
solvers from HSL and elsewhere. Currently supported solvers include
MA27/SILS, HSL\_MA57, HSL\_MA77, HSL\_MA86,
HSL\_MA87 and HSL\_MA97 from HSL,
SSIDS from SPRAL,
PARDISO both from the Pardiso Project and Intel's MKL
and WSMP from the IBM alpha Works, as
well as POTR, SYTR and SBTR from LAPACK.
Note that
** the solvers themselves do not form part of this package and
must be obtained separately.**
Dummy instances are provided for solvers that are unavailable.
Also note that additional flexibility may be obtained by calling the
solvers directly rather that via this package.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

August 2009, C interface December 2021.

## Terminology

The solvers used each produce an $L D L^T$ factorization of
$A$ or a perturbation thereof, where $L$ is a permuted
lower triangular matrix and $D$ is a block diagonal matrix with
blocks of order 1 and 2. It is convenient to write this factorization in
the form
$A + E = P L D L^T P^T,$ where
$P$ is a permutation matrix and $E$ is any diagonal
perturbation introduced.

## sls\_solvers Supported external solvers

The key features of the external solvers supported by sls are
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

## Method

Variants of sparse Gaussian elimination are used.

The solver SILS is available as part of GALAHAD and relies on
the HSL Archive package MA27. To obtain HSL Archive packages, see

http://hsl.rl.ac.uk/archive/ .

The solvers
HSL\_MA57,
HSL\_MA77,
HSL\_MA86,
HSL\_MA87
and
HSL\_MA97, the ordering packages
MC61 and HSL\_MC68, and the scaling packages
HSL\_MC64 and MC77
are all part of HSL 2011.
To obtain HSL 2011 packages, see

http://hsl.rl.ac.uk

The solver SSIDS is from the SPRAL sparse-matrix collection,
and is available as part of GALAHAD.

The solver PARDISO is available from the Pardiso Project;
version 4.0.0 or above is required.
To obtain PARDISO, see

http://www.pardiso-project.org/ .

The solver MKL PARDISO is available as part of Intel's oneAPI Math Kernel
Library (oneMKL).
To obtain this version of PARDISO, see

https://software.intel.com/content/www/us/en/develop/tools/oneapi.html .

The solver WSMP is available from the IBM alpha Works;
version 10.9 or above is required.
To obtain WSMP, see

http://www.alphaworks.ibm.com/tech/wsmp .

The solvers POTR, SYTR and PBTR,
are available as
S/DPOTRF/S,
S/DSYTRF/S and S/DPBTRF/S
as part of LAPACK. Reference versions
are provided by GALAHAD, but for good performance
machined-tuned versions should be used.

Explicit sparsity re-orderings are obtained by calling the HSL package
HSL\_MC68.
Both this, HSL\_MA57 and PARDISO rely optionally
on the ordering package MeTiS (version 4) from the Karypis Lab.
To obtain METIS, see

http://glaros.dtc.umn.edu/gkhome/views/metis/ .

Bandwidth, Profile and wavefront reduction is supported by
calling HSL's MC61.

## Reference

The methods used are described in the user-documentation for

HSL 2011, A collection of Fortran codes for large-scale scientific
 computation (2011). http://www.hsl.rl.ac.uk

and papers

O. Schenk and K. G&auml;rtner,
“Solving Unsymmetric Sparse Systems of Linear Equations with PARDISO”.
Journal of Future Generation Computer Systems \b, 20(3) (2004) 475--487,

O. Schenk and K. G&auml;rtner,
“On fast factorization pivoting methods for symmetric indefinite systems”.
Electronic Transactions on Numerical Analysis \b 23 (2006) 158--179,
and

A. Gupta,
“WSMP: Watson Sparse Matrix Package Part I - direct
solution of symmetric sparse systems”.
IBM Research Report RC 21886, IBM T. J. Watson Research Center,
NY 10598, USA (2010).


## Call order

To solve a given problem, functions from the sls package must be called
in the following order:

- sls\_initialize - provide default control parameters and set up initial data structures
- sls\_read\_specfile (optional) - override control values by reading replacement values from a file
- sls_analyse\_matrix - set up matrix data structures
 and analyse the structure to choose a suitable order for factorization
- sls\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- sls_factorize\_matrix - form and factorize the
matrix $A$
- one of
- sls\_solve_system - solve the linear system of
equations $Ax=b$
- sls_partial\_solve_system - solve a linear system
$Mx=b$ involving one of the matrix factors $M$ of $A$
- sls\_information (optional) - recover information about the solution and solution process
- sls\_terminate - deallocate data structures

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
