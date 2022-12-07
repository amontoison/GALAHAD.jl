# Introduction

## Purpose

Compute the the **solution to an extended system of $n + m$
sparse real linear equations in $n + m$ unknowns,**
$\mbox{(1)}\;\; \mat{cc}{ A & B \\ C & D } \vect{x_1 \\ x_2} =\vect{b_1 \\ b_2}$
 \n
 (1)( AB ) ( x_1 ) = ( b_1 )
( CD ) ( x_2 ) ( b_2 )
 \n
in the case where the $n$ by $n$ matrix $A$ is nonsingular
and solutions to the systems
$A x=b \;\mbox{and}\; A^T y=c$
 \n
 A x=bandA^T y=c
 \n
may be obtained from an external source, such as an existing
factorization.The subroutine uses reverse communication to obtain
the solution to such smaller systems.The method makes use of
the Schur complement matrix
$S = D - C A^{-1} B.$
 \n
 S = D - C A^{-1} B.$
 \n
The Schur complement is stored and factorized as a dense matrix
and the subroutine is thus appropriate only if there is
sufficient storage for this matrix. Special advantage is taken
of symmetry and definiteness in the coefficient matrices.
Provision is made for introducing additional rows and columns
to, and removing existing rows and columns from, the extended
matrix.

Currently, only the control and inform parameters are exposed;
these are provided and used by other GALAHAD packages with C interfaces.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

March 2005, C interface January 2022.

## Method

The subroutine galahad_factorize forms the Schur complement
$S = D - C A^{-1} B$ of $ A$
in the extended matrix by repeated reverse communication to
obtain the columns of
$A^{-1} B$. The Schur complement or its negative is then factorized
into its QR or, if possible, Cholesky factors.

The subroutine galahad\_solve solves the extended system using
the following well-known scheme:
 -# Compute the solution to $ A u=b_1$;
 -# Compute $x_2$ from $ S x_2=b_2-C u$;
 -# Compute the solution to $ A v=B x_2$; and
 -# Compute $x_1 = u - v$.

The subroutines galahad_append and galahad_delete compute the factorization
of the Schur complement after a row and column have been appended
to, and removed from, the extended matrix, respectively.
The existing factorization is updated
to obtain the new one; this is normally more efficient than
forming the factorization from scratch.

## Call order

To solve a given problem, functions from the scu package must be called
in the following order:

- scu\_initialize - provide default control parameters and
set up initial data structures
- scu\_read\_specfile (optional) - override control values
by reading replacement values from a file
- scu_form_and_factorize - form and factorize the
 Schur-complement matrix $S$
- scu\_solve_system - solve the block system (1)
- scu_add_rows_and_cols (optional) - update the factors of
 the Schur-complement matrix when rows and columns are added to (1).
- scu_delete_rows_and_cols (optional) - update the factors of
 the Schur-complement matrix when rows and columns are removed from (1).
- scu\_information (optional) - recover information about
the solution and solution process
- scu\_terminate - deallocate data structures

