# Introduction

## Purpose

Given an under-determined set of linear equations/constraints
$a_i^T x = b_i^{}$, $i = 1, \ldots, m$ involving
$n \geq m$ unknowns $x$, this package
**determines whether the constraints are consistent, and if
so how many of the constraints are dependent**; a list of dependent
constraints, that is, those which may be removed without changing the
solution set, will be found and the remaining $a_i$ will be linearly
independent.Full advantage is taken of any zero coefficients in the
vectors $a_i$.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

August 2006, C interface January 2021

## Method

A choice of two methods is available. In the first, the matrix
$K = \mat{cc}{ \alpha I & A^T \\ A & 0 }$
is formed and factorized for some small $\alpha > 0$ using the
GALAHAD package SLS---the
factors $K = P L D L^T P^T$ are used to determine
whether $A$ has dependent rows. In particular, in exact arithmetic
dependencies in $A$ will correspond to zero pivots in the block
diagonal matrix $D$.

The second choice of method finds factors
$A = P L U Q$ of the rectangular matrix $A$
using the GALAHAD package ULS.
In this case, dependencies in $A$ will be reflected in zero diagonal
entries in $U$ in exact arithmetic.

The factorization in either case may also be used to
determine whether the system is consistent.

## Call order

To solve a given problem, functions from the fdc package must be called
in the following order:

- fdc\_initialize - provide default control parameters and set up initial data structures
- fdc\_read\_specfile (optional) - override control values by reading replacement values from a file
- fdc_find_dependent_rows - find the number of dependent
rows and, if there are any, whether the constraints are
independent
- fdc\_terminate - deallocate data structures

## fdc_array_indexing Array indexing

Both C-style (0 based)and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; add 1 to input integer arrays if fortran-style indexing is
used, and beware that return integer arrays will adhere to this.

