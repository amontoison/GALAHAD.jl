# Introduction

## Purpose

The ugo package aims to find the **global minimizer of a univariate
twice-continuously differentiable function $f(x)$ of a single variable
over the finite interval $x^l \leq x \leq x^u$.** Function and
derivative values may be provided either via a subroutine call,
or by a return to the calling program. Second derivatives may be used
to advantage if they are available.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

July 2016, C interface August 2021.

## Method

The algorithm starts by splitting the interval $[x^l,x^u]$ into a
specified number of subintervals $[x_i,x_{i+1}]$ of equal length,
and evaluating $f$ and its derivatives at each $x_i$. A surrogate
(approximating) lower bound function is constructed on each subinterval
using the function and derivative values at each end, and an estimate of
the first- and second-derivative Lipschitz constant. This surrogate is
minimized, the true objective evaluated at the best predicted point,
and the corresponding interval split again at this point.
Any interval whose surrogate lower bound value exceeds an evaluated
actual value is discarded. The method continues until only one interval
of a maximum permitted width remains.

## References

Many ingredients in the algorithm are based on the paper

D. Lera and Ya. D. Sergeyev (2013),
“Acceleration of univariate global optimization algorithms working with
Lipschitz functions and Lipschitz first derivatives”
SIAM J. Optimization Vol. 23, No. 1, pp. 508–529,

but adapted to use second derivatives.

# Call order

To solve a given problem, functions from the ugo package must be called
in the following order:

- ugo\_initialize - provide default control parameters and
set up initial data structures
- ugo\_read\_specfile (optional) - override control values
by reading replacement values from a file
- ugo\_import - set up problem data structures and fixed
values
- ugo\_reset\_control (optional) - possibly change control
parameters if a sequence of problems are being solved
- solve the problem by calling one of
 - ugo\_solve_direct - solve using function calls to
 evaluate function and derivative values, or
 - ugo\_solve\_reverse - solve returning to the
 calling program to obtain function and derivative values
- ugo\_information (optional) - recover information about
the solution and solution process
- ugo\_terminate - deallocate data structures
