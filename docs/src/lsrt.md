# Introduction

## Purpose

Given a real $m$ by $n$ matrix $A$, a real $m$ vector $b$ and
scalars $\sigma>0$ and $p \geq 2$, this package finds an
**approximate minimizer of the regularised linear-least-squares
objective function
$\frac{1}{2}\| A x - b\|_2^2 + \frac{1}{p} \sigma \| x \|_2^p$.
**
This problem commonly occurs as a subproblem in nonlinear
optimization calculations involving cubic regularisation,
and may be used to regularise the solution
of under-determined or ill-conditioned linear least-squares problems.
The method may be suitable for large $m$ and/or $n$ as no factorization
involving $A$ is required. Reverse communication is used to obtain
matrix-vector products of the form $u + A v$ and
$v + A^T u$.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

November 2007, C interface December 2021.

## Terminology

The required solution $x$ necessarily satisfies the optimality condition
$A^T ( A x - b ) + \lambda x = 0$, where the multiplier
$\lambda = \sigma \|x\|_2^{p-2}$.

## Method

The method is iterative. Startingwith the vector $u_1 = b$, a
bi-diagonalisation process is used to generate the vectors $v_k$ and
$u_k+1$ so that the $n$ by $k$ matrix
$V_k = ( v_1 \ldots v_k)$
and the $m$ by $(k+1)$ matrix $U_k = ( u_1 \ldots u_{k+1})$
together satisfy
$A V_k = U_{k+1} B_k \;\mbox{and}\; b = \|b\|_2 U_{k+1} e_1$
\n
\n
where $B_k$ is $(k+1)$ by $k$ and lower bi-diagonal,
$U_k$ and$V_k$ have orthonormal columns
and $e_1$ is the first unit vector.The solution sought is of the
form $x_k = V_k y_k$, where $y_k$
solves the bi-diagonal regularised least-squares problem
$(1) \;\;\; \min \| B_k y - \|b\| e_1 \|_2 + \frac{1}{p} \sigma \|y \|_2^p.$
\n
 (1) min || B_k y - ||b|| e_1 ||_2+ 1/p sigma || y||^p_2.
\n
To minimize (1), the optimality conditions
$( B_k^T ( B_k^{} y(\lambda) - \|b\| e_1^{} ) + \lambda y(\lambda) = 0,$
\n
\n
where $\lambda = \sigma \|y(\lambda)\|_2^{p-2} \|$,
are used as the basis of an iteration.
The vector $y(\lambda)$ is equivalently the solution to the
regularised least-squares problem
$(2) \;\;\; \min\left \| \vect{ B_k \\ \lambda^{\frac{1}{2}} I } y - \|b\| e_1^{} \right \|_2.$
\n
(2)min||B_k y - ||b|| e_1 ||
||lambda^{1/2} y||
\n
Thus, given an estimate $\lambda \geq 0$, (2) may be efficiently
solved to give $y(\lambda)$.
It is then simply a matter of adjusting $\lambda$
(for example by a Newton-like process) to solve the scalar nonlinear equation
$(3) \;\;\; \theta(\lambda) \equiv\| y(\lambda) \|_2^{p-2} - \frac{\lambda}{\sigma} = 0.$
\n
 (3) theta(lambda) = || y(lambda) ||_2^{p-2} - lambda/sigma = 0.
\n
In practice (3) is reformulated, and a more rapidly converging
iteration is used. Having found$y_k$, a second pass in which
$x_k = V_k y_k$ is regenerated is needed---this need only be done
once $x_k$ has implicitly deemed to be sufficiently close to optimality.
As this second pass is an additional expense, a record is kept of the
optimal objective function values for each value of $k$, and the second
pass is only performed so far as to ensure a given fraction of the
final optimal objective value. Large savings may be made in the second
pass by choosing the required fraction to be significantly smaller than one.

Special code is used in the special case $p=2$, as in this case
a single pass suffices.

## Reference

A complete description of the un- and quadratically-regularised
cases is given by

C. C. Paige and M. A. Saunders,
LSQR: an algorithm for sparse linear equations and sparse leastsquares.
ACM Transactions on Mathematical Software, 8(1):43--71, 1982

and

C. C. Paige and M. A. Saunders,
ALGORITHM 583: LSQR: an algorithm for sparse linear equations and
sparse least squares.
ACM Transactions on Mathematical Software, 8(2):195--209, 1982.

Additional details on the Newton-like process needed to determine
$\lambda$ and other details are described in

C. Cartis, N. I. M. Gould and Ph. L. Toint,
Trust-region and other regularisation of linear
least-squares problems.
BIT 49(1):21-53 (2009).

## Call order

To solve a given problem, functions from the lsrt package must be called
in the following order:

- lsrt\_initialize - provide default control parameters and set up initial data structures
- lsrt\_read\_specfile (optional) - override control values by reading replacement values from a file
- lsrt\_import\_control - import control parameters prior to
solution
- lsrt\_solve_problem - solve the problem by reverse
communication, a sequence of calls are made under control of a status
parameter, each exit either asks the user to provide additional
informaton and to re-enter, or reports that either the solution has
been found or that an error has occurred
- lsrt\_information (optional) - recover information about the solution and solution process
- lsrt\_terminate - deallocate data structures

See Section~\ref{examples} for an example of use.
See the <a href="examples.html">examples tab</a> for an illustration of use.
See the examples section for an illustration of use.