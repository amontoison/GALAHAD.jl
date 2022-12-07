# Introduction

## Purpose

Given a real $m$ by $n$ matrix $A$, a real
$m$ vector $b$ and a scalar $\Delta>0$, this package finds an
** approximate minimizer of $\| A x - b\|_2$, where the vector
$x$ is required to satisfy the “trust-region”
constraint $\|x\|_2 \leq\Delta$.**
This problem commonly occurs as a trust-region subproblem in nonlinear
optimization calculations, and may be used to regularize the solution
of under-determined or ill-conditioned linear least-squares problems.
The method may be suitable for large $m$ and/or $n$ as no
factorization involving $A$ is required. Reverse communication is used
to obtain matrix-vector products of the form $u + A v$ and
$v + A^T u$.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

November 2007, C interface December 2021.

## Terminology

The required solution $x$ necessarily satisfies the optimality condition
$A^T ( A x - b ) + \lambda x = 0$, where $\lambda \geq 0$
is a Lagrange multiplier corresponding to the trust-region constraint
$\|x\|_2\leq\Delta$.

## Method

The method is iterative. Startingwith the vector $u_1 = b$, a
bi-diagonalisation process is used to generate the vectors $v_k$ and
$u_k+1$ so that the $n$ by $k$ matrix $V_k = ( v_1 \ldots v_k)$
and the $m$ by $(k+1)$ matrix $U_k = ( u_1 \ldots u_{k+1})$
together satisfy
$A V_k = U_{k+1} B_k \;\mbox{and}\; b = \|b\| U_{k+1} e_1,$
\n
 A V_k = U_{k+1} B_k and b = ||b|| U_{k+1} e_1,
\n
where $B_k$ is $(k+1)$ by $k$ and lower bi-diagonal, $U_k$ and
$V_k$ have orthonormal columns and $e_1$ is the first unit vector.
The solution sought is of the form $x_k = V_k y_k$, where $y_k$
solves the bi-diagonal least-squares trust-region problem
$(1) \;\;\; \min \| B_k y - \|b\| e_1 \|_2 \;\mbox{subject to}\; \|y\|_2 \leq \Delta.$
\n
 (1)min || B_k y - \|b\| e_1 ||_2 subject to ||y||_2 <= Delta.
\n

If the trust-region constraint is inactive, the solution $y_k$
may be found, albeit indirectly, via the LSQR algorithm of Paige and Saunders
which solves the bi-diagonal least-squares problem
$ \min \| B_k y - \|b\| e_1 \|_2$
\n
 min || B_k y - ||b|| e_1 ||_2
\n
using a QR factorization of $B_k$. Only the most recent
$v_k$ and $u_{k+1}$
are required, and their predecessors discarded, to compute $x_k$ from
$x_{k-1}$. This method has the important property that the iterates
$y$ (and thus $x_k$) generated increase in norm with $k$. Thus as
soon as an LSQR iterate lies outside the trust-region, the required solution
to (1) and thus to the original problem must lie on the boundary of the
trust-region.

If the solution is so constrained, the simplest strategy is to interpolate
the last interior iterate with the newly discovered exterior one to find the
boundary point---the so-called Steihaug-Toint point---between them.
Once the solution is known to lie on the trust-region boundary,
further improvement may be made by solving
$ \min \| B_k y - \|b\| e_1 \|_2 \;\mbox{subject to}\;|\|y\|_2 = \Delta,$
\n
 min || B_k y - ||b|| e_1 ||_2 subject to ||y||_2 = Delta,
\n
for which the optimality conditions require that $y_k = y(\lambda_k)$
where $\lambda_k$ is the positive root of
$B_k^T ( B_k^{} y(\lambda) - \|b\| e_1^{} ) + \lambday(\lambda) = 0 \;\mbox{and}\;\|y(\lambda)\|_2 = \Delta$
\n
B_k^T ( B_k y(lambda) - ||b|| e_1 ) + lambda y(lambda) = 0
and ||y(lambda)||_2 = Delta
\n
The vector $y(\lambda)$ is equivalently the solution to the
regularized least-squares problem
$\min\left \| \vect{ B_k \\ \lambda^{\frac{1}{2}} I } y - \|b\| e_1^{} \right \|$
\n
min||B_k y - ||b|| e_1 ||
 ||lambda^{1/2} y||
\n
and may be found efficiently. Given$y(\lambda)$, Newton's method
is then used to find $\lambda_k$ as the positive root of
$\|y(\lambda)\|_2 = \Delta$. Unfortunately, unlike when the solution
lies in the interior of the trust-region, it is not known how to recur
$x_k$ from $x_{k-1}$ given $y_k$, and a second pass in which
$x_k = V_k y_k$ is regenerated is needed---this need only be done
once $x_k$ has implicitly deemed to be sufficiently close to optimality.
As this second pass is an additional expense, a record is kept of the
optimal objective function values for each value of $k$, and the second
pass is only performed so far as to ensure a given fraction of the
final optimal objective value. Large savings may be made in the second
pass by choosing the required fraction to be significantly smaller than one.

## Reference

A complete description of the unconstrained case is given by

C. C. Paige and M. A. Saunders,
LSQR: an algorithm for sparse linear equations and sparse leastsquares.
ACM Transactions on Mathematical Software, 8(1):43--71, 1982

and

C. C. Paige and M. A. Saunders,
ALGORITHM 583: LSQR: an algorithm for sparse linear equations and
sparse least squares.
ACM Transactions on Mathematical Software, 8(2):195--209, 1982.

Additional details on how to proceed once the trust-region constraint is
encountered are described in detail in

C. Cartis, N. I. M. Gould and Ph. L. Toint,
Trust-region and other regularisation of linear
least-squares problems.
BIT 49(1):21-53 (2009).

## Call order

To solve a given problem, functions from the lstr package must be called
in the following order:

- lstr\_initialize - provide default control parameters and set up initial data structures
- lstr\_read\_specfile (optional) - override control values by reading replacement values from a file
- lstr\_import\_control - import control parameters prior to
solution
- lstr\_solve_problem - solve the problem by reverse
communication, a sequence of calls are made under control of a status
parameter, each exit either asks the user to provide additional
informaton and to re-enter, or reports that either the solution has
been found or that an error has occurred
- lstr\_information (optional) - recover information about the solution and solution process
- lstr\_terminate - deallocate data structures

See Section~\ref{examples} for an example of use.
See the <a href="examples.html">examples tab</a> for an illustration of use.
See the examples section for an illustration of use.

