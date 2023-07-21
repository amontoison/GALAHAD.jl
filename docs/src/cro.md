# Introduction

## Purpose

Provides a **crossover** from a solution
to the **convex quadratic programming problem**
$\mbox{minimize}\;\; q(x) = \frac{1}{2} x^T H x + g^T x + f $
\[
minimize q(x) := 1/2 x^T H x + g^T x + f
\]
subject to the general linear constraints
$c_i^l\leqa_i^Tx\leq c_i^u, \;\;\; i = 1, \ldots , m,$
\n
 c_i^l \[<=] a_i^Tx \[<=] c_i^u, i = 1, ... , m,
\n
and the simple bound constraints
$x_j^l\leqx_j \leq x_j^u, \;\;\; j = 1, \ldots , n,$
\n
 x_j^l \[<=] x_j \[<=] x_j^u, j = 1, ... , n,
\n
found by an interior-point method to one in which the
**matrix of defining active constraints/variables is of full rank.**
Here, the $n$ by $n$ symmetric, positive-semi-definite matrix
$H$, the vectors $g$, $a_i$, $c^l$, $c^u$, $x^l$,
$x^u$, the scalar $f$ are given. In addition a solution $x$ along
with optimal Lagrange multipliers $y$ for the general constraints
and dual variables $z$ for the simple bounds must be provided
(see Section~\ref{galmethod}). These will be adjusted as necessary.
Any of the constraint bounds $c_i^l$, $c_i^u$, $x_j^l$
and $x_j^u$ may be infinite.
Full advantage is taken of any zero coefficients in the matrix $H$
or the matrix $A$ of vectors $a_i$.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

August 2010, C interface January 2022.

## Terminology

Any required solution $x$ necessarily satisfies
the primal optimality conditions
$\mbox{(1a) $\hspace{66mm} A x = c\hspace{66mm}$}$
\n
(1a) A x = c
\n
and
$\mbox{(1b) $\hspace{52mm} c^l \leq c \leq c^u, \;\; x^l \leq x \leq x^u,\hspace{52mm}$} $
\n
(1b) c^l \[<=] c \[<=] c^u, x^l \[<=] x \[<=] x^u,
\n
the dual optimality conditions
$\mbox{(2a) $\hspace{58mm} H x + g = A^T y + z\hspace{58mm}$}$
\n
(2a) H x + g = A^T y + z
\n
where
$\mbox{(2b) $\hspace{24mm} y = y^l + y^u, \;\; z = z^l + z^u, \,\,
 y^l \geq 0 , \;\;y^u \leq 0 , \;\;
 z^l \geq 0 \;\; \mbox{and} \;\; z^u \leq 0,\hspace{24mm}$} $
\n
 (2b) y = y^l + y^u, z = z^l + z^u, y^l \[>=] 0, y^u \[<=] 0,
z^l \[>=] 0 and z^u \[<=] 0,
\n
and the complementary slackness conditions
$\mbox{(3) $\hspace{12mm}
( A x - c^l )^T y^l = 0,\;\;( A x - c^u )^T y^u = 0,\;\;
(x -x^l )^T z^l = 0 \;\;\mbox{and} \;\; (x -x^u )^T z^u = 0,\hspace{12mm} $}$
\n
(3) (A x - c^l)^T y^l = 0, (A x - c^u)^T y^u = 0,
(x -x^l)^T z^l = 0 and (x -x^u)^T z^u = 0,
\n
where the vectors $y$ and $z$ are known as the Lagrange multipliers
for the general linear constraints, and the dual variables for the bounds,
respectively, and where the vector inequalities hold component-wise.

## Method

Denote the active constraints by $A_A x = c_A$ and the active bounds by
$I_A x = x_A$. Then any optimal solution satisfies the linear system
$\left(\begin{array}{ccc}H & - A_A^T & - I^T_A \\ A_A & 0 & 0 \\ I_A & 0 & 0 \end{array}\right) \left(\begin{array}{c}x \\ y_A \\ z_A\end{array}\right) =
\left(\begin{array}{c}- g \\ c_A \\ x_A\end{array}\right).$
\n
 ( H - A_A^T - I_A^T ) (x) ( - g )
 ( A_A 0 0 ) ( y_A ) = ( c_A ),
 ( I_A 0 0 ) ( z_A ) ( x_A )
\n
where $y_A$ and $z_A$ are the corresponding active Lagrange
multipliers and dual variables respectively. Consequently the difference
between any two solutions $(\Delta x, \Delta y, \Delta z)$ must satisfy
$\mbox{(4)}\;\; \left(\begin{array}{ccc}H & - A_A^T & - I^T_A \\ A_A & 0 & 0 \\ I_A & 0 & 0 \end{array}\right) \left(\begin{array}{c}\Delta x \\ \Delta y_A \\ \Delta z_A\end{array}\right) = 0.$
\n
( H - A_A^T - I_A^T ) (Delta x)
(4) ( A_A 0 0 ) ( Delta y_A ) = 0,
( I_A 0 0 ) ( Delta z_A )
\n
Thus there can only be multiple solution if the coefficient matrix $K$
of (4) is singular. The algorithm used in CRO
exploits this. The matrix $K$ is checked for singularity
using the GALAHAD package ULS. If $K$ is
non singular, the solution is unique and the solution input by the user
provides a linearly independent active set. Otherwise $K$ is singular,
and partitions $A_A^T = ( A_B^T \;\; A_N^T)$ and
$I_A^T = ( I_B^T \;\; I_N^T)$ are found so that
$\left(\begin{array}{ccc}H & - A_B^T & - I_B^T \\ A_B & 0 & 0 \\ I_B & 0 & 0 \end{array}\right)$
\n
 ( H - A_B^T - I_B^T )
 ( A_B 0 0 )
 ( I_B 0 0 )
\n
is non-singular and the **non-basic** constraints $A_N^T$
and $I_N^T$ are linearly dependent on the **basic** ones
$( A_B^T \;\; I_B^T)$. In this case (4) is equivalent to
$\mbox{(5)}\;\; \left(\begin{array}{ccc}H & - A_B^T & - I_B^T \\ A_B & 0 & 0 \\ I_B & 0 & 0 \end{array}\right) = \left(\begin{array}{c}A_N^T \\ 0 \\ 0\end{array}\right) \Delta y_N + \left(\begin{array}{c}I_N^T \\ 0 \\ 0\end{array}\right) \Delta z_N$
\n
( H - A_B^T - I_B^T ) (Delta x)
(5) ( A_B 0 0 ) ( Delta y_A ) =
( I_B 0 0 ) ( Delta z_A )

( A_N^T ) ( I_N^T )
( 0 ) Delta y_N + ( 0 ) Delta z_N.
( 0 ) ( 0 )
\n
Thus, starting from the user's $(x, y, z)$
and with a factorization of the coefficient matrix of (5)
found by the GALAHAD package SLS, the alternative solution
$(x + \alpha x, y + \alpha y, z + \alpha z)$,
featuring
$(\Delta x, \Delta y_B, \Delta z_B)$
from (5)in which successively one of the components of $\Delta y_N$
and $\Delta z_N$ in turn is non zero, is taken.
The scalar $\alpha$ at each stage
is chosen to be the largest possible that guarantees (2.b);
this may happen when a non-basic multiplier/dual variable reaches zero,
in which case the corresponding constraint is disregarded, or when this
happens for a basic multiplier/dual variable, in which case this constraint is
exchanged with the non-basic one under consideration and disregarded.
The latter corresponds to changing the basic-non-basic partition
in (5), and subsequent solutions may be found by updating
the factorization of the coefficient matrix in (5)
following the basic-non-basic swap using the GALAHAD package SCU.

## Reference

## Call order

To solve a given problem, functions from the cro package must be called
in the following order:

- cro\_initialize - provide default control parameters and set up initial data structures
- cro\_read\_specfile (optional) - override control values by reading replacement values from a file
- cro_crossover_solution - move from a primal-dual soution
to a full rank one
- cro\_terminate - deallocate data structures

## fdc_array_indexing Array indexing

Both C-style (0 based)and fortran-style (1-based) indexing is allowed.
Choose control.f_indexing as false for C style and true for
fortran style; add 1 to input integer arrays if fortran-style indexing is
used, and beware that return integer arrays will adhere to this.

