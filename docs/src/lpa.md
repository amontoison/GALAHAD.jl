# Introduction

## Purpose

This package uses the ** simplex method**
to solve the **linear programming problem**
$\mbox{minimize}\;\; q(x) = g^T x + f $
\n
minimize q(x) := g^T x + f
\n
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
where the vectors $g$, $w$, $x^{0}$,
$a_i$, $c^l$, $c^u$, $x^l$,
$x^u$ and the scalar $f$ are given.
Any of the constraint bounds $c_i^l$, $c_i^u$,
$x_j^l$ and $x_j^u$ may be infinite.
Full advantage is taken of any zero coefficients in the matrix
$A$ whose rows are the transposes of the vectors $a_i$.

**N.B.** The package is simply a sophisticated interface to the
HSL package LA04, and requires that a user has obtained the latter.
** LA04 is not included in GALAHAD**
but is available without charge to recognised academics, see
http://www.hsl.rl.ac.uk/catalogue/la04.html. If LA04
is unavailable, the GALAHAD interior-point linear programming
package LPB is an alternative.

## Authors

N. I. M. Gould and J. K. Reid, STFC-Rutherford Appleton Laboratory,
England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

October 2018, C interface September 2021.

## Terminology

The required solution $x$ necessarily satisfies
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
$\mbox{(2a) $\hspace{3mm} g = A^T y + z$}$
\n
(2a) g = A^T y + z
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
where the vectors $y$ and $z$ are
known as the Lagrange multipliers for
the general linear constraints, and the dual variables for the bounds,
respectively, and where the vector inequalities hold component-wise.

The so-called dual to this problem is another linear program
$- \mbox{minimize} \;\; c^{lT} y^l + c^{uT} y^u + x^{lT} z^l + x^{uT} z^u + f \;\; \mbox{subject to the constraints (2a) and (2b)}$
\n
- minimize c^{lT} y^l + c^{uT} y^u + x^{lT} z^l + x^{uT} z^u + f
subject to the constraints (2a) and (2b)
\n
that uses the same data. The solution to the two problems, it is exists,
is the same, but if one is infeasible, the other is unbounded. It can be
more efficient to solve the dual, particularly if $m$ is much larger
than $n$.

## Method

The bulk of the work is peformed by the HSL package LA04. The
main subbroutine from this package requires that the input problem
be transformed into the “standard form”
\[\begin{array}{rl}\mbox{minimize} & g^{\prime T} x^{\prime} \\(4)\;\; \mbox{subject to} & A^{\prime} x^{\prime} = b \\ &l_i \leq x^{\prime}_i \leq u_i, \;\;(i\leq k) \\ \mbox{and} & x^{\prime}_l \geq 0, \;\; (i \geq l) \end{array} \]
$ \begin{array}{rl}\mbox{minimize} & g^{\prime T} x^{\prime} \\(4)\;\; \mbox{subject to} & A^{\prime} x^{\prime} = b \\ &l_i \leq x^{\prime}_i \leq u_i, \;\;(i\leq k) \\ \mbox{and} & x^{\prime}_l \geq 0, \;\; (i \geq l) \end{array} $
\n
 minimize g'^T x'
(4)subject to A' x' = b
l_i <= x'_i <= u_i, for i <= k
and x'_l >= 0, for i >= l
\n
by introducing slack an surpulus variables, reordering and
removing fixed variables and free constraints. The resulting
problem involves $n'$ unknowns and $m'$ general constraints.
In order to deal with the possibility that the general constraints
are inconsistent or not of full rank,
LA04 introduces additional “artifical” variables $v$ and replaces
the constraints of (4) by
$(5) \;\; A' x' + v = b$
and gradually encourages $v$ to zero as a first solution phase.

Once a selection of $m'$ independent **non-basic** variables
is made, the constraints (5) determine the remaining $m'$
dependent **basic** variables. The **simplex method** is a
scheme for systematically adjusting the choice of basic and non-basic
variables until a set which defines an optimal solution of (4) is
obtained. Each iteration of the simplex method requires the solution
of a number of sets of linear equations whose coefficient matrix is
the **basis** matrix $B$, made up of the columns of
$[A'$$I]$ corresponding to the basic variables, or its transpose
$B^T$. As the basis matrices for consecutive iterations are
closely related, it is normally advantageous to update (rather than
recompute) their factorizations as the computation proceeds.If an
initial basis is not provided by the user, a set of basic variables
which provide a (permuted) triangular basis matrix is found by the
simple crash algorithm of Gould and Reid (1989), and initial
steepest-edge weights are calculated.

Phases one (finding a feasible solution) and two (solving (4)
of the simplex method are applied, as appropriate, with the choice of
entering variable as described by Goldfarb and Reid (1977) and the
choice of leaving variable as proposed by Harris (1973).
Refactorizations of the basis matrix are performed whenever doing so
will reduce the average iteration time or there is insufficient memory
for its factors.The reduced cost for the entering variable is
computed afresh. If it is found to be of a different sign from the
recurred value or more than 10\% different in magnitude, a fresh
computation of all the reduced costs is performed.Details of the
factorization and updating procedures are given by Reid (1982).
Iterative refinement is encouraged for the basic solution and for the
reduced costs after each factorization of the basis matrix and when
they are recomputed at the end of phase 1.

## References

D. Goldfarb and J. K. Reid (1977).
A practicable steepest-edge simplex algorithm.
Mathematical Programming **12** 361-371.

N. I. M. Gould and J. K. Reid (1989)
New crash procedures for large systems of linear constraints.
Mathematical Programming **45** 475-501.

P. M. J. Harris (1973).
Pivot selection methods of the Devex LP code.
Mathematical Programming **5** 1-28.

J. K. Reid (1982)
A sparsity-exploiting variant of the Bartels-Golub
decomposition for linear-programming bases.
Mathematical Programming **24** 55-69.

## Call order

To solve a given problem, functions from the lpa package must be called
in the following order:

- lpa\_initialize - provide default control parameters and set up initial data structures
- lpa\_read\_specfile (optional) - override control values by reading replacement values from a file
- lpa\_import - set up problem data structures and fixed values
- lpa\_reset\_control (optional) - possibly change control parameters if a sequence of problems are being solved
- lpa\_solve_lp - solve the linear program
- lpa\_information (optional) - recover information about the solution and solution process
- lpa\_terminate - deallocate data structures

##  Unsymmetric matrix storage formats

The unsymmetric $m$ by $n$ constraint matrix $A$ may be presented
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
