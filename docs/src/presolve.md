# Introduction

## Purpose

Presolving aims to **improve the formulation of a given optimization
problem by applying a sequence of simple transformations**, and thereby
to produce a \a reduced problem in a \a standard \a form that should be
simpler to solve.This reduced problem may then be passed to an
appropriate solver.Once the reduced problem has been solved, it is
then \a restored to recover the solution for the original formulation.

This package applies presolving techniques to a **linear**
$\mbox{minimize}\;\; l(x) = g^T x + f $
\n
minimize l(x) := g^T x + f
\n
or **quadratic program**
$\mbox{minimize}\;\; q(x) = \frac{1}{2} x^T H x + g^T x + f $
\n
minimize q(x) := 1/2 x^T H x + g^T x + f
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
where the $n$ by $n$ symmetric matrix $H$,
the vectors $g$, $a_i$, $c^l$, $c^u$, $x^l$,
$x^u$ and the scalar $f$ are given.
Any of the constraint bounds $c_i^l$, $c_i^u$,
$x_j^l$ and $x_j^u$ may be infinite.

In addition, bounds on the Lagrange multipliers $y$ associated with
the general linear constraints and on the dual variables $z$ associated
with the simple bound constraints
$ y_{i}^{l}\leqy_{i}\leqy_{i}^{u}, \;\;\;i = 1, \ldots , m,$
\n
 y_j^i \[<=] y_i \[<=] y_i^u, i = 1, ... , m,
\n
and
$z_{i}^{l}\leqz_{i}\leqz_{i}^{u}, \;\;\;i = 1, \ldots , n,$
\n
 z_j^l \[<=] z_j \[<=] z_j^u, j = 1, ... , n,
\n
are also provided, where the $m$-dimensional vectors $y^l$ and
$y^u$, as well as the $n$-dimensional vectors $x^l$ and $x^u$
are given.Any component of $c^l$, $c^u$, $x^l$, $x^u$,
$y^l$, $y^u$, $z^l$ or $z^u$ may be infinite.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England
and Ph. L. Toint, University of Namur, Belgium

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

March 2002, C interface March 2022.

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
for2 the general linear constraints, and the dual variables for the bounds,
respectively, and where the vector inequalities hold component-wise.

## Method

The purpose of presolving is to exploit these equations in order to reduce
the problem to the standard form defined as follows:
- The variables are ordered so that their bounds appear in the order
$\begin{array}{lccccc}
\mbox{free}&&& x &&\\
\mbox{non-negativity}& 0& \leq & x &&\\
\mbox{lower} & x^l & \leq & x & &\\
\mbox{range} & x^l & \leq & x & \leq & x^u\\
\mbox{upper} &&& x & \leq & x^u \\
\mbox{non-positivity}&&& x & \leq &0
\end{array}$
\n
free x
non-negativity 0<= x
lower x^l <= x
range x^l <= x<= x^u
upperx<= x^u
non-positivity x<=0
\n
Fixed variables are removed. Within each category, the variables
are further ordered so that those with non-zero diagonal Hessian
entries occur before the remainder.
- The constraints are ordered so that their bounds appear in the order
$\begin{array}{lccccc}
\mbox{non-negativity} & 0& \leq & A x &&\\
\mbox{equality} & c^l & =& A x && \\
\mbox{lower} & c^l & \leq & A x && \\
\mbox{range} & c^l & \leq & A x & \leq & c^u\\
\mbox{upper} &&& A x & \leq & c^u \\
\mbox{non-positivity} &&& A x & \leq & 0\\
\end{array}$
\n
non-negativity 0<= A x
equalityc^l= A x
lower c^l <= A x
range c^l <= A x <= c^u
upperA x <= c^u
non-positivity A x <=0
\n
Free constraints are removed.
- In addition, constraints may be removed or bounds tightened, to reduce the
size of the feasible region or simplify the problem if this is possible, and
bounds may be tightened on the dual variables and the multipliers
associatedwith the problem.

The presolving algorithm proceeds by applying a (potentially long) series of
simple transformations to the problem, each transformation introducing a
further simplification of the problem. These involve the removal of empty and
singleton rows, the removal of redundant and forcing primal constraints, the
tightening of primal and dual bounds, the exploitation of linear singleton,
linear doubleton and linearly unconstrained columns, the merging dependent
variables, row sparsification and split equalities. Transformations are
applied in successive passes, each pass involving the following actions:

-# remove empty and singletons rows,
-# try to eliminate variables that are linearly unconstrained,
-# attempt to exploit the presence of linear singleton columns,
-# attempt to exploit the presence of linear doubleton columns,
-# complete the analysis of the dual constraints,
-# remove empty and singletons rows,
-# possibly remove dependent variables,
-# analyze the primal constraints,
-# try to make $A$ sparser by combining its rows,
-# check the current status of the variables, dual variables and multipliers.

All these transformations are applied to the structure of the original
problem, which is only permuted to standard form after all transformations are
completed. <em>Note that the Hessian and Jacobian of the resulting reduced
problem are always stored in sparse row-wise format.</em> The reduced problem
is then solved by a quadratic or linear programming solver, thus ensuring
sufficiently small primal-dual feasibility and complementarity. Finally, the
solution of the simplified problem is re-translated in the
variables/constraints/format of the original problem formulation by a
\a restoration phase.

If the number of problem transformations exceeds
\p control.transf_buffer_size,the transformation buffer size,
then they are saved in a “history” file, whose
name may be chosen by specifying the control.transf_file_name control
parameter,When this is the case, this file
is subsequently reread by \p presolve_restore_solution. It must not be
altered by the user.

Overall, the presolving process follows one of the two sequences:

$\fbox{initialize} \rightarrow \left[ \fbox{apply transformations}
 \rightarrow \mbox{(solve problem)}
 \rightarrow \fbox{restore} \right] \rightarrow \fbox{terminate}$
or
$\fbox{initialize} \rightarrow \left[ \fbox{read specfile}
 \rightarrow \fbox{apply transformations}
 \rightarrow \mbox{(solve problem)}
 \rightarrow \fbox{restore} \right] \rightarrow \fbox{terminate}$
 (ignore garbled doxygen phrase)
\n
 --------------[-------------------------
 | initialize | -> [ | apply transformations | -> (solve problem) ->
 --------------[-------------------------
----------- ]-------------
| restore | ] -> | terminate |
----------- ]-------------
 or
 --------------[ ------------------------------------------
 | initialize | -> [ | read specfile | -> | apply transformations | ->
 --------------[ ------------------------------------------
 ----------- ]-------------
(solve problem) -> | restore | ] -> | terminate |
 ----------- ]-------------
\n

where the procedure's control parameter may be modified by
reading the specfile, and where (solve problem) indicates that the reduced
 problem is solved. Each of the
“boxed” steps in these sequences corresponds to calling a specific
routine of the package In the diagrams above, brackated subsequence of
steps means that they can be repeated with problem having the same
structure. The value of the \p problem.new_problem_structure
must be true on entry of \p presolve_apply_to_problem on the
first time it is used in this repeated subsequence. Such a subsequence must
be terminated by a call to\p presolve\_terminate before presolving is
applied to a problem with a different structure.

Note that the values of the multipliers and dual variables (and thus of
their respective bounds) depend on the functional form assumed for the
Lagrangian function associated with the problem.This form is given by
$L(x,y,z) = q x) - y\_{sign} * y^T (Ax-c) - z\_{sign} * z,$
(considering only active constraints $A x = c$), where the parameters
y_{sign} and z_{sign} are +1 or -1 and can be chosen by the user.
Thus, if $y_{sign}$ = +1, the multipliers associated to active constraints
originally posed as inequalities are non-negative if the inequality is a lower
bound and non-positive if it is an upper bound. Obvioulsy they are not
constrained in sign for constraints originally posed as equalities. These
sign conventions are reversed if $y_{sign}$ = -1.
Similarly, if $z_{sign}$ = +1}, the dual variables associated to active
bounds are non-negative if the original bound is an lower bound, non-positive
if it is an upper bound, or unconstrained in sign if the variables is fixed;
and this convention is reversed in $z\_{sign}$ = -1}. The values of
$z_{sign}$ and $y_{sign}$ may be chosen by setting the corresponding
components of the \p control structure to \p 1 or \p -1.

## Reference

The algorithm is described in more detail in

N. I. M. Gould and Ph. L. Toint (2004).
Presolving for quadratic programming.
Mathematical Programming **100**(1), pp 95--132.

## Call order

To solve a given problem, functions from the presolve package must be called
in the following order:

- presolve\_initialize - provide default control parameters and
set up initial data structures
- presolve\_read\_specfile (optional) - override control values
by reading replacement values from a file
- presolve\_import_problem - import the problem data and report
the dimensions of the transformed problem
- presolve_transform_problem - apply the presolve algorithm
to transform the data
- presolve_restore_solution - restore the solution from
 that of the transformed problem
- presolve\_information (optional) - recover information about
the solution and solution process
- presolve\_terminate - deallocate data structures

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

##  Symmetric matrix storage formats

Likewise, the symmetric $n$ by $n$ objective Hessian matrix
$H$ may be presented
and stored in a variety of formats. But crucially symmetry is exploited
by only storing values from the lower triangular part
(i.e, those entries that lie on or below the leading diagonal).

### Dense storage format

The matrix $H$ is stored as a compactdense matrix by rows, that is,
the values of the entries of each row in turn are
stored in order within an appropriate real one-dimensional array.
Since $H$ is symmetric, only the lower triangular part (that is the part
$h_{ij}$ for $0 \leq j \leq i \leq n-1$) need be held.
In this case the lower triangle should be stored by rows, that is
component $i \ast i / 2 + j$of the storage array H_val
will hold the value $h_{ij}$ (and, by symmetry, $h_{ji}$)
for $0 \leq j \leq i \leq n-1$.

###  Sparse co-ordinate storage format

Only the nonzero entries of the matrices are stored.
For the $l$-th entry, $0 \leq l \leq ne-1$, of $H$,
its row index i, column index j
and value $h_{ij}$, $0 \leq j \leq i \leq n-1$,are stored as
the $l$-th components of the integer arrays H_row and
H_col and real array H_val, respectively, while the number of nonzeros
is recorded as H_ne = $ne$.
Note that only the entries in the lower triangle should be stored.

###  Sparse row-wise storage format

Again only the nonzero entries are stored, but this time
they are ordered so that those in row i appear directly before those
in row i+1. For the i-th row of $H$ the i-th component of the
integer array H_ptr holds the position of the first entry in this row,
while H_ptr(n) holds the total number of entries plus one.
The column indices j, $0 \leq j \leq i$, and values
$h_{ij}$ of theentries in the i-th row are stored in components
l = H_ptr(i), $\ldots$, H_ptr(i+1)-1 of the
integer array H_col, and real array H_val, respectively.
Note that as before only the entries in the lower triangle should be stored.
For sparse matrices, this scheme almost always requires less storage than
its predecessor.

### symmetric\_matrix_diagonal Diagonal storage format

If $H$ is diagonal (i.e., $H_{ij} = 0$ for all
$0 \leq i \neq j \leq n-1$) only the diagonals entries
$H_{ii}$, $0 \leq i \leq n-1$ need
be stored, and the first n components of the array H_val may be
used for the purpose.

### symmetric\_matrix_scaled_identity Multiples of the identity storage format

If $H$ is a multiple of the identity matrix, (i.e., $H = \alpha I$
where $I$ is the n by n identity matrix and $\alpha$ is a scalar),
it suffices to store $\alpha$ as the first component of H_val.

### symmetric\_matrix_identity The identity matrix format

If $H$ is the identity matrix, no values need be stored.

### symmetric\_matrix_zero The zero matrix format

The same is true if $H$ is the zero matrix.

