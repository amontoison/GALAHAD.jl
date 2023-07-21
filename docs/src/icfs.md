# Introduction

## Purpose

Given a symmetric matrix $\bmA$, this package
** computes a symmetric, positive-definite approximation
$\bmL \bmL^T$ using an incomplete Cholesky factorization**;
the resulting matrix $\bmL$ is lower triangular.
Subsequently, the solution $\bmx$ to the either of the linear systems
$\bmL \bmx = \bmb$ and $\bmL^T \bmx = \bmb$
may be found for a given vector $\bmb$.

## Authors

C.-J, Lin and J. J. Moré, Argonne National Laboratory,

C interface, additionally N. I. M. Gould and J. Fowkes, STFC-Rutherford
Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

 ## Originally released

May 1998, C interface December 2022.
