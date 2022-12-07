# Introduction

## Purpose

Given a sequence of vectors
$\{s_k\}$ and $\{y_k\}$ \mbox{and scale factors} $\{\delta_k\}$,
{s<sub>k</sub>} and {y<sub>k</sub>} and scalars {&#948<sub>k</sub>},
{s_k} and {y_k} and scalars {delta_k},
**obtain the product of a limited-memory secant
approximation $H_k$ (or its inverse) with a given vector**,
using one of a variety of well-established formulae.

Currently, only the control and inform parameters are exposed;
these are provided and used by other GALAHAD packages with C interfaces.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

July 2014, C interface January 2022.

## Method

Given a sequence of vectors
$\{s_k\}$ and $\{y_k\}$ \mbox{and scale factors} $\{\delta_k\}$,
{s<sub>k</sub>} and {y<sub>k</sub>} and scalars {&#948<sub>k</sub>},
{s_k} and {y_k} and scalars {delta_k},
a limited-memory secant approximation $H_k$ is chosen
 so that $H_{\max(k-m,0)} = \delta_k I$, $H_{k-j} s_{k-j} = y_{k-j}$
 and $\| H_{k-j+1} - H_{k-j}\|$ is “small” for
$j = \min(k-1,m-1), \ldots, 0$.
Different ways of quantifying “small” distinguish different methods,
but the crucial observation is that it is possible to construct
$H_k$ quickly from ${s_k}$, ${y_k}$ and $\delta_k$,
 and to apply it and its inverseto a given vector $v$.
 It is also possible to apply similar formulae to the “shifted” matrix
$H_k + \lambda_k I$ that occurs in trust-region methods.

## Reference

The basic methods are those given by

R. H. Byrd, J. Nocedal and R. B. Schnabel (1994)
Representations of quasi-Newton matrices and their use in
limited memory methods.
Mathenatical Programming, **63(2)** 129-156,

with obvious extensions.
