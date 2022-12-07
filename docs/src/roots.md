# Introduction

## Purpose

Use classical formulae together with Newton’s method to find all the real
roots of a real polynomial.

Currently, only the control and inform parameters are exposed;
these are provided and used by other GALAHAD packages with C interfaces.

## Authors

N. I. M. Gould, STFC-Rutherford Appleton Laboratory, England.

C interface, additionally J. Fowkes, STFC-Rutherford Appleton Laboratory.

Julia interface, additionally A. Montoison and D. Orban, Polytechnique Montréal.

## Originally released

April 2005, C interface January 2022.

## Method

Littlewood and Ferrari's algorithms are used to find estimates of the
real roots of cubic and quartic polynomials, respectively; a
stabilized version of the well-known formula is used in the quadratic
case. Newton's method is used to further refine the computed roots if
necessary. Madsen and Reid's method is used for polynomials whose
degree exceeds four.
