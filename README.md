# Poisson binomial model of dispersal

This repository contains R scripts associated with the paper:

XXXXXX

The main purpose of these scripts is to estimate the parameters of 16 plausible dispersal models given experimental data involving flour beetles in a 2-patch metacommunity.

We assume that the number of beetles in both patches (1 and 2) prior to dispersal is known and model the number of beetles in patch 1 following dispersal as the random variable Z = X + Y, where X is the number of beetles that did not disperse from patch 1 and Y is the number that did disperse from patch 2. Here, X and Y are treated as binomial random variables with probabilities p1 and p2. Therefore, Z is a Poisson binomial random variable. The probability mass function for this distribution is known and we use the package [poibin](https://cran.r-project.org/web/packages/poibin/index.html) to solve it.

The probabilities p1 and p2 are calculated from the 16 dispersal models.

*Note: Assuming a closed system, the number of beetles after dispersal in patch 1 exactly determines the number of beetles in patch 2. so it only necessary to model the number of beetles in a single patch*

## nll.R

This function calculates the negative log-likelihood of the data given the parameters and the dispersal model being tested.

## models/

<<<<<<< HEAD
This folder contains the 16 dispersal models considered in the main text. 

In addition, it includes an alternate, mathematically equivalent version of the full model ("modelalt5") that conceptualizes the system as one in which densities of beetles equilibrate between the two patches. See the supplementary text for a description of this model and its interpretation.
=======
This folder contains the 16 dispersal models; that is, the 16 different ways of calculating p1 and p2.
>>>>>>> c13bf3425eae1bd4ab1f261e0272d0e68ad2d2c9

## analysis-fitting.R

This script uses base function "optim" and likelihood function "nll.R" to determine which of the models has the best fit to the data.

## analysis-bootstrap.R

This script uses non-parametric bootstrapping to estimate the error in the estimates of the best-fit model.

## check-nll.R

This script can be used to check the statistical power of our likelihood function. It estimates the parameter values of the different dispersal models given simulated data

## simulations/

This folder contains scripts for simulating each of the 16 dispersal models

## data/

This folder contains the data used for optimization and boostrapping

