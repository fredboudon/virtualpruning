from math import exp
import numpy as np
from numpy.random import binomial, poisson, uniform, normal, negative_binomial
import sys


def realisation(mean, sd = None, minval = 0, maxval = sys.maxsize, rfunc = normal):
    assert maxval > minval
    if sd is None:
        val = rfunc(mean, 1)
    else:
        val = rfunc(mean, sd, 1)
    count = 0
    while (val < minval) or (val > maxval):
        count += 1
        if count >= 1000:
            raise ValueError(mean, sd, maxval, minval)
        if sd is None:
            val = rfunc(mean, 1)
        else:
            val = rfunc(mean, sd, 1)
    return val

def binomial_proba_from_latent(latent):
    return exp(latent)/(1+exp(latent))

def binomial_proba(intercept, coefs, factors):
    latent = intercept + np.array(np.array(coefs) * np.array(factors)).sum()
    return binomial_proba_from_latent(latent)

def binomial_realization(proba):
    return bool( binomial(1,proba) )

def poisson_proba_from_latent(latent):
    return exp(latent)

def poisson_proba(intercept, coefs, factors):
    latent = intercept + np.array(np.array(coefs) * np.array(factors)).sum()
    return poisson_proba_from_latent(latent)

def poisson_realization(proba, maxval = sys.maxsize, minval = 0):
    return int(realisation(proba,None,minval,maxval,poisson))


def negativebinomial_mu_from_latent(latent):
    return exp(latent)

def negativebinomial_realization(mu, theta, maxval = sys.float_info.max, minval = sys.float_info.min):
    proba = theta/(theta + mu)
    if proba < 0 or proba > 1:
        raise ValueError(proba, theta, mu)
    return int(realisation(theta,proba,minval,maxval,negative_binomial))


def normal_realization(mean, sd, maxval = sys.float_info.max, minval = sys.float_info.min):
    return int(realisation(mean,sd,minval,maxval,normal))

