from numpy.random import binomial, poisson, uniform, normal

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
    assert maxval > minval
    val = int( poisson(proba, 1) )
    count = 0
    while (val < minval) or (val > maxval):
        count += 1
        if count >= 1000:
            raise ValueError(proba, maxval, minval)
        val = int( poisson(proba, 1) )
    return val

def normal_realization(mean, sd, maxval = sys.float_info.max, minval = sys.float_info.min):
    assert maxval > minval
    val = normal(mean, sd)
    count = 0
    while (val < minval) or (val > maxval):
        count += 1
        if count >= 1000:
            raise ValueError(mean, sd, maxval, minval, val)
        val = normal(mean, sd)
    return val
