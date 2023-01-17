import numpy as np
from regrowth_base import *
import mangoG3 as mg
from pruning import intensity_level

def burst_pruned(diameter, intensity, Zeta_mean):
    if np.isnan(Zeta_mean):
        Zeta_mean = 0
    intercept = -9.85 
    latent = intercept + 0.39 * diameter + 7.60 * intensity + 4.39 * Zeta_mean
    probavalue = binomial_proba_from_latent(latent)

    return binomial_realization(probavalue)

def burst_unpruned(diameter, intensity, Zeta_mean):
    if np.isnan(Zeta_mean):
        Zeta_mean = 0
    intercept = -7.62 
    latent = intercept + 0.32 * diameter + 5.71 * intensity + 4.28 * Zeta_mean

    probavalue = binomial_proba_from_latent(latent)
    return binomial_realization(probavalue)

################################################################################# Ajout 2022 Intensité du débourrement

def nb_daughter_pruned(diameter, intensity, Zeta_mean):
    import numpy as np
    if np.isnan(zeta_mean):
        Zeta_mean = 0
    intercept = - 4.12
    probavalue = poisson_proba_from_latent(intercept + 0.32 * diameter 
                                                     + 1.19 * intensity 
                                                     + 4.12 * Zeta_mean
                                                     -0.25 * diameter * Zeta_mean)
    try:
        return poisson_realization( probavalue, 10)+1  
    except ValueError as ve:
        print(diameter)
        raise ve

def nb_daughter_unpruned(diameter, intensity, Zeta_mean):
    intercept = -8.15
    probavalue = poisson_proba_from_latent(intercept + 0.23 * diameter 
                                                     + 2.12 * intensity 
                                                     + 6.30 * Zeta_mean)
    try:
        return poisson_realization( probavalue, 10)+1 
    except ValueError as ve:
        print(intensity, diameter, Zeta_mean)
        raise ve

################################################################################# Ajout 2022 Dynamique du débourrement
def flush_selection(intensity):
    intercept = 2.44
    latent = intercept + 6.02 * intensity
    firstflushproba = binomial_proba_from_latent(latent)       
    return binomial_realization(firstflushproba)

def burstdelay(intensity):
    
    firstflush = flush_selection(intensity)
    if firstflush :
        minval = 9
        maxval = 34
        mean = 24.3
        sd = 2.91
    else:
        minval = 38
        maxval = 41
        mean = 40
        sd = 0.73
    return  normal_realization(mean, sd, maxval, minval)
 #################################################################################


from datetime import date, timedelta
from importlib import reload
import mortality ; reload(mortality)
from mortality import gu_mortality_post_regrowth

def growth_pruned_gu(mtg, vid, intensity, Zeta_mean, pruningdate, mortalityenabled = True):
    if np.isnan(Zeta_mean):
        Zeta_mean = 0
    diameter = mg.get_gu_diameter(mtg, vid)
    veggrowth = burst_pruned(diameter, intensity, Zeta_mean)
    if veggrowth:
        nbdaughter = nb_daughter_pruned(diameter, intensity, Zeta_mean)
        assert nbdaughter >= 1
        if mortalityenabled:
            nbdaughter = gu_mortality_post_regrowth(mtg, vid, nbdaughter)
        if nbdaughter > 0:
            ilevel = intensity_level(intensity)
            totalleafarea = total_leafarea_pruned(diameter)
            individualleafarea = individual_leafarea_pruned(ilevel)
            burstdate = pruningdate + timedelta(days = burstdelay(intensity))
            return create_daughters(mtg, vid, 0, nbdaughter, burstdate, totalleafarea, individualleafarea)       


def growth_unpruned_gu(mtg, vid, intensity, Zeta_mean, pruningdate, mortalityenabled = True):
    if np.isnan(Zeta_mean):
        Zeta_mean = 0
    diameter = mg.get_gu_diameter(mtg, vid)
    veggrowth = burst_unpruned(diameter, intensity, Zeta_mean)
    if veggrowth:
        nbdaughter = nb_daughter_unpruned(diameter, intensity, Zeta_mean)
        assert nbdaughter >= 1
        if mortalityenabled:
            nbdaughter = gu_mortality_post_regrowth(mtg, vid, nbdaughter)
        if nbdaughter > 0:
            apical_bud = binomial_realization(0.52) # On le laisse ?
            ilevel = intensity_level(intensity)
            totalleafarea = total_leafarea_unpruned(ilevel, diameter, apical_bud)
            individualleafarea = individual_leafarea_unpruned(ilevel)
            burstdate = pruningdate + timedelta(days = burstdelay(intensity))
            return create_daughters(mtg, vid, int(apical_bud), nbdaughter-int(apical_bud), burstdate, totalleafarea, individualleafarea)



def growth(mtg, Zeta_mean = None, intensity = None, 
           pruningdate = date(2021,2,24), maxdiamunpruned = 10, mortalityenabled = True):
    if intensity is None:
        from pruning import continuous_intensity_from_pruned
        intensity = continuous_intensity_from_pruned(mtg)
    from copy import deepcopy
    newmtg = deepcopy(mtg)
    pruned = mtg.property('pruned')
    terminals = mg.get_all_terminal_gus(mtg)
    nbterminals = len(terminals)
    newids = []
    print("Should examine", nbterminals, "terminal GUs.")
    nbpruned, nbunpruned, nbignored = 0,0,0
    for vid in terminals:
        if vid in pruned:
            lnewids = growth_pruned_gu(newmtg, vid, 
                                       intensity, 
                                       Zeta_mean.get(vid,0),
                                       pruningdate,
                                       mortalityenabled=mortalityenabled)
            nbpruned += 1
        elif mg.get_gu_diameter(mtg, vid) <= maxdiamunpruned and not 'A' in mtg.property('Taille').get(vid,''):
            lnewids = growth_unpruned_gu(newmtg, vid, 
                                         intensity, 
                                         Zeta_mean.get(vid,0), 
                                         pruningdate,
                                         mortalityenabled=mortalityenabled)
            nbunpruned += 1
        else :
            nbignored += 1
        if lnewids:
            newids += lnewids
    print("Processed", nbpruned, "pruned terminal GU and", nbunpruned, "unpruned terminal GU and ", nbignored, "ignored.")
    return newmtg, newids

