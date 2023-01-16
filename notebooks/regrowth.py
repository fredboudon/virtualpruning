import numpy as np
from regrowth_base import *
import mangoG3 as mg
from pruning import intensity_level

def burst_pruned(diameter, intensity, TrPPFD_mean):
    if np.isnan(TrPPFD_mean):
        TrPPFD_mean = 0
    intercept = -8.020 
    latent = intercept + 0.380 * diameter + 8.366 * intensity + 20.137 * TrPPFD_mean
    probavalue = binomial_proba_from_latent(latent)

    return binomial_realization(probavalue)

def burst_unpruned(diameter, intensity, TrPPFD_mean):
    if np.isnan(TrPPFD_mean):
        TrPPFD_mean = 0
    intercept = -5.581 
    latent = intercept + 0.303 * diameter + 6.339 * intensity + 4.335 * TrPPFD_mean

    probavalue = binomial_proba_from_latent(latent)
    return binomial_realization(probavalue)

################################################################################# Ajout 2022 Intensité du débourrement

def nb_daughter_pruned(diameter, intensity, TrPPFD_min, Zeta_8H):
    import numpy as np
    if np.isnan(TrPPFD_min):
        TrPPFD_min = 0
    if np.isnan(Zeta_8H):
        Zeta_8H = 0
    intercept = - 1.51
    probavalue = poisson_proba_from_latent(intercept + 0.14 * diameter 
                                                     + 1.24 * intensity 
                                                     + 10.74 * TrPPFD_min 
                                                     + 0.56 * Zeta_8H 
                                                     -0.65 * diameter * TrPPFD_min)
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

def growth_pruned_gu(mtg, vid, intensity, TrPPFD_mean, TrPPFD_min, Zeta_8H, pruningdate, mortalityenabled = True):
    if np.isnan(TrPPFD_mean):
        TrPPFD_mean = 0
    if np.isnan(TrPPFD_min):
        TrPPFD_min = 0
    if np.isnan(Zeta_8H):
        Zeta_8H = 0
    diameter = mg.get_gu_diameter(mtg, vid)
    veggrowth = burst_pruned(diameter, intensity, TrPPFD_mean)
    if veggrowth:
        nbdaughter = nb_daughter_pruned(diameter, intensity, TrPPFD_min, Zeta_8H)
        assert nbdaughter >= 1
        if mortalityenabled:
            nbdaughter = gu_mortality_post_regrowth(mtg, vid, nbdaughter)
        if nbdaughter > 0:
            ilevel = intensity_level(intensity)
            totalleafarea = total_leafarea_pruned(diameter)
            individualleafarea = individual_leafarea_pruned(ilevel)
            burstdate = pruningdate + timedelta(days = burstdelay(intensity))
            return create_daughters(mtg, vid, 0, nbdaughter, burstdate, totalleafarea, individualleafarea)       


def growth_unpruned_gu(mtg, vid, intensity, TrPPFD_mean, Zeta_mean, pruningdate, mortalityenabled = True):
    if np.isnan(TrPPFD_mean):
        TrPPFD_mean = 0
    if np.isnan(Zeta_mean):
        TrPPFD_mean = 0
    diameter = mg.get_gu_diameter(mtg, vid)
    veggrowth = burst_unpruned(diameter, intensity, TrPPFD_mean)
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



def growth(mtg, TrPPFD_mean = None, TrPPFD_min = None, Zeta_mean = None, Zeta_8H = None, intensity = None, 
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
                                       TrPPFD_mean.get(vid,0), 
                                       TrPPFD_min.get(vid,0),
                                       Zeta_8H.get(vid,0), 
                                       pruningdate,
                                       mortalityenabled=mortalityenabled)
            nbpruned += 1
        elif mg.get_gu_diameter(mtg, vid) <= maxdiamunpruned and not 'A' in mtg.property('Taille').get(vid,''):
            lnewids = growth_unpruned_gu(newmtg, vid, 
                                         intensity, 
                                         TrPPFD_mean.get(vid,0), 
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

