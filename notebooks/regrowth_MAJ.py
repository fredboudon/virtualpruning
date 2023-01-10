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
    if np.isnan(TrPPFD_min):
        TrPPFD_min = 0
    if np.isnan(Zeta_8H):
        Zeta_8H = 0
    intercept = - 1.840
    probavalue = poisson_proba_from_latent(intercept + 0.163 * diameter + 2.459 * intensity + 19.018 * TrPPFD_min + 0.622 * Zeta_8H +  
                                           (-0.753 * diameter * TrPPFD_min) )
    try:
        return poisson_realization( probavalue, 10)+1  #### A remplacer par rbinom ?
    except ValueError as ve:
        print(diameter)
        raise ve

def nb_daughter_unpruned(diameter, intensity, Zeta_mean):
    intercept = -8.151
    probavalue = poisson_proba_from_latent(intercept + 0.239 * diameter + 2.164 * intensity + 6.247 * Zeta_mean)
    try:
        return poisson_realization( probavalue, 10)+1 ## remplacer par rbinom ? 
    except ValueError as ve:
        print(intensity, diameter, Zeta_mean)
        raise ve

################################################################################# Ajout 2022 Dynamique du débourrement
def flush_selection(intensity):
    intercept = -2.44
    latent = intercept + -6.02 * intensity
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

def growth_pruned_gu(mtg, vid, intensity, TrPPFD_mean, TrPPFD_min, Zeta_8H, pruningdate):
    if np.isnan(TrPPFD_mean):
        TrPPFD_mean = 0
    if np.isnan(TrPPFD_min):
        TrPPFD_min = 0
    if np.isnan(Zeta_8H):
        Zeta_8H = 0
    diameter = mg.get_gu_diameter(mtg, vid)
    veggrowth = burst_pruned(diameter, intensity, TrPPFD_mean)
    if veggrowth:

        nbdaughter = nb_daughter_pruned(diameter, intensity, TrPPFD_mean, Zeta_8H) # Pourquoi  +1 ?
        assert nbdaughter >= 1
        ilevel = intensity_level(intensity)
        totalleafarea = total_leafarea_pruned(intensity, diameter)
        individualleafarea = individual_leafarea_pruned(ilevel)
        burstdate = pruningdate + timedelta(days = burstdelay(intensity))
        return create_daughters(mtg, vid, 0, nbdaughter, burstdate, totalleafarea, individualleafarea)       


def growth_unpruned_gu(mtg, vid, intensity, TrPPFD_mean, Zeta_mean, pruningdate):
    if np.isnan(TrPPFD_mean):
        TrPPFD_mean = 0
    if np.isnan(Zeta_mean):
        TrPPFD_mean = 0
    diameter = mg.get_gu_diameter(mtg, vid)
    apical_bud = binomial_realization(0.52) # On le laisse ?
    veggrowth = burst_unpruned(diameter, intensity, TrPPFD_mean)
    if veggrowth:
        nbdaughter = nb_daughter_unpruned(diameter, intensity, Zeta_mean)
        assert nbdaughter >= 1
        ilevel = intensity_level(intensity)
        totalleafarea = total_leafarea_unpruned(intensity)
        individualleafarea = individual_leafarea_unpruned(ilevel)
        burstdate = pruningdate + timedelta(days = burstdelay(intensity))
        return create_daughters(mtg, vid, int(apical_bud), nbdaughter-int(apical_bud), burstdate, totalleafarea, individualleafarea)



def growth(mtg, TrPPFD_mean = None, TrPPFD_min = None, Zeta_mean = None, Zeta_min = None, Zeta_8H = None, intensity = None, 
           pruningdate = date(2021,2,24), maxdiamunpruned = 10):
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
            lnewids = growth_pruned_gu(newmtg, vid, intensity, TrPPFD_mean.get(vid,0), TrPPFD_min.get(vid,0),
            Zeta_8H.get(vid,0), pruningdate)
            nbpruned += 1
        elif mg.get_gu_diameter(mtg, vid) <= maxdiamunpruned and not 'A' in mtg.property('Taille').get(vid,''):
            nbunpruned += 1
            lnewids = growth_unpruned_gu(newmtg, vid, intensity, TrPPFD_mean.get(vid,0), Zeta_mean.get(vid,0), pruningdate)
        else :
            nbignored += 1
        if lnewids:
            newids += lnewids
    print("Processed", nbpruned, "pruned terminal GU and", nbunpruned, "unpruned terminal GU and ", nbignored, "ignored.")
    return newmtg, newids

