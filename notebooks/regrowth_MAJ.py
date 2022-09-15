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

################################################################################# Ajout 2022 Mortalité des UCs mères
def Death_pruned_2022(diameter, intensity, severity, Zeta_min):
    if np.isnan(Zeta_min):
        Zeta_min = 0
    intercept = 0.888 
    Coef_severity = { n1 : 0.091, 
                      n2 : 0.198,
                      n3 : 0.232 }
    linear = intercept + (-0.227) * diameter + 2.324 * intensity + (-2.697) * Zeta_min + Coef_severity[severity] 
    probavalue = exp(linear)/(1+exp(linear))

    return binomial_realization(probavalue)

def Death_unpruned_2022(diameter, intensity, severity, Zeta_min):
    intercept = -1.597 
    linear = intercept + (-3.696) * Zeta_min 
    probavalue = exp(linear)/(1+exp(linear))

    return binomial_realization(probavalue)
################################################################################# Ajout 2022 Mortalité des UCs filles
## Occurence
def DeathOccurency_GUf_unpruned_2022(ToT_GUf):
    intercept = - 7.641
    linear = intercept + 1.503 * ToT_GUf
    probavalue = exp(linear)/(1+exp(linear))

    return binomial_realization(probavalue)

def DeathOccurency_GUf_pruned_2022(ToT_GUf):
    intercept = - 2.362
    linear = intercept + 0.473 * ToT_GUf
    probavalue = exp(linear)/(1+exp(linear))

    return binomial_realization(probavalue)

## Intensité
def nb_daughter_deathGUf_pruned_2022(ToT_GUf):
    intercept = -1.514
    probavalue = exp(intercept + 0.205 * ToT_GUf)
    try:
        return poisson_realization(probavalue, 10)+1 ## remplacer par rbinom ? 
    except ValueError as ve:
        print(ToT_GUf, PPFDcum_day_max)
        raise ve
## Pas assez d'individus pour voir l'effet de facteurs sur l'intensite de la mortalité
## des GUfilles des UCs non taillées

################################################################################# 

################################################################################# Ajout 2022 Mortalité des UCs mères
def Kill_pruned_GU_mother(mtg, vid, intensity, severity, pruningdate):
    diameter = get_gu_diameter(mtg, vid)
    # Metter une condition, car il n'y a que sur les UCs taillées qu'il faut tuer des UCs !
    Death = Death_pruned_2022(diameter, intensity, severity, TrPPFD_day_min)
    if Death:
        ## Ajouter le code pour faire disparaitre l'UC de la maquette 
        return # retourner la maquette sans les UCs mortes pour lancer les lois de débourrement

################################################################################# Ajout 2022 Mortalité des UCs filles
def Kill_unpruned_GUf_2022(mtg, ToT_GUf):
    # 1) Sélection des UCfilles
    
    Death_GUf = DeathOccurency_GUf_unpruned_2022(ToT_GUf) # 2) Application de l'occurence de leur mortalité
    ## Prédiction du nombre de filles qui meurt
    if Death_GUf :
        nbdaughter = 1
        ## + Code pour faire disparaitre l'UC fille
        return        

def Kill_pruned_GUf_2022(mtg, ToT_GUf, PPFDcum_day_max):
    # 1) Sélection des UCfilles
    
    Death_GUf = DeathOccurency_GUf_pruned_2022(ToT_GUf) # 2) Application de l'occurence de leur mortalité
    ## Prédiction du nombre de filles qui meurt
    if Death_GUf :
        nbdaughter = nb_daughter_deathGUf_pruned_2022(ToT_GUf, PPFDcum_day_max)
        ## + Code pour faire disparaitre les UC filles
        return    

################################################################################# Ajout 2022 Chute des UCs feuilles
def leaf_fall_2022(rank):
    Coef_rank = { 1 : 0.06, 
                  2 : 0.23,
                  4 : 0.33}
    probavalue = Coef_rank[rank] 

    return binomial_realization(probavalue)

def prop_leaf_fall_2022(rank):
    Coef_rank = { 1 : 0.225, 
                 2 : 0.180,
                 4 : 0.142}
    probavalue = Coef_rank[rank] 

    return binomial_realization(probavalue) 

## N.B Pour l'intensité de la chute des feuille, on a travaillé dans le modèle avec une loi binomiale, puisque la variable réponse "y" est une matrice entre
# y1 = (nb feuille mois n-1) - (nb feuille mois n) : succès
# y2 = nb noeuds - {(nb feuille mois n-1) - (nb feuille mois n)} : Echec

def leaf_death(mtg, vid):
    # Ecrire une fonction qui définie le rang
    Leaf_fall = nb_leaf_fall_2022(rank)
    if Leaf_fall:
        prop_leaf_fall = prop_leaf_fall(intensity, diameter, apical_bud)
        # Soustraire la proportion d'UC qui chute
        return # retourner la maquette sans les UCs mortes pour lancer les lois de débourrement

## Etape finale, combiner tout les processus ? 