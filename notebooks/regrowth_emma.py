import numpy as np
from regrowth_base import *
import mangoG3 as mg
from pruning import intensity_level

def burst_pruned(diameter, intensity, TrPPFD_day_min, zeta_day_12):
    intercept = -8.278 
    latent = intercept + 0.391 * diameter + 8.357 * intensity + 20.137 * TrPPFD_day_min + 1.122 * zeta_day_12
    probavalue = binomial_proba_from_latent(latent)

    return binomial_realization( probavalue)

def burst_unpruned(diameter, intensity, zeta_day_min):
    if np.isnan(zeta_day_min):
        zeta_day_min = 0
    intercept = -5.996 
    latent = intercept + 0.315 * diameter + 6.074 * intensity + 3.729 * zeta_day_min

    probavalue = binomial_proba_from_latent(latent)
    return binomial_realization(probavalue)

################################################################################# Ajout 2022 Intensité du débourrement

def nb_daughter_pruned(diameter, intensity, TrPPFD_day_min):
    intercept = - 1.840
    probavalue = poisson_proba_from_latent(intercept + 0.163 * diameter + 2.459 * intensity + 19.018 * TrPPFD_day_min + (-0.834 * diameter*TrPPFD_day_min) + (-13.037 * intensity*TrPPFD_day_min))
    try:
        return poisson_realization( probavalue, 10)+1  #### A remplacer par rbinom ?
    except ValueError as ve:
        print(diameter)
        raise ve

def nb_daughter_unpruned(diameter, intensity, zeta_day_min):
    intercept = -4.111
    probavalue = poisson_proba_from_latent(intercept + 0.237 * diameter + 2.231 * intensity + 2.967 * zeta_day_min)
    try:
        return poisson_realization(probavalue, 10)+1 ## remplacer par rbinom ? 
    except ValueError as ve:
        print(intensity, diameter, zeta_day_min)
        raise ve

################################################################################# Ajout 2022 Dynamique du débourrement
def flush_selection(intensity):
    intercept = 1.409
    latent = intercept + 6.232 * intensity
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
        sd = 0.75
    return  normal_realization(mean, sd, maxval, minval)
 #################################################################################


from datetime import date, timedelta

def growth_pruned_gu(mtg, vid, intensity, TrPPFD_day_min, zeta_day_12, pruningdate):
    diameter = mg.get_gu_diameter(mtg, vid)
    veggrowth = burst_pruned(diameter, intensity, TrPPFD_day_min, zeta_day_12)
    if veggrowth:

        nbdaughter = nb_daughter_pruned(diameter, intensity, TrPPFD_day_min) # Pourquoi  +1 ?
        assert nbdaughter >= 1
        ilevel = intensity_level(intensity)
        totalleafarea = total_leafarea_pruned(ilevel, diameter)
        individualleafarea = individual_leafarea_pruned(ilevel)
        burstdate = pruningdate + timedelta(days = burstdelay(intensity))
        return create_daughters(mtg, vid, 0, nbdaughter, burstdate, totalleafarea, individualleafarea)       


def growth_unpruned_gu(mtg, vid, intensity, zeta_day_min, pruningdate):
    if np.isnan(zeta_day_min):
        zeta_day_min = 0
    diameter = mg.get_gu_diameter(mtg, vid)
    apical_bud = binomial_realization(0.52) # On le laisse ?
    veggrowth = burst_unpruned(diameter, intensity, zeta_day_min)
    if veggrowth:
        nbdaughter = nb_daughter_unpruned(diameter, intensity, zeta_day_min)
        assert nbdaughter >= 1
        ilevel = intensity_level(intensity)
        totalleafarea = total_leafarea_unpruned(ilevel, diameter, apical_bud)
        individualleafarea = individual_leafarea_unpruned(ilevel)
        burstdate = pruningdate + timedelta(days = burstdelay(intensity))
        return create_daughters(mtg, vid, int(apical_bud), nbdaughter-int(apical_bud), burstdate, totalleafarea, individualleafarea)



def growth(mtg, TrPPFD_day_min = None, zeta_day_12 = None, zeta_day_min = None, intensity = None, 
           pruningdate = date(2021,2,1), maxdiamunpruned = 10):
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
            lnewids = growth_pruned_gu(newmtg, vid, intensity, TrPPFD_day_min.get(vid,0), zeta_day_12.get(vid,0), pruningdate)
            nbpruned += 1
        elif mg.get_gu_diameter(mtg, vid) <= maxdiamunpruned and not 'A' in mtg.property('Taille').get(vid,''):
            nbunpruned += 1
            lnewids = growth_unpruned_gu(newmtg, vid, intensity, zeta_day_min.get(vid,0), pruningdate)
        else :
            nbignored += 1
        if lnewids:
            newids += lnewids
    print("Processed", nbpruned, "pruned terminal GU and", nbunpruned, "unpruned terminal GU and ", nbignored, "ignored.")
    return newmtg, newids


################################################################################# Ajout 2022 Mortalité des UCs mères
def Death_pruned_2022(diameter, intensity, severity, TrPPFD_day_min):
    intercept = 0.456 
    Coef_severity = { n1 : -0.703, 
                      n2 : 0.244,
                      n3 : 0.456 }
    linear = intercept + (-0.237) * diameter + 1.989 * intensity + (-16.27) * TrPPFD_day_min + Coef_severity[severity] 
    probavalue = exp(linear)/(1+exp(linear))

    return binomial_realization(probavalue)

### Il n'y a pas de facteurs significatifs sur la mortalité des UCs mères non taillées (unprunned)
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
def nb_daughter_deathGUf_pruned_2022(ToT_GUf, PPFDcum_day_max):
    intercept = -0.649
    probavalue = exp(intercept + 0.198 * ToT_GUf + (-0.142) * PPFDcum_day_max)
    try:
        return poisson_realization(probavalue, 10)+1 ## remplacer par rbinom ? 
    except ValueError as ve:
        print(ToT_GUf, PPFDcum_day_max)
        raise ve
## Pas assez d'individus pour voir l'effet de facteurs sur l'intensite de la mortalité
## des GUfilles des UCs non taillées

################################################################################# Ajout 2022 Croissance 

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

