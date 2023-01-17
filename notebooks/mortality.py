from randomgeneration import *
import mangoG3 as mg
from pruning import n1, n2, n3

################################################################################# Ajout 2022 Mortalité des UCs mères
def pruned_gu_mortality_1(diameter, intensity, severity, Zeta_mean):
    if np.isnan(Zeta_mean):
        Zeta_mean = 0
    intercept = 1.40 
    Coef_severity = { n1 : 0.091, 
                      n2 : 0.198,
                      n3 : 0.232 }
    linear = intercept + (-0.25) * diameter + 2.29 * intensity + (-1.85) * Zeta_mean + Coef_severity[severity] 
    probavalue = binomial_proba_from_latent(linear)

    return binomial_realization(probavalue)

def pruned_gu_mortality_2(diameter, intensity, severity, zone):
    intercept = 0.021 
    Coef_severity = { n1 : 0.09, 
                      n2 : 0.20,
                      n3 : 0.24 }
    Coef_zone = { I : 0.11, 
                  B : 0.25,
                  H : 0.17 }
    linear = intercept + (-0.24) * diameter + 1.63 * intensity + Coef_zone[zone] + Coef_severity[severity] 
    probavalue = binomial_proba_from_latent(linear)

    return binomial_realization(probavalue)

def pruned_gu_mortality_3():
    intercept = -1.32 
    probavalue = binomial_proba_from_latent(intercept)

    return binomial_realization(probavalue)

def unpruned_gu_mortality_1(Zeta_mean):
    if np.isnan(Zeta_mean):
        Zeta_mean = 0
    intercept = -1.67 
    linear = intercept + (-1.71) * Zeta_mean 
    probavalue = binomial_proba_from_latent(linear)

    return binomial_realization(probavalue)

def unpruned_gu_mortality_2():
    intercept = -2.79 
    probavalue = binomial_proba_from_latent(intercept)

    return binomial_realization(probavalue)

################################################################################# Ajout 2022 Mortalité des UCs filles
## Occurence
def unpruned_gu_daughter_mortality_occurence(nb_daughter_gus): 
    intercept = - 7.641
    linear = intercept + 1.503 * nb_daughter_gus
    probavalue = binomial_proba_from_latent(linear)

    return binomial_realization(probavalue)

def unpruned_gu_daughter_mortality_intensity(nb_daughter_gus): 
    if nb_daughter_gus <= 1 : return nb_daughter_gus
    return poisson_realization(poisson_proba_from_latent(-0.64), min(3,nb_daughter_gus-1))+1

def pruned_gu_daughter_mortality_occurence(nb_daughter_gus):
    intercept = - 2.362
    linear = intercept + 0.473 * nb_daughter_gus
    probavalue = binomial_proba_from_latent(linear)

    return binomial_realization(probavalue)


def pruned_gu_daughter_mortality_intensity(nb_daughter_gus):
    if nb_daughter_gus <= 1 : return nb_daughter_gus
    intercept = -1.514
    probavalue = poisson_proba_from_latent(intercept + 0.205 * (nb_daughter_gus-1))
    return poisson_realization(probavalue, nb_daughter_gus-1)+1

################################################################################# 

def remove_gu(mtg, vid):
    toremove = [vid]
    if mtg.edge_type(mtg.parent(vid)) == '+':
        toremove.append(mtg.parent(vid))
    for lvid in toremove:
        mtg.remove_vertex(lvid)
        for propname, propvalues in mtg.properties().items():
            if lvid in propvalues:
                del propvalues[lvid]


def gu_mortalities_post_pruning(mtg, Zeta_mean = None, intensity = None, inplace = False):
   if intensity is None:
       from pruning import continuous_intensity_from_pruned
       intensity = continuous_intensity_from_pruned(mtg)

   if Zeta_mean is None:
        from lightestimation import light_variables_mortality
        from mtgplot import representation
        prunedrepr = representation(mtg, wood = False, leaves=True)
        Zeta_min = light_variables_mortality(prunedrepr)


   if not inplace:
        from copy import deepcopy
        newmtg = deepcopy(mtg)
   else:
        newmtg = mtg
   prunedgus = mtg.property('pruned')

   terminals = mg.get_all_terminal_gus(mtg)
   print("Should examine", len(terminals), " GUs.")

   toremove = []
   for vid in terminals:
        if vid in prunedgus:
            severity =  prunedgus[vid]
            occurence = pruned_gu_mortality(mg.get_gu_diameter(mtg, vid), intensity, severity, Zeta_mean.get(vid,0))
            if occurence:
                toremove.append(vid)
        else:
            occurence = unpruned_gu_mortality(Zeta_mean.get(vid,0))
            if occurence:
                toremove.append(vid)
   
   for vid in toremove:
        assert(newmtg.nb_children(vid) == 0)
        remove_gu(newmtg, vid)
   
   return newmtg, toremove

def gu_mortalities_post_regrowth(mtg, newids, inplace = False):
   if intensity is None:
       from pruning import continuous_intensity_from_pruned
       intensity = continuous_intensity_from_pruned(mtg)

   if not inplace:
        from copy import deepcopy
        newmtg = deepcopy(mtg)
   else:
        newmtg = mtg
   
   prunedgus = mtg.property('pruned')

   toconsider = set([mtg.parent(vid) for vid in newids])

   print("Should examine", len(toconsider), " GUs.")

   toremove = []
   for vid in terminals:
        if vid in prunedgus:
            occurence = pruned_gu_daughter_mortality_occurence(mtg.nb_children(vid))
            if occurence:
                daughters = mg.get_children(mtg,vid)
                nbdeaddaughters = pruned_gu_daughter_mortality_intensity(mtg.nb_children(vid))
                toremove += np.random.choice(daughters, size=nbdeaddaughters, replace=False)
        else:
            occurence = unpruned_gu_daughter_mortality_occurence(mtg.nb_children(vid))
            if occurence:
                daughters = mg.get_children(mtg,vid)
                nbdeaddaughters = unpruned_gu_daughter_mortality_intensity(mtg.nb_children(vid))
                toremove += np.random.choice(daughters, size=nbdeaddaughters, replace=False)

   for vid in toremove:
        assert(newmtg.nb_children(vid) == 0)
        remove_gu(newmtg, vid)

   return newmtg, toremove


def gu_mortality_post_regrowth(mtg, motherid, nb_daughter_gus):
    """ Return the number of daughter gus after mortality post_regrowth """
    if motherid in mtg.property('pruned'):
        occurence = pruned_gu_daughter_mortality_occurence(nb_daughter_gus)
        if occurence:
            nbdeaddaughters = pruned_gu_daughter_mortality_intensity(nb_daughter_gus)
            return nb_daughter_gus - nbdeaddaughters
    else:
        occurence = unpruned_gu_daughter_mortality_occurence(nb_daughter_gus)
        if occurence:
            nbdeaddaughters = unpruned_gu_daughter_mortality_intensity(nb_daughter_gus)
            return nb_daughter_gus - nbdeaddaughters
    return nb_daughter_gus


################################################################################# Ajout 2022 Chute des UCs feuilles
def leaf_fall(rank):
    rank_coef = { 1 : 0.06, 
                  2 : 0.23,
                  4 : 0.33}
    return binomial_realization(rank_coef[rank])

def prop_leaf_fall(rank):
    Coef_rank = { 1 : 0.225, 
                 2 : 0.180,
                 4 : 0.142}
    probavalue = Coef_rank[rank] 

    return binomial_realization(probavalue) 

## N.B Pour l'intensité de la chute des feuille, on a travaillé dans le modèle avec une loi binomiale, puisque la variable réponse "y" est une matrice entre
# y1 = (nb feuille mois n-1) - (nb feuille mois n) : succès
# y2 = nb noeuds - {(nb feuille mois n-1) - (nb feuille mois n)} : Echec

def leaf_mortality(mtg, vid):
    # Ecrire une fonction qui définie le rang
    Leaf_fall = nb_leaf_fall(rank)
    if Leaf_fall:
        prop_leaf_fall = prop_leaf_fall(intensity, diameter, apical_bud)
        # Soustraire la proportion d'UC qui chute
        return # retourner la maquette sans les UCs mortes pour lancer les lois de débourrement

## Etape finale, combiner tout les processus ? 