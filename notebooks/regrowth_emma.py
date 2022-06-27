from numpy.random import binomial, poisson, uniform, normal
import sys
import numpy as np
def binomial_proba(intercept, slope, factor):
    linear = intercept +slope*factor
    return exp(linear)/(1+exp(linear))

def binomial_realization(proba):
    return bool( binomial(1,proba) )

def poisson_proba(intercept, slope, factor):
    return exp(intercept+slope*factor)

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

def lateral_directions(maindir, angle, nb):
    from math import pi
    from openalea.plantgl.all import direction, Matrix3
    assert nb >= 1
    rotdir = direction(maindir.anOrthogonalVector())
    rotmat = Matrix3.axisRotation(rotdir, radians(angle))
    v0 = rotmat * maindir
    deltahangle = 2*pi/nb
    result = [v0] + [Matrix3.axisRotation(maindir, i * deltahangle) *v0 for i in range(1,nb)]
    return result

from mangoG3 import *

def create_daughters(mtg, vid, apical, nblateral, burstdate, totalleafarea = None,  individualleafarea = None):
    parentdirection = get_gu_normed_direction(mtg, vid)
    topposition = get_gu_top_position(mtg, vid)
    diam = get_gu_diameter(mtg, vid)
    newgus = []
    # newconnections = []
    lengths = []
    if not totalleafarea is None:
        set_gu_property(mtg, vid, "RegeneratedLeafArea", totalleafarea)

    if apical:
        apicaldaughter = mtg.add_child(vid, edge_type = '<', UnitType = 'U', label = 'S'+str(int(mtg.label(vid)[1:])+1))        
        l = gu_length(eApical, gu_position(mtg, vid) )
        set_gu_top_position(mtg, apicaldaughter, topposition + parentdirection * l )
        set_gu_diameter(mtg, apicaldaughter, diam)
        newgus.append(apicaldaughter)
        lengths.append(l)
        set_gu_property(mtg, apicaldaughter, "Regrowth", True)
        set_gu_property(mtg, apicaldaughter, "BurstDate", burstdate)
        
    branching_angle = 60
    if nblateral > 0:
        for latiter, latdirection in enumerate(lateral_directions(parentdirection, branching_angle, nblateral)):
            lateralconnection = mtg.add_child(vid, edge_type = '+', label = 'S1', Position = topposition, UnitType = 'U')
            lateraldaughter   = mtg.add_child(lateralconnection, edge_type = '<', label = 'S2', UnitType = 'U', Diameter=diam/10.)
            l = gu_length(eLateral, gu_position(mtg, vid) )
            set_gu_top_position(mtg, lateraldaughter, topposition + latdirection * l )
            set_gu_diameter(mtg, lateraldaughter, diam)
            # newconnections.append(lateralconnection)
            newgus.append(lateraldaughter)
            lengths.append(l)
            set_gu_property(mtg, lateraldaughter, "Regrowth", True)
            set_gu_property(mtg, lateraldaughter, "BurstDate", burstdate)
            #print(get_gu_top_position(mtg, lateraldaughter),get_gu_bottom_position(mtg, lateraldaughter), get_gu_diameter(mtg, lateraldaughter))

    totlength = sum(lengths)
    if not totalleafarea is None and totalleafarea > 0:
        unitla = totalleafarea / totlength
        for vid, l in zip(newgus, lengths):
            leafarea = unitla*l
            set_gu_property(mtg, vid, "LeafArea", leafarea)
            if not individualleafarea is None:
                set_gu_property(mtg, vid, "NbLeaf", int(round(leafarea / individualleafarea)))
    
    return newgus
        

def gu_length(position, motherposition):
    gu_length_distrib = { (eApical, eApical)   : ( 18.14 , 4.14 ) ,
                          (eApical, eLateral)  : ( 13.79 , 4.03 ) ,
                          (eLateral, eApical)  : ( 12.59 , 3.38 ) ,
                          (eLateral, eLateral) : ( 12.59 , 3.38 ) }
    mean, sd = gu_length_distrib[(position, motherposition)]
    return normal_realization(mean, sd, 25, 5)

def gu_nb_leaf(position):
    leaf_nb_distrib = { eApical  : ( 0.59, 5.5), eLateral : ( 0.62, 0.36) }
    return normal_realization(*leaf_nb_distrib[position])
   
def burst_pruned_2017(intensity, diameter):
    probas = { T1 : (-2.479643,0.29399), 
               T3 : (-1.284197,0.29399) }
    intercept, slope = probas[intensity]
    probavalue = binomial_proba(intercept,slope,diameter)
    return binomial_realization( probavalue)

 ################################################################################# Ajout 2022 Occurence du débourrement

def burst_pruned_2022(diameter, intensity, TrPPFD_day_min, zeta_day_12):
    intercept = -8.278 
    linear = intercept + 0.391 * diameter + 8.357 * intensity + 20.137 * TrPPFD_day_min + 1.122 * zeta_day_12
    probavalue = exp(linear)/(1+exp(linear))

    return binomial_realization( probavalue)

def burst_unpruned_2022(diameter, intensity, zeta_day_min):
    intercept = -5.996 
    linear = intercept + 0.315 * diameter + 6.074 * intensity + 3.729 * zeta_day_min
    probavalue = exp(linear)/(1+exp(linear))
        
    return binomial_realization( probavalue)

################################################################################# Ajout 2022 Intensité du débourrement

# interaction_1 = diameter:TrPPFD_day_min 
# interaction_2 = intensite:TrPPFD_day_min

def nb_daughter_pruned_2022(diameter, intensity, TrPPFD_day_min, interaction_1, interaction_2):
   intercept = - 1.840
   probavalue = exp(intercept + 0.163 * diameter + 2.459 * intensity + 19.018 * TrPPFD_day_min + (-0.834 * interaction_1) + (-13.037 * interaction_2))
    try:
        return poisson_realization( probavalue, 10)+1  #### A remplacer par rbinom ?
    except ValueError as ve:
        print(diameter)
        raise ve

def nb_daughter_unpruned_2022(diameter, intensity, zeta_day_min):
    intercept = -4.111
    probavalue = exp(intercept + 0.237 * diameter + 2.231 * intensity + 2.967 * zeta_day_min)
    try:
        return poisson_realization(probavalue, 10)+1 ## remplacer par rbinom ? 
    except ValueError as ve:
        print(intensity, diameter, zeta_day_min)
        raise ve

################################################################################# Ajout 2022 Dynamique du débourrement
 def flush_selection_2022(intensity):
    intercept = 1.409
        linear = intercept + 6.232 * intensity
        firstflushproba = exp(linear)/(1+exp(linear))       
    return binomial_realization(firstflushproba)

def burstdelay_2022(intensity):
    
    firstflush = flush_selection_2022(intensity)
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

def total_leafarea_pruned(intensity, diameter):
    """ Total leaf area generated from a pruned gu in dm2 """
    probas = (0,0.6351646,4.978825)
    intercept, slope,sd = probas
    probavalue = slope*diameter+intercept
    minval=0.1 
    maxval=28
    return  normal_realization(probavalue, sd, maxval, minval)

def individual_leafarea_pruned(intensity):
    """ Individual leaf area of GU borned from a pruned GU in dm2 """
    probas = { T1 : (47.1,25.9), 
               T3 : (62.3,21.1) }
    probavalue, sd = probas[intensity]
    minval=6 
    maxval=150
    return  normal_realization(probavalue, sd, maxval, minval)/100.

 def total_leafarea_unpruned(intensity, diameter, apical_bud):
    """ Total leaf area generated from a pruned gu in dm2 """
    if apical_bud :
        probas = { T0 : (-0.3099147, 0.6698266,2.205946),
                   T1 : (0.863612 ,0.6698266,2.342514), 
                   T3 : (2.579996, 0.6698266,2.458084) }
    else:
        probas = { T0 : (-1.913821, 0.6698266,2.361798),
                   T1 : (-0.740294 ,0.6698266,2.342514), 
                   T3 : (2.579996, 0.6698266,1.290018) }
        
    intercept, slope,sd = probas[intensity]
    probavalue = slope*diameter+intercept
    minval=0.1
    maxval=20
    return  normal_realization(probavalue, sd, maxval, minval)

def individual_leafarea_unpruned(intensity):
    """ Individual leaf area of GU borned from an unpruned GU in dm2 """
    probas = { T0 : (37.9,18.7), 
               T1 : (49.6,26.4), 
               T3 : (70.4,29.03) }
    probavalue, sd = probas[intensity]
    minval=6
    maxval=150
    return  normal_realization(probavalue, sd, maxval, minval)/100.

def burstdelay_pruned_2017(intensity,  severity):
    minval = 10
    
    if intensity == T1 :
        
        flushproba = { 1  : 0.28, 2 : 0.69, 3 : 0.64}
        firstflush = binomial_realization(flushproba[severity])
        
        if firstflush : # Correspond to first flush
            probas = { 1 : (14.7, 1.5) ,
                       2 : (19.4, 2.9), 
                       3 : (17.8, 2.1) }
        else:
            probas = { 1 : (28.2, 3.9),
                       2 : (28.3, 0.7), 
                       3 : (28.3, 2.0) }
        mean, sd = probas[severity]
        maxval = 40
        
    if intensity== T3 :
        mean, sd = 16.5, 2.6
        maxval = 30
    
    return normal_realization(mean, sd, maxval, minval)

   
def flush_selection_unpruned_2017(intensity):
    if intensity == T0 :
        firstflushproba = 0.85
    elif intensity == T1 :
        firstflushproba = 0.48        
    return binomial_realization(firstflushproba)

def burstdelay_unpruned_2017(intensity):
    minval = 10
    maxval = 35
    
    if intensity == T0 :
        firstflush = flush_selection_unpruned(intensity)
        if firstflush :
            mean, sd = 13.4, 3.0
        else:
            mean, sd = 27.5, 1.1
    if intensity == T1:
        firstflush = flush_selection_unpruned(intensity)
        if firstflush :
            mean, sd = 13.5, 2.0
        else:
            mean, sd = 27.5, 1.1
        
    if intensity == T3:
        mean, sd = 16.5, 2.6
        maxval = 25
        
    return  normal_realization(mean, sd, maxval, minval)


from datetime import date, timedelta

def growth_pruned_gu_2017(mtg, vid, intensity, severity, pruningdate):
    diameter = get_gu_diameter(mtg, vid)
    veggrowth = burst_pruned(intensity, diameter)
    if veggrowth:
        nbdaughter = nb_daughter_pruned(intensity, diameter)+1
        assert nbdaughter >= 1
        totalleafarea = total_leafarea_pruned(intensity, diameter)
        individualleafarea = individual_leafarea_pruned(intensity)
        burstdate = pruningdate + timedelta(days = burstdelay_pruned(intensity,  severity))
        return create_daughters(mtg, vid, 0, nbdaughter, burstdate, totalleafarea, individualleafarea)       


def growth_unpruned_gu_2017(mtg, vid, intensity, pruningdate):
    diameter = get_gu_diameter(mtg, vid)
    apical_bud = binomial_realization(0.52)
    veggrowth = burst_unpruned(intensity, diameter, apical_bud)
    if veggrowth:
        nbdaughter = nb_daughter_unpruned(intensity, diameter, apical_bud)
        assert nbdaughter >= 1
        totalleafarea = total_leafarea_unpruned(intensity, diameter, apical_bud)
        individualleafarea = individual_leafarea_unpruned(intensity)
        burstdate = pruningdate + timedelta(days = burstdelay_unpruned(intensity))
        return create_daughters(mtg, vid, int(apical_bud), nbdaughter-int(apical_bud), burstdate, totalleafarea, individualleafarea)

import util ; reload(util)
from util import *


def growth(mtg, listidpruned, intensity, pruningdate = date(2017,2,1), maxdiamunpruned = 10):
    newmtg = deepcopy(mtg)
    listidpruned = dict([(ni, set(ids)) for ni,ids in listidpruned.items() ])
    terminals = get_all_terminal_gus(mtg)
    nbterminals = len(terminals)
    newids = []
    print("Should examine", nbterminals, "terminal GUs.")
    nbpruned, nbunpruned, nbignored = 0,0,0
    for count, vid in enumerate(terminals):
        pruned = False
        lnewids = None
        for severity, ids in listidpruned.items():
            if vid in ids:
                # we consider a pruned gu:
                lnewids = growth_pruned_gu(newmtg, vid, intensity, severity, pruningdate)
                nbpruned += 1
                pruned = True
                break
        if not pruned and get_gu_diameter(mtg, vid) <= maxdiamunpruned and not 'A' in mtg.property('Taille').get(vid,''):
            nbunpruned += 1
            lnewids = growth_unpruned_gu(newmtg, vid, intensity, pruningdate)
        elif not pruned:
            nbignored += 1
        if lnewids:
            newids += lnewids
    print("Processed", nbpruned, "pruned terminal GU and", nbunpruned, "unpruned terminal GU and ", nbignored, "ignored.")
    return newmtg, newids


class GrowthColoring: 
    def __init__(self, newids):
        self.newids = set(newids)
    def prepare_turtle(self, turtle):
        from openalea.plantgl.all import Material
        turtle.setMaterial(1,Material((45,65,15))) # ,transparency=0.8))
        turtle.setMaterial(10,Material((200,200,0)))
        turtle.setMaterial(11,Material((200,0,0)))
    def set_mtg(self, mtg):
        self.mtg = mtg            
        self.colors = { False : 1, True : 2}
        self.mindate = min(self.mtg.property('BurstDate').values()) #date(2017,2,1)
        self.maxdate = max(self.mtg.property('BurstDate').values()) #date(2017,6,1)
        print(self.mindate, self.maxdate)
        self.deltadate = float((self.maxdate - self.mindate).days)
    def __call__(self, turtle, vid):
        if vid in self.newids:
            d = (get_gu_property(self.mtg, vid, "BurstDate")-self.mindate).days
            turtle.interpolateColors(10,11,d/self.deltadate)
        else:
            turtle.setColor(1)


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
def growth_pruned_gu_2022(mtg, vid, intensity, TrPPFD_day_min, zeta_day_12, pruningdate):
    # interaction_1 = diameter:TrPPFD_day_min 
    # interaction_2 = intensite:TrPPFD_day_min
    diameter = get_gu_diameter(mtg, vid)
    veggrowth = burst_pruned_2022(diameter, intensity, TrPPFD_day_min, zeta_day_12)
    if veggrowth:
        nbdaughter = nb_daughter_pruned_2022(diameter, intensity, TrPPFD_day_min, interaction_1, interaction_2) # Pourquoi  +1 ?
        assert nbdaughter >= 1
        totalleafarea = total_leafarea_pruned(intensity, diameter)
        individualleafarea = individual_leafarea_pruned(intensity)
        burstdate = pruningdate + timedelta(days = burstdelay_2022(intensity))
        return create_daughters(mtg, vid, 0, nbdaughter, burstdate, totalleafarea, individualleafarea)       


def growth_unpruned_gu_2022(mtg, vid, intensity, zeta_day_min, pruningdate):
    diameter = get_gu_diameter(mtg, vid)
    apical_bud = binomial_realization(0.52) # On le laisse ?
    veggrowth = burst_unpruned_2022(diameter, intensity, zeta_day_min)
    if veggrowth:
        nbdaughter = nb_daughter_unpruned_2022(diameter, intensity, zeta_day_min)
        assert nbdaughter >= 1
        totalleafarea = total_leafarea_unpruned(intensity, diameter, apical_bud)
        individualleafarea = individual_leafarea_unpruned(intensity)
        burstdate = pruningdate + timedelta(days = burstdelay_2022(intensity))
        return create_daughters(mtg, vid, int(apical_bud), nbdaughter-int(apical_bud), burstdate, totalleafarea, individualleafarea)

################################################################################# Ajout 2022 Mortalité des UCs mères
def Kill_pruned_GU_mother(mtg, vid, intensity, severity, pruningdate):
    diameter = get_gu_diameter(mtg, vid)
    # Mettre une condition, car il n'y a que sur les UCs taillées qu'il faut tuer des UCs !
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

def Kill_leaf(mtg, vid, rank):
    # Ecrire une fonction qui définie le rang
    Leaf_fall = nb_leaf_fall_2022(rank)
    if Leaf_fall:
        prop_leaf_fall = prop_leaf_fall(intensity, diameter, apical_bud)
        # Soustraire la proportion d'UC qui chute
        return # retourner la maquette sans les UCs mortes pour lancer les lois de débourrement

## Etape finale, combiner tout les processus ? 

def growth_2022(mtg, listidpruned, intensity, pruningdate = date(2021,2,1), maxdiamunpruned = 10):
    newmtg = deepcopy(mtg)
    listidpruned = dict([(ni, set(ids)) for ni,ids in listidpruned.items() ])
    terminals = get_all_terminal_gus(mtg)
    nbterminals = len(terminals)
    newids = []
    print("Should examine", nbterminals, "terminal GUs.")
    nbpruned, nbunpruned, nbignored = 0,0,0
    for count, vid in enumerate(terminals):
        pruned = False
        lnewids = None
        for severity, ids in listidpruned.items():
            if vid in ids:
                # we consider a pruned gu:
                lnewids = growth_pruned_gu(newmtg, vid, intensity, severity, pruningdate)
                nbpruned += 1
                pruned = True
                break
        if not pruned and get_gu_diameter(mtg, vid) <= maxdiamunpruned and not 'A' in mtg.property('Taille').get(vid,''):
            nbunpruned += 1
            lnewids = growth_unpruned_gu(newmtg, vid, intensity, pruningdate)
        elif not pruned:
            nbignored += 1
        if lnewids:
            newids += lnewids
    print("Processed", nbpruned, "pruned terminal GU and", nbunpruned, "unpruned terminal GU and ", nbignored, "ignored.")
    return newmtg, newids