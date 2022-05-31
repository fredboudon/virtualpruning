from random import randint
from mangoG3 import *

meandepth = lambda depths : int(round(np.mean(depths)))
maxoccurence = lambda depths :max(depths,key=depths.count)
hasdepth = lambda depths : [i for i in range(1,4) if i in depths] 

diameterrange = { 1 : (6.9, 1.6) , 2 : (9.5, 2.5), 3 : (13.1, 3.4)}
n1, n2, n3 = 1,2,3


def valid_diameter(mtg, vid, depth):
    mean, std = diameterrange[depth]
    return mean - std <= get_gu_diameter(mtg, vid) <= mean + std


def determine_potential_cutpoints(mtg, depthfunc = max, diameterconstraints = False):
    """ Determine the ids of the potential n1, n2, n3 cut points 
        from constraint of depth and also diameters observed  experimentally
        :param mtg: the structure of the tree with characteristic
        :return: 3 lists of ids representing n1, n2 and n3
    """
    ids_ni = { n1 : [] , n2 : [], n3 : []}
    depths = gus_depth_from_terminal(mtg, depthfunc, 1)
    for vid in get_all_gus(mtg):            
        depth = depths[vid]
        if depth in [1,2,3] and (not diameterconstraints or valid_diameter(mtg, vid, depth)):
                ids_ni[depth].append(vid)
            
    return ids_ni


def check_cutpoint_diameter_validity(mtg, potential_cutpoints ):
    for o, nis in enumerate(potential_cutpoints, 1):
        for ni in nis:
            assert valid_diameter(g, ni, o)


def tag_pruning(mtg, listidpruned):
    pruningproperty = {}
    if not type(listidpruned) == dict:
        listidpruned = dict(enumerate(listidpruned,1))
    for order, pruneds in listidpruned.items():
        for p in pruneds:
            for vid in get_descendants_gus_from_ancestor(mtg, p):
                if True:
                    if vid in pruningproperty:
                        if order < pruningproperty[vid]:
                            raise ValueError(vid, order, pruningproperty[vid])
                    else:
                        pruningproperty[vid] = order                    
    return pruningproperty

class PruningColoring: 
    def __init__(self):
        pass
    def set_mtg(self,mtg):
        self.mtg = mtg            
        self.colors = { None : 2, 1 : (0,0,255), 2: (255,255,0), 3: 3}
    def __call__(self, turtle, vid):
        color = self.colors[self.mtg.property('pruning').get(vid)]
        if type(color) == int:
            turtle.setColor(color)
        else:
            from openalea.plantgl.all import Material
            turtle.setCustomAppearance(Material(color))
            
def assign_pruning(mtg, listidpruned, checkvalidity = True):
    mtg.property('pruning').clear()
    mtg.property('pruning').update(tag_pruning(mtg, listidpruned))


def plot_pruning(mtg, listidpruned, leaves = True, checkvalidity = True):
    assign_pruning(mtg, listidpruned, checkvalidity)
    import mtgplot as mp
    sc = mp.representation(mtg, colorizer=PruningColoring, gc = False, leaves=leaves)
    return mp.display(sc)


def define_pruning(mtg, nbcuts, potential_cutpoints):
    """ Determine at which point the pruning should occurs.
        :params nbcuts: int or tuple. Nb of cuts. 
        If a single value is given, it applies this number of cut for n1, n2 et n3. 
        If a tuple of 3 value is given,  it is used for n1, n2 et n3 respectivelly.
    """
    import numpy as np
    if type(nbcuts) == int:
        nbcuts = (nbcuts, nbcuts, nbcuts)
    
    assert len(nbcuts) == 3
    
    if potential_cutpoints is None:
        allnis = determine_potential_cutpoints(mtg)
    else:
        import copy
        allnis = copy.deepcopy(potential_cutpoints)
    listidpruned = { n1 : [] , n2 : [], n3 : []}
    
    removedterminals = set()
    for order,nis in reversed(list(allnis.items())):
        if nbcuts[order-1] > 0:
            if len(nis) < nbcuts[order-1]:
                raise ValueError('Cannot make so many cuts for n%i: %i. Maximum is %i.' %(order, nbcuts[order-1], len(nis)))
            nblcuts = 0
            np.random.shuffle(nis)
            for ni in nis:
                mterminals = set(get_terminal_gus_from_ancestor(mtg, get_ancestor(mtg, ni, 3-order) if order < 3 else ni))
                if len(mterminals & removedterminals) == 0:
                    listidpruned[order].append(ni)
                    removedterminals |= mterminals
                    nblcuts += 1
                    if nblcuts == nbcuts[order-1] : break
            if nblcuts != nbcuts[order-1]:
                raise ValueError('Cannot make %i cuts for n%i. Achieved only %i.' %(nbcuts[order-1], order, nblcuts))
    
    return listidpruned

def tag_selected_nodes(mtg, listidpruned):
    pruningproperty = {}
    for order, pruneds in listidpruned.items():
        for vid in pruneds:
            if vid in pruningproperty and order < pruningproperty[vid]:
                if checkvalidity:
                    raise ValueError(vid, order, pruningproperty[vid])
            else:
                pruningproperty[vid] = order
    return pruningproperty

def apply_pruning(mtg, listprunedids, inplace = False):
    if not inplace:
        from copy import deepcopy
        mtg = deepcopy(mtg)
    mtg.property('cuted').clear()
    mtg.property('cuted').update(tag_selected_nodes(mtg, listprunedids))
    pnames = mtg.property_names()
    for order, prunedids in listprunedids.items():
        for pid in prunedids:
            descendants = list(reversed(mtg.Descendants(pid)))
            for d in descendants:
                if d != pid:
                    for pname in pnames:
                        if d in mtg.property(pname):
                            del mtg.property(pname)[d]
                    mtg.remove_vertex(d)
            if mtg.nb_children(pid) != 0 :
                raise ValueError(pid, mtg.children(pid))
    return mtg

def continuous_intensity(mtg, listprunedgu): 
    from mangoG3 import volume, get_gu_diameter
    from allometry import gu_biomass
    mvolume = volume(mtg)
    diameters = [[get_gu_diameter(mtg,vid) for vid in nis] for nis in listprunedgu.values()]
    biomasses = [sum(list(map(gu_biomass, diam)))/1000 for diam in diameters]
    biomass = sum(biomasses)

    return biomass / mvolume 


# Threshold are expressed in kg.m-3
T1threshold= 0.20
T2threshold=0.40
T3threshold = 0.60

T0, T1, T2, T3 = "T0", "T1", "T2", "T3"

def intensity_level(continuous_intensity): 
    if continuous_intensity == 0: return T0
    elif continuous_intensity <= T1threshold: return T1
    elif continuous_intensity <= T2threshold: return T2
    elif continuous_intensity <= T3threshold: return T3
    raise ValueError(continuous_intensity)
    


