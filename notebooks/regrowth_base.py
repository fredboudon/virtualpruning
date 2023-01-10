from numpy.random import binomial, poisson, uniform, normal
import sys
from math import exp, radians

from pruning import T0, T1, T2, T3
from mangoG3 import eApical, eLateral


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


def create_daughters(mtg, vid, apical, nblateral, burstdate, totalleafarea = None,  individualleafarea = None):
    import mangoG3 as mg3
    parentdirection = mg3.get_gu_normed_direction(mtg, vid)
    topposition = mg3.get_gu_top_position(mtg, vid)
    diam = mg3.get_gu_diameter(mtg, vid)
    newgus = []
    # newconnections = []
    lengths = []
    if not totalleafarea is None:
        mg3.set_gu_property(mtg, vid, "RegeneratedLeafArea", totalleafarea)

    if apical:
        apicaldaughter = mtg.add_child(vid, edge_type = '<', UnitType = 'U', label = 'S'+str(int(mtg.label(vid)[1:])+1))        
        l = gu_length(mg3.eApical, mg3.gu_position(mtg, vid) )
        mg3.set_gu_top_position(mtg, apicaldaughter, topposition + parentdirection * l )
        mg3.set_gu_diameter(mtg, apicaldaughter, diam)
        newgus.append(apicaldaughter)
        lengths.append(l)
        mg3.set_gu_property(mtg, apicaldaughter, "Regrowth", True)
        mg3.set_gu_property(mtg, apicaldaughter, "BurstDate", burstdate)
        
    branching_angle = 60
    if nblateral > 0:
        for latiter, latdirection in enumerate(lateral_directions(parentdirection, branching_angle, nblateral)):
            lateralconnection = mtg.add_child(vid, edge_type = '+', label = 'S1', Position = topposition, UnitType = 'U')
            lateraldaughter   = mtg.add_child(lateralconnection, edge_type = '<', label = 'S2', UnitType = 'U', Diameter=diam/10.)
            l = gu_length(mg3.eLateral, mg3.gu_position(mtg, vid) )
            mg3.set_gu_top_position(mtg, lateraldaughter, topposition + latdirection * l )
            mg3.set_gu_diameter(mtg, lateraldaughter, diam)
            # newconnections.append(lateralconnection)
            newgus.append(lateraldaughter)
            lengths.append(l)
            mg3.set_gu_property(mtg, lateraldaughter, "Regrowth", True)
            mg3.set_gu_property(mtg, lateraldaughter, "BurstDate", burstdate)
            #print(get_gu_top_position(mtg, lateraldaughter),get_gu_bottom_position(mtg, lateraldaughter), get_gu_diameter(mtg, lateraldaughter))

    totlength = sum(lengths)
    if not totalleafarea is None and totalleafarea > 0:
        unitla = totalleafarea / totlength
        for vid, l in zip(newgus, lengths):
            leafarea = unitla*l
            mg3.set_gu_property(mtg, vid, "LeafArea", leafarea)
            if not individualleafarea is None:
                mg3.set_gu_property(mtg, vid, "NbLeaf", int(round(leafarea / individualleafarea)))
    
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


class GrowthColoring: 
    def __init__(self):
        pass
    def prepare_turtle(self, turtle):
        from openalea.plantgl.all import Material
        turtle.setMaterial(1,Material((45,65,15))) # ,transparency=0.8))
        turtle.setMaterial(10,Material((200,200,0)))
        turtle.setMaterial(11,Material((200,0,0)))
    def set_mtg(self, mtg):
        self.mtg = mtg
        self.newids = set([vid for vid,regrowth in mtg.property('Regrowth').items() if regrowth == True])
        self.colors = { False : 1, True : 2}
        self.mindate = min(self.mtg.property('BurstDate').values()) #date(2017,2,1)
        self.maxdate = max(self.mtg.property('BurstDate').values()) #date(2017,6,1)
        #print(self.mindate, self.maxdate)
        self.deltadate = float((self.maxdate - self.mindate).days)
    def __call__(self, turtle, vid):
        if vid in self.newids:
            d = (self.mtg.property('BurstDate')[vid]-self.mindate).days
            turtle.interpolateColors(10,11,d/self.deltadate)
        else:
            turtle.setColor(1)



def plot_growth(mtg, **options):
    import mtgplot as mp
    options.setdefault('leaves',True)
    options.setdefault('gc',False)
    return mp.plot_tree(mtg, colorizer = GrowthColoring, **options)


def plot_growth_dynamic(regrowth, **options):
    from ipywidgets import interact, interactive, fixed, interact_manual
    import ipywidgets as widgets
    import datetime
    from IPython.display import display
    import mtgplot as mp
    options.setdefault('leaves',True)
    options.setdefault('gc',False)
    mindate = min(regrowth.property('BurstDate').values())
    maxdate = max(regrowth.property('BurstDate').values())
    print(mindate,'---',maxdate)
    nbdays = (maxdate-mindate).days
    sw = mp.plot_tree(regrowth, colorizer = GrowthColoring, todate = maxdate, **options)
    display(sw)
    x=widgets.IntSlider(min=0, max=nbdays, step=1, value=nbdays)
    def plot_dyn(x):
        print(mindate+datetime.timedelta(days=x))
        sc = mp.representation(regrowth, colorizer = GrowthColoring,  todate = mindate+datetime.timedelta(days=x), **options)
        sw.set_scenes([sc]) 
    im = interact_manual(plot_dyn, x=x)
    im.widget.children[0].description = 'Display'
    return im