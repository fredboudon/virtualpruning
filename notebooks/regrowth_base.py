from pruning import T0, T1, T2, T3
from mangoG3 import eApical, eLateral
from pruning import n1, n2, n3

from randomgeneration import *

def lateral_directions(maindir, angle, nb):
    from math import pi, radians
    from openalea.plantgl.all import direction, Matrix3
    assert nb >= 1
    rotdir = direction(maindir.anOrthogonalVector())
    rotmat = Matrix3.axisRotation(rotdir, radians(angle))
    v0 = rotmat * maindir
    deltahangle = 2*pi/nb
    result = [v0] + [Matrix3.axisRotation(maindir, i * deltahangle) *v0 for i in range(1,nb)]
    return result

def create_daughters_unpruned(mtg, vid, apical, nblateral, burstdate, intensity):
    import mangoG3 as mg3
    parentdirection = mg3.get_gu_normed_direction(mtg, vid)
    topposition = mg3.get_gu_top_position(mtg, vid)
    diam = mg3.get_gu_diameter(mtg, vid)
    newgus = []
    # newconnections = []
    lengths = []

    if apical:
        apicaldaughter = mtg.add_child(vid, edge_type = '<', UnitType = 'U', label = 'S'+str(int(mtg.label(vid)[1:])+1))        
        l = gu_length_unpruned(mg3.eApical)
        nb_leaf = gu_nb_leaf_unpruned(l, mg3.eApical)    
        leaf_area = gu_leaf_area_unpruned(intensity, mg3.eApical)
        mg3.set_gu_top_position(mtg, apicaldaughter, topposition + parentdirection * l )
        mg3.set_gu_diameter(mtg, apicaldaughter, diam)
        newgus.append(apicaldaughter)
        mg3.set_gu_property(mtg, apicaldaughter, "Regrowth", True)
        mg3.set_gu_property(mtg, apicaldaughter, "BurstDate", burstdate)
        mg3.set_gu_property(mtg, apicaldaughter, "NbLeaf", nb_leaf)
        mg3.set_gu_property(mtg, apicaldaughter, "UnitLeafArea", leaf_area)
        
    branching_angle = 60
    if nblateral > 0:
        for latiter, latdirection in enumerate(lateral_directions(parentdirection, branching_angle, nblateral)):
            lateralconnection = mtg.add_child(vid, edge_type = '+', label = 'S1', Position = topposition, UnitType = 'U')
            lateraldaughter   = mtg.add_child(lateralconnection, edge_type = '<', label = 'S2', UnitType = 'U', Diameter=diam/10.)
            l = gu_length_unpruned(mg3.eLateral)
            nb_leaf = gu_nb_leaf_unpruned(l, mg3.eLateral)    
            leaf_area = gu_leaf_area_unpruned(intensity, mg3.eLateral)
            mg3.set_gu_top_position(mtg, lateraldaughter, topposition + latdirection * l )
            mg3.set_gu_diameter(mtg, lateraldaughter, diam)
            newgus.append(lateraldaughter)
            mg3.set_gu_property(mtg, lateraldaughter, "Regrowth", True)
            mg3.set_gu_property(mtg, lateraldaughter, "BurstDate", burstdate)
            mg3.set_gu_property(mtg, lateraldaughter, "NbLeaf", nb_leaf)
            mg3.set_gu_property(mtg, lateraldaughter, "UnitLeafArea", leaf_area)
    
    return newgus

def create_daughters_pruned(mtg, vid, apical, nblateral, burstdate, severity, intensity):
    import mangoG3 as mg3
    parentdirection = mg3.get_gu_normed_direction(mtg, vid)
    topposition = mg3.get_gu_top_position(mtg, vid)
    diam = mg3.get_gu_diameter(mtg, vid)
    newgus = []
    # newconnections = []
    lengths = []

    if apical:
        apicaldaughter = mtg.add_child(vid, edge_type = '<', UnitType = 'U', label = 'S'+str(int(mtg.label(vid)[1:])+1))        
        l = gu_length_pruned(severity)
        nb_leaf = gu_nb_leaf_pruned(l, severity)    
        leaf_area = gu_leaf_area_pruned(intensity, severity)
        mg3.set_gu_top_position(mtg, apicaldaughter, topposition + parentdirection * l )
        mg3.set_gu_diameter(mtg, apicaldaughter, diam)
        newgus.append(apicaldaughter)
        mg3.set_gu_property(mtg, apicaldaughter, "Regrowth", True)
        mg3.set_gu_property(mtg, apicaldaughter, "BurstDate", burstdate)
        mg3.set_gu_property(mtg, apicaldaughter, "NbLeaf", nb_leaf)
        mg3.set_gu_property(mtg, apicaldaughter, "UnitLeafArea", leaf_area)
        
    branching_angle = 60
    if nblateral > 0:
        for latiter, latdirection in enumerate(lateral_directions(parentdirection, branching_angle, nblateral)):
            lateralconnection = mtg.add_child(vid, edge_type = '+', label = 'S1', Position = topposition, UnitType = 'U')
            lateraldaughter   = mtg.add_child(lateralconnection, edge_type = '<', label = 'S2', UnitType = 'U', Diameter=diam/10.)
            l = gu_length_pruned(severity)
            nb_leaf = gu_nb_leaf_pruned(l, severity)    
            leaf_area = gu_leaf_area_pruned(intensity, severity)
            mg3.set_gu_top_position(mtg, lateraldaughter, topposition + latdirection * l )
            mg3.set_gu_diameter(mtg, lateraldaughter, diam)
            # newconnections.append(lateralconnection)
            newgus.append(lateraldaughter)
            mg3.set_gu_property(mtg, lateraldaughter, "Regrowth", True)
            mg3.set_gu_property(mtg, lateraldaughter, "BurstDate", burstdate)
            mg3.set_gu_property(mtg, lateraldaughter, "NbLeaf", nb_leaf)
            mg3.set_gu_property(mtg, lateraldaughter, "UnitLeafArea", leaf_area)
    
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

## Proposition modèle pour générer la longueur d'une UC fille ###########################################
def gu_length_unpruned(position):
    """ GU daugther length """
    probas = { eApical : (6.31, 2.95),
               eLateral : (4.61, 2.15)}

    probavalue, sd = probas[position]
    minval = 0.1 
    maxval = 18.50
    return normal_realization(probavalue, sd, maxval, minval)

def gu_length_pruned(severity):
    """ GU daugther length """
    probas = { n1 : (8.62, 4.15),
               n2 : (10.15, 4.55),
               n3 : (11.74, 4.89) }
    probavalue, sd = probas[severity]
    minval = 0.1 
    maxval = 28.50
    return normal_realization(probavalue, sd, maxval, minval)

## Proposition modèle pour générer le nombre de feuilles d'une UC fille ##############################

def gu_nb_leaf_unpruned(length, position):
    """ GU daugther number of leaf """
    intercept = 4.44
    sd = 2.5
    minval= 0 
    maxval= 18
    Coef_position = { eApical : 0., 
                      eLateral : -2.66 }
    probavalue = intercept + 0.61 * length + Coef_position[position] 
    return round(normal_realization(probavalue, sd, maxval, minval))

def gu_nb_leaf_pruned(length, severity):
    """ GU daugther number of leaf """
    intercept = 2.65
    sd = 3.29
    minval= 0 
    maxval= 17
    Coef_severity = { n1 : 0, 
                      n2 : -0.22,
                      n3 : -0.70 }
    probavalue = intercept + 0.46 * length + Coef_severity[severity] 
    return round(normal_realization(probavalue, sd, maxval, minval))

## Proposition modèle pour générer la surface individuelle de feuille pour une UC fille ##############################
def gu_leaf_area_unpruned(intensity, position):
    """ GU daugther number of leaf """
    intercept = 38.41
    sd = 20
    minval= 4.88
    maxval= 111.16
    
    Coef_position = { eApical : 0, 
                      eLateral : -16.58 }
        
    probavalue = intercept + 42.91 * intensity + Coef_position[position] 
    return normal_realization(probavalue, sd, maxval, minval)

def gu_leaf_area_pruned(intensity, severity):
    """ GU daugther number of leaf """
    intercept = 28.68
    sd = 15
    minval= 4.44
    maxval= 115.75

    Coef_severity = { n1 : 0, 
                      n2 : 1.34,
                      n3 : 3.45 }
        
    probavalue = intercept + 30.45 * intensity + Coef_severity[severity] 
    return normal_realization(probavalue, sd, maxval, minval)


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