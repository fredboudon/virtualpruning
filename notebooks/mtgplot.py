from openalea.mtg import MTG
import importlib
import openalea.mtg.plantframe as opf; importlib.reload(opf) 
from openalea.mtg.plantframe.plantframe import PlantFrame
from openalea.mtg.plantframe.dresser import DressingData
from openalea.plantgl.all import *
from randomgeneration import realisation
from mangoG3 import eApical, eLateral, gu_position, get_parent
from math import exp

def pos_prop(mtg):
    xx,yy,zz = mtg.property("XX"),mtg.property("YY"),mtg.property("ZZ")
    return dict ([(node,(x,-yy[node],-zz[node])) for node,x in list(xx.items())])

from math import radians, degrees

def orientation_prop(mtg):
    aa,bb,cc = mtg.property("AA"),mtg.property("BB"),mtg.property("CC")
    return dict ([(node,(radians(a),radians(bb[node]),radians(cc[node]))) for node,a in list(aa.items())])


matrixmethod = True


class HeigthColoring:
    def __init__(self):
        pass

    def set_mtg(self, mtg):
        self.heights = dict([(v,mtg.Height(v)) for v in mtg.vertices(scale=mtg.max_scale())])
        self.maxh = float(max(self.heights.values()))

    def __call__(self, turtle, vid):
        turtle.interpolateColors(1,2,self.heights[vid]/self.maxh)

class LeafColoring: 
    def __init__(self):
        pass

    def set_mtg(self, mtg):
        self.mtg = mtg
    def __call__(self, turtle, vid):
        turtle.setColor(2 if self.mtg.property('NbLeaf').get(vid,0) > 0 else 1)

class ClassColoring:
    def __init__(self):
        pass

    def set_mtg(self, mtg):
        self.mtg = mtg
        self.unittype = self.mtg.property('UnitType')
        self.deadunits = self.mtg.property('Dead')
        self.colors = { 'B' : 7, 'D' : 1, 'O' : 4, 'U' : 2}
        self.black = Material('BLACKMAT',(0,0,0))
    def __call__(self, turtle, vid):
        typechange =  (self.unittype[self.mtg.parent(vid)] != self.unittype[vid]) if self.mtg.parent(vid) else False
        gcon = turtle.getParameters().isGeneralizedCylinderOn()
        if typechange and gcon:
            turtle.stopGC()
        if vid in self.deadunits:
            turtle.setCustomAppearance(self.black)
        else:
            turtle.setColor(self.colors[self.unittype.get(vid,'B')])
        if typechange  and gcon:
            turtle.startGC()

class BlackColoring:
    def __init__(self):
        pass

    def set_mtg(self, mtg):
        self.mtg = mtg
        self.unittype = self.mtg.property('UnitType')
        self.color = Material('BLACKMAT',(0,0,0))
    def __call__(self, turtle, vid):
        turtle.setCustomAppearance(self.color)


def tofloat(v):
    import datetime
    if isinstance(v,datetime.timedelta):
        return float(v.days)
    else:
        return float(v)

def return_default_if_none(value, default):
    if value is None : return default
    else : return value

class PropertyGradientColoring: 
    def __init__(self, prop, listids = None, tofloat = tofloat, minvalue = None, maxvalue = None, cmap='jet', defaultvalue = None):
        import mangoG3 as mg
        self.prop = prop
        self.tofloat = tofloat
        self.listids = set(listids) if not listids is None else listids
        self.cmap = PglMaterialMap(min(self.prop) if minvalue is None else minvalue, max(self.prop) if maxvalue is None else maxvalue,cmap)
        self.defaultvalue = defaultvalue
        self.made = set()
        self.minvalue = minvalue
        self.maxvalue = maxvalue

    def get_prop(self):
        if type(self.prop) == str:
            return lambda v : mg.get_gu_property(self.mtg, v, self.prop, self.defaultvalue)
        elif hasattr(self.prop,'__getitem__'):
            return lambda v : self.prop.get(v, self.defaultvalue)
        elif callable(self.prop):
            if self.defaultvalue is None:
                return self.prop
            else:
                return lambda v : return_default_if_none(self.prop(v), self.defaultvalue)
        else :
            return lambda v : self.prop

    def prepare_turtle(self, turtle):
        pass
        #from openalea.plantgl.all import Material
        #turtle.setMaterial(1, self.defaultcolor) # ,transparency=0.8))
        #turtle.setMaterial(10,self.mincolor)
        #turtle.setMaterial(11,self.maxcolor)

    def done(self, turtle):
        pass #print self.listids.symmetric_difference(self.made)

    def set_mtg(self, mtg):
        import mangoG3 as mg
        self.mtg = mtg
        f = self.get_prop()
        if not self.listids is None:
            toconsider = [f(vid) for vid in self.listids]  
        else:
            toconsider = self.prop.values
        #if self.minvalue is None:
        #    self.minvalue = min(toconsider)
        #if self.maxvalue is None:
        #    self.maxvalue = max(toconsider)
        #print(self.minvalue, self.maxvalue)
        #self.deltavalue = self.tofloat(self.maxvalue - self.minvalue)


    def __call__(self, turtle, vid):
        import mangoG3 as mg
        f = self.get_prop()
        if self.listids is None or vid in self.listids :
            v = f(vid)
            if not v is None:
                turtle.setCustomAppearance(self.cmap(v))
            else:
                return False
        else:
            return False

def leafsmb():

    leafdiam = QuantisedFunction(NurbsCurve2D(Point3Array([Vector3(0,0.0,1),Vector3(0.239002,1.00091,1),Vector3(0.485529,0.991241,1),Vector3(0.718616,1.00718,1),Vector3(0.877539,0.231273,1),Vector3(1,0.0,1)])))
    leafpath = NurbsCurve2D(Point3Array([(-0.5, 0, 1),(-0.145022, -0.0735931, 1),(0.0844156, -0.212121, 1),(0.123377, -0.497835, 1)]))    
    leafsection = NurbsCurve2D(Point3Array([Vector3(-0.5,0.15,1),Vector3(-0.5,0.1,1),Vector3(-0.1,-0.1,1),Vector3(0,0.2,1),Vector3(0.1,-0.1,1),Vector3(0.5,0.1,1),Vector3(0.5,0.15,1)]))
    leafsection.stride = 1
    length = 1

    turtle = PglTurtle()
    turtle.push()
    turtle.down(60)
    turtle.startGC().sweep(leafpath,leafsection,length,length/3.,length*0.24,leafdiam).stopGC()
    turtle.pop()

    leafsmb = turtle.getScene()[0].geometry

    t = Tesselator()
    leafsmb.apply(t)
    leafsmb = t.result

    surfs = surfaces(leafsmb.indexList, leafsmb.pointList)
    toremove = [i for i,s in enumerate(surfs) if s < 1e-5]
    leafsmb.indexList = leafsmb.indexList.opposite_subset(toremove)
    assert leafsmb.isValid()

    #print(abs(surface(leafsmb) - 0.18))
    assert abs(surface(leafsmb) - 0.18) < 0.1

    leafsmb.name = 'leaf'
    return leafsmb

axis, north = (0,0,1), -(90-53)

def rotate_scene(sc, axis = axis, angle = north):
    return Scene([Shape(AxisRotated(axis,angle,sh.geometry),sh.appearance,sh.id,sh.parentId) for sh in sc])

def internode_length_distribution(nb_internodes, gu_length):
  """ Internode length distribution """
  if nb_internodes <= 1 : return [ gu_length ]
  lengths = [exp(-2.64 * i / float(nb_internodes-1)) for i in range(nb_internodes)]
  scaling = gu_length/ sum(lengths)
  return [l*scaling for l in lengths]

def length_before_first_leaf(position, final_length_gu):
  #length of space before the first leaf
  if position == eApical:
    from numpy.random import gamma
    # LEPF = gauss(2.63,1.72)
    LEPF = min(realisation(2.007, 0.763, 0, 8, gamma)[0],final_length_gu)

  else: # Lateral case
    #length of space before the first leaf depend of GU's length
    LEPF = final_length_gu * 0.38 + 0.88
  return LEPF

def representation(mtg, focus = None, 
                        colorizer = ClassColoring(), 
                        leaves = False, 
                        wood = True, 
                        gc = True,
                        sensors = False,
                        fieldoriented = True, 
                        leafreorientation = 30,
                        sensorsize = 5,
                        todate = None):
    import inspect
    if inspect.isclass(colorizer):
        colorizer = colorizer()

    posproperty = pos_prop(mtg)
    orientations = orientation_prop(mtg)

    div10 = lambda x : abs(x/10.) if x else x
    minus = lambda x : -x if x else x

    topdia = lambda x: mtg.property('Diameter').get(x)
    diameters = {}
    for lvid in mtg.vertices(scale=3):
        td = topdia(lvid)
        if td is None :
            if not mtg.parent(lvid) is None:
                diameters[lvid] = diameters[mtg.parent(lvid)]
        else:
            diameters[lvid] = td

    zz = lambda x: minus(mtg.property('ZZ').get(x))
    yy = lambda x: minus(mtg.property('YY').get(x))
    TopDiameter = lambda v: diameters.get(v)

    pf = PlantFrame(mtg, 
                    # BottomDiameter=botdia, 
                    TopDiameter=TopDiameter, 
                    YY = yy, 
                    ZZ = zz, 
                    origin=posproperty[mtg.roots(3)[0]],
                    scale=3
                    )

    #diameters = pf.compute_diameters()
    #pf.points = dict([(node,Vector3(pos)-pf.origin) for node,pos in pf.points.items()])
    pf.points = mtg.property('Position') # dict([(node,Vector3(pos)) for node,pos in pf.points.items()])

    colorizer.set_mtg(mtg)

    leaf_length_distrib = { eApical  : ( 17.06 , 2.7) ,
                            eLateral : ( 14.87 , 2.7) }
    leaf_area_length_ratio   = 2.3594
    ePruned, eUnPruned = 0,1
    leaf_area_length_ratio_regrowth   = { ePruned : 2.61, eUnPruned : 2.74}

    InternodesLength = mtg.property('InternodesLength')
    UnitLeafArea = mtg.property('UnitLeafArea')
    Dead = mtg.property('Dead')
    
    todraw = None
    if not focus is None:
        if type(focus) == str:
            vids = [lvid for lvid,rem in list(mtg.property('Id').items()) if rem == focus]
            if len(vids) >= 1:
                todraw = set()
                for v in vids:
                    print('Display', v)
                    todraw |= set(mtg.Descendants(v))
                    todraw |= set(mtg.Ancestors(v))
                focus = set(vids)
            else:
                print('Cannot find', focus)
                focus = None
        else:
            todraw = set(mtg.Descendants(focus))
            todraw |= set(mtg.Ancestors(focus))
            focus = set([focus])

    if todate:
        import pandas as pd
        todate = pd.Timestamp(todate)

    def leaf(turtle, length):
        turtle.push()
        turtle.down(60)
        # turtle.setColor(2)
        turtle.startGC().sweep(leafpath,leafsection,length,length/10.,length*0.24,leafdiam).stopGC()
        turtle.pop()

    def plantframe_visitor(g, v, turtle):
        if todraw and not v in todraw: return
        if todate :
            digitdate = g.property('DigitDate').get(v)
            if digitdate and digitdate > todate: return
            burstdate = g.property('BurstDate').get(v)
            if burstdate and burstdate > todate: return


        from random import gauss
        radius = diameters.get(v)
        pt = pf.points.get(v)

        unittype = g.property('UnitType').get(v,'B')
        dead = Dead.get(v,False)

        if pt:
            if sensors:
                turtle.setId(Shape.NOID-1)
            else:
                turtle.setId(v)
            if 'M' in unittype:
                pass
            elif 'C' in unittype:
                turtle.lineTo(pt)
            else:
                if focus and v in focus:
                    turtle.setColor(3)
                    colortodraw = True
                else:
                    colortodraw = colorizer(turtle, v)
                    if colortodraw is None: 
                        colortodraw = True
                if g.edge_type(v) == '<':
                    nbleaf = g.property('NbLeaf').get(v,0)
                    if not colortodraw :
                        turtle.move(pt)
                        turtle.setWidth(radius)                        
                    elif not leaves or nbleaf == 0 or dead:
                        if wood:
                            turtle.lineTo(pt, radius)
                        else:
                            turtle.pinpoint(pt)
                            turtle.move(pt)
                            turtle.setWidth(radius)
                    else:
                        turtle.pinpoint(pt)
                        parent = g.parent(v)
                        parentpos = pf.points[parent]
                        length = norm(pt-parentpos)

                        position = eApical if gu_position(g, v) else eLateral

                        internodes_length = InternodesLength.get(v)
                        if internodes_length is None:
                            LEPF = length_before_first_leaf(position, length)
                            internodes_length = [LEPF] + internode_length_distribution(nbleaf-1, length-LEPF)
                            InternodesLength[v] = internodes_length

                        unitleafarea = UnitLeafArea.get(v)
                        if unitleafarea is None:
                            unitleafarea = gauss(*leaf_length_distrib[position])*leaf_area_length_ratio
                            UnitLeafArea[v] = unitleafarea
                        lalratio = leaf_area_length_ratio
                        if v in mtg.property('Regrowth'):
                            lalratio = leaf_area_length_ratio_regrowth[mtg.property('pruned').get(get_parent(g,v)) is None]

                        parentradius = diameters.get(parent)
                        segdiaminc = (radius-parentradius)/nbleaf
                        for i, ilength in enumerate(internodes_length):
                            if wood:
                                turtle.F(ilength,parentradius+segdiaminc*(i+1))
                            else:
                                turtle.f(ilength)
                            turtle.rollR(144)
                            turtle.push()
                            if not leafreorientation is None:
                                turtle.down(90)
                                turtle.rollToVert()
                                turtle.rollToHorizontal()
                                turtle.up(90-leafreorientation)
                            positionratio = 1
                            if nbleaf > 2 and i >= nbleaf-2:
                                if i == nbleaf-2 : positionratio = 0.52
                                else: positionratio = 0.38
                            turtle.surface('leaf', unitleafarea*positionratio/lalratio)
                            turtle.pop()
                        turtle.setWidth(radius)
                    if sensors and nbleaf > 0:
                        turtle.push()
                        turtle.setId(v)
                        turtle.setHead((0,0,1),((1,0,0)))
                        turtle.f(0.5)
                        turtle.down(90)
                        turtle.f(-sensorsize/2)
                        turtle.setColor(20)
                        turtle.setWidth(sensorsize/2)
                        turtle.quad(sensorsize)
                        turtle.pop()
                else:
                    if gc : turtle.stopGC()
                    if g.edge_type(v) == '+':
                        parent = mtg.parent(v)
                        grandparent = mtg.parent(parent)
                        parentpos = pf.points[parent]
                        grandparentpos = pf.points[grandparent]
                        bpt, l, u = closestPointToSegment(pt, grandparentpos, parentpos)
                        turtle.move(bpt)
                        turtle.setWidth(radius)
                        if gc : turtle.startGC()
                        #if norm(bpt-pt) > 1e-3:
                        #    turtle.lineTo(pt)
                    else:
                        turtle.move(pt)
                        turtle.setWidth(radius)
                        if gc : turtle.startGC()

    turtle = PglTurtle()
    turtle.setSurface('leaf', leafsmb())
    turtle.setColorAt(7,(30,22,7))
    turtle.setColorAt(1,(65,45,15))
    turtle.setColorAt(20,(0,0,0))
    if hasattr(colorizer,'prepare_turtle'):
        colorizer.prepare_turtle(turtle)
    sc = pf.plot(origins=[(0,0,0)],visitor=plantframe_visitor,gc=gc, turtle = turtle, display = False)
    if hasattr(colorizer,'done'):
        colorizer.done(turtle)
    if fieldoriented:
        sc = rotate_scene(sc, axis, north)
    return sc


def retrievedates(mtg):
    import numpy as np
    import pandas as pd
    dates = np.unique(list(mtg.property('DigitDate').values()))+np.unique(list(mtg.property('BurstDate').values()))
    dates.sort()
    dates = np.unique([pd.Timestamp(d.year,d.month,d.day,23,59) for d in dates])
    return dates

def treecentroid(mtg, date = None):
    from openalea.plantgl.all import Vector3, Point3Array
    xx = mtg.property('XX')
    yy = mtg.property('YY')
    zz = mtg.property('ZZ')
    if date:
        date = date.to_pydatetime().date()
        centroid = Point3Array([(xx[vid],-yy[vid],-zz[vid]) for vid in mtg.vertices(scale=3) if mtg.property('DigitDate')[vid].to_pydatetime().date() == date]).getCenter()
    else:
        centroid = Point3Array([(xx[vid],-yy[vid],-zz[vid]) for vid in mtg.vertices(scale=3)]).getCenter()

    return centroid


plot3D = True

def display(scene):
    from pgljupyter import SceneWidget
    from openalea.plantgl.all import BoundingBox
    bbx = BoundingBox(scene)
    sw = norm(bbx.getSize())
    return SceneWidget(scene, size_world=sw, size_display = (800,600))

def plot_tree(mtg, **kwd):
    if plot3D :
        sc = representation(mtg, **kwd)
        return display(sc)
 

def plot_terminal_diameter(mtg, leaves = True, minvalue = None, maxvalue = None):
    if plot3D :
        import mangoG3 as mg
        vids = [vid for vid in mg.get_all_terminal_gus(mtg) if not mg.was_previously_pruned(mtg,vid)]
        colorizer= PropertyGradientColoring(lambda v : mg.get_gu_diameter(mtg,v), vids, minvalue = minvalue, maxvalue = maxvalue)
        sc = representation(mtg, colorizer=colorizer , gc = False, leaves=leaves )
        # Viewer.display(sc)
        return display(sc)


def projectview(scene, shvalue):
    from openalea.plantgl.scenegraph.colormap import PglMaterialMap
    if isinstance(shvalue, dict):
        cmap = PglMaterialMap(min(shvalue.values()),max(shvalue.values()))
    else:
        cmap = PglMaterialMap(min(shvalue),max(shvalue))
    return Scene([Shape(sh.geometry, cmap(shvalue[sh.id]),sh.id,sh.parentId) for sh in scene])

def plot_projection(scene, shvalue):
    return display(projectview(scene, shvalue))