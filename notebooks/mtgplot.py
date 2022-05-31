from openalea.mtg import MTG
import importlib
import openalea.mtg.plantframe as opf; importlib.reload(opf) 
from openalea.mtg.plantframe.plantframe import PlantFrame
from openalea.mtg.plantframe.dresser import DressingData
from openalea.plantgl.all import *

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
        self.colors = { 'B' : 7, 'D' : 1, 'O' : 4, 'U' : 2}
    def __call__(self, turtle, vid):
        typechange =  (self.unittype[self.mtg.parent(vid)] != self.unittype[vid]) if self.mtg.parent(vid) else False
        gcon = turtle.getParameters().isGeneralizedCylinderOn()
        if typechange and gcon:
            turtle.stopGC()
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

class PropertyGradientColoring: 
    def __init__(self, prop, listids = None, tofloat = tofloat, minvalue = None, maxvalue = None, mincolor = Material((200,200,0)), maxcolor = Material((200,0,0)), defaultcolor = Material((45,65,15))):
        import mangoG3 as mg
        self.prop = prop
        self.tofloat = tofloat
        self.listids = set(listids) if listids else listids
        self.mincolor = mincolor
        self.maxcolor = maxcolor
        self.defaultcolor = defaultcolor
        self.made = set()
        self.minvalue = minvalue
        self.maxvalue = maxvalue

    def get_prop(self):
        if type(self.prop) == str:
            return lambda v : mg.get_gu_property(self.mtg, v, self.prop)
        elif hasattr(self.prop,'__getitem__'):
            return lambda v : self.prop[v]
        elif callable(self.prop):
            return self.prop
        else :
            return lambda v : self.prop

    def prepare_turtle(self, turtle):
        from openalea.plantgl.all import Material
        turtle.setMaterial(1, self.defaultcolor) # ,transparency=0.8))
        turtle.setMaterial(10,self.mincolor)
        turtle.setMaterial(11,self.maxcolor)

    def done(self, turtle):
        pass #print self.listids.symmetric_difference(self.made)

    def set_mtg(self, mtg):
        import mangoG3 as mg
        self.mtg = mtg
        f = self.get_prop()
        if not self.listids is None:
            toconsider = [f(vid) for vid in self.listids]  
        else:
            toconsider = list(prop.values())
        if self.minvalue is None:
            self.minvalue = min(toconsider)
        if self.maxvalue is None:
            self.maxvalue = max(toconsider)
        #print(self.minvalue, self.maxvalue)
        self.deltavalue = self.tofloat(self.maxvalue - self.minvalue)


    def __call__(self, turtle, vid):
        import mangoG3 as mg
        f = self.get_prop()
        if self.listids is None or vid in self.listids :
            v = f(vid)
            if self.minvalue <= v <= self.maxvalue:
                self.made.add(vid)
                value = self.tofloat(v-self.minvalue)/self.deltavalue
                turtle.interpolateColors(10,11,value)
            else:
                turtle.setColor(1)
        else:
            turtle.setColor(1)

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

def representation(mtg, focus = None, colorizer = ClassColoring(), leaves = False, wood = True, gc = True, todate = None):
    import inspect
    if inspect.isclass(colorizer):
        colorizer = colorizer()

    posproperty = pos_prop(mtg)
    orientations = orientation_prop(mtg)

    div10 = lambda x : abs(x/10.) if x else x
    minus = lambda x : -x if x else x

    topdia = lambda x: mtg.property('Diameter').get(x)
    diameters = {}
    for vid in mtg.vertices(scale=3):
        td = topdia(vid)
        if td is None :
            if not mtg.parent(vid) is None:
                diameters[vid] = diameters[mtg.parent(vid)]
        else:
            diameters[vid] = td

    zz = lambda x: minus(mtg.property('ZZ').get(x))
    yy = lambda x: minus(mtg.property('YY').get(x))
    TopDiameter = lambda v: diameters.get(v)

    pf = PlantFrame(mtg, 
                    # BottomDiameter=botdia, 
                    TopDiameter=TopDiameter, 
                    YY = yy, 
                    ZZ = zz, 
                    origin=posproperty[mtg.roots(3)[0]]
                    )

    #diameters = pf.compute_diameters()
    #pf.points = dict([(node,Vector3(pos)-pf.origin) for node,pos in pf.points.items()])
    pf.points = mtg.property('Position') # dict([(node,Vector3(pos)) for node,pos in pf.points.items()])

    colorizer.set_mtg(mtg)

    meanleaflength = 20
    leaf_length_distrib = { '<'  : ( 17.06 , 2.7) ,
                            '+' : ( 14.87 , 2.7) }
    
    todraw = None
    if not focus is None:
        if type(focus) == str:
            vids = [vid for vid,rem in list(mtg.property('Id').items()) if rem == focus]
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
        if todate and g.property('DigitDate')[v] > todate: return

        from random import gauss
        radius = diameters.get(v)
        pt = pf.points.get(v)

        unittype = g.property('UnitType').get(vid,'B')
        if pt:
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
                    if colortodraw is None: colortodraw = True
                if g.edge_type(v) == '<':
                    nbleaf = mtg.property('NbLeaf').get(v,0)
                    if not colortodraw :
                        turtle.moveTo(pt, radius)                        
                    elif not leaves or nbleaf == 0:
                        if wood:
                            turtle.lineTo(pt, radius)
                        else:
                            turtle.pinpoint(pt)
                            turtle.move(pt)
                            turtle.setWidth(radius)
                    else:
                        turtle.pinpoint(pt)
                        parent = mtg.parent(v)
                        parentpos = pf.points[parent]
                        length = norm(pt-parentpos)
                        seglength = length/nbleaf
                        parentradius = diameters.get(parent)
                        segdiaminc = (radius-parentradius)/nbleaf
                        for i in range(nbleaf):
                            if wood:
                                turtle.F(seglength,parentradius+segdiaminc*(i+1))
                            else:
                                turtle.f(seglength)
                            turtle.rollR(144)
                            turtle.surface('leaf', gauss(*leaf_length_distrib[mtg.edge_type(mtg.parent(v))]))
                            #leaf(turtle, gauss(*leaf_length_distrib[mtg.edge_type(mtg.parent(v))]) )
                        turtle.setWidth(radius)

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
    if hasattr(colorizer,'prepare_turtle'):
        colorizer.prepare_turtle(turtle)
    sc = pf.plot(origins=[(0,0,0)],visitor=plantframe_visitor,gc=gc, turtle = turtle, display = False)
    if hasattr(colorizer,'done'):
        colorizer.done(turtle)
    return sc


def retrievedates(mtg):
    import numpy as np
    import pandas as pd
    dates = np.unique(list(mtg.property('DigitDate').values()))
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
    cmap = PglMaterialMap(min(shvalue),max(shvalue))
    return Scene([Shape(sh.geometry, cmap(shvalue[sh.id]),sh.id,sh.parentId) for sh in scene])

def plot_projection(scene, shvalue):
    return display(projectview(scene, shvalue))