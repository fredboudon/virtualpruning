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
        print(self.minvalue, self.maxvalue)
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


def representation(mtg, focus = None, colorizer = ClassColoring(), leaves = False, gc = True, todate = None):
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

    leafdiam = QuantisedFunction(NurbsCurve2D(Point3Array([Vector3(0,0.0846264,1),Vector3(0.239002,1.00091,1),Vector3(0.485529,0.991241,1),Vector3(0.718616,1.00718,1),Vector3(0.877539,0.231273,1),Vector3(1,0.00332359,1)])))
    leafpath = NurbsCurve2D(Point3Array([(-0.5, 0, 1),(-0.145022, -0.0735931, 1),(0.0844156, -0.212121, 1),(0.123377, -0.497835, 1)]))    
    leafsection = NurbsCurve2D(Point3Array([Vector3(-0.508209,0.16873,1),Vector3(-0.515031,0.138195,1),Vector3(-0.198373,-0.0924227,1),Vector3(-0.00298323,0.188761,1),Vector3(0.0897461,-0.106293,1),Vector3(0.555704,0.0979703,1),Vector3(0.545047,0.12817,1)]))

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
                        turtle.lineTo(pt, radius)
                    else:
                        turtle.pinpoint(pt)
                        parent = mtg.parent(v)
                        parentpos = pf.points[parent]
                        length = norm(pt-parentpos)
                        seglength = length/nbleaf
                        parentradius = diameters.get(parent)
                        segdiaminc = (radius-parentradius)/nbleaf
                        for i in range(nbleaf):
                            turtle.F(seglength,parentradius+segdiaminc*(i+1))
                            turtle.rollR(144)
                            turtle.surface('leaf', gauss(*leaf_length_distrib[mtg.edge_type(mtg.parent(v))]))
                            #leaf(turtle, gauss(*leaf_length_distrib[mtg.edge_type(mtg.parent(v))]) )

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
                        if norm(bpt-pt) > 1e-3:
                            turtle.lineTo(pt)
                    else:
                        turtle.move(pt)
                        turtle.setWidth(radius)
                        if gc : turtle.startGC()


    turtle = PglTurtle()
    leaf(turtle, 1)
    leafsmb = turtle.getScene()[0].geometry
    leafsmb.name = 'leaf'

    turtle = PglTurtle()
    turtle.setSurface('leaf', leafsmb)
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

def camerapositions(mtg, dist = 1200):
    centroid = treecentroid(mtg)
    result = []
    for d in retrievedates(mtg):
        dcentroid = treecentroid(mtg,d)
        dcentroid.z = centroid.z
        dir = direction(dcentroid-centroid)
        camposition = centroid + dir*dist
        result.append((d,camposition))
    tomerge = []
    for i in range(len(result)-1):
        j = i+1
        ok = True
        while angle(result[i][1], result[j][1],(0,0,1)) < 0:
            ok = False
            j = j+1
        if ok is False:
            tomerge.append(list(range(i,j)))
            i = j
    for atm in tomerge:
        npos = Point3Array([result[did][1] for did in atm]).getCenter()
        for did in atm: 
            result[did] = (result[did][0],npos)
    return dict(result), centroid


def cameraangles(mtg, dist = 1200):
    cc, pos = camerapositions(mtg, dist)
    cc = list(cc.items())
    cc.sort()
    ref = Vector3(cc[0][1])-pos
    cc2 = [(d,angle(ref,Vector3(cci)-pos,(0,0,1))) for d,cci in cc]
    return dict(cc2),pos,Vector3(cc[0][1])

def interpolatecamerapositions(mtg, dist = 1200, nbposition = 25):
    cc, pos = camerapositions(mtg, dist)
    cc = list(cc.items())
    cc.sort()
    ref = Vector3(cc[0][1])-pos
    res = []
    for i in range(len(cc)):
        d,cci = cc[i]
        d,ccj = cc[(i+1)%len(cc)]
        for x in range(nbposition):
            ccd = (cci*(nbposition-x) + ccj*x)/nbposition
            ccd = pos + direction(ccd-pos)*dist
            res.append(ccd)
    res.append(cc[0][1])
    return res, pos

def bezinterpolatecamerapositions(mtg, dist = 1200, nbposition = 25):
    res, pos = interpolatecamerapositions(mtg, dist = dist, nbposition = 4)
    nbc = NurbsCurve(Point4Array(res,1))
    nbc.stride = ((len(res)-1)*nbposition/4)-1
    d = Discretizer()
    nbc.apply(d)
    return d.result.pointList, pos

def plot_camerapositions(mtg, dist = 1200, interpolation = True):
    if interpolation:
        cc1, pos = bezinterpolatecamerapositions(mtg, dist)
    else:
        cc, pos = camerapositions(mtg, dist)
        cc = list(cc.items())
        cc.sort()
        cc1 = [b for a,b in cc]
    Viewer.add(Polyline(cc1))

def timeplot(mtg, colorizer = ClassColoring, leaves = False, gc = True):
    from random import seed
    from openalea.plantgl.all import Viewer
    campos, centroid = camerapositions(mtg)
    for d in retrievedates(mtg):
        seed(0)
        Viewer.redrawPolicy = False
        sc = plot(mtg, colorizer=colorizer, leaves=leaves, gc=gc, todate = d)
        Viewer.redrawPolicy = True
        Viewer.camera.lookAt(campos[d], centroid)


header = '''
/*
 * A povray file generated with GEOM.
 * Example of use : povray -Ifile.pov -Ofile.png +FN +H600 +W800 +A.
 * File Generated with PlantGL 3D Viewer.
 */


#ifndef (__camera_definition__)
#declare __camera_definition__ = true;

#declare Location = %s;
#declare Right = %s;
#declare NRight = %s;
#declare LookAt = %s;
#declare RotAngle = %f;
#declare Direction = %s;
#declare Date = "%s";

camera {
   perspective
    location Location
    sky <0,0,1>
    right Right
    look_at LookAt
}

light_source {
     <27,16,800>
    color rgb 0.9
}



background { color rgb <1,1,1> }

plane { <0,0,1> 0 
        clipped_by { plane {<-1,0,0> 1500     rotate <0,0,RotAngle> } }
        texture {
            pigment {
               color rgb 1
            }
            finish {
              ambient 0.4
              diffuse 1
              specular 0
            }

       }
}


#end // __camera_definition__

'''
datedef = '''

object {
    text {
    ttf "crystal.ttf" Date 0.01, 0
    pigment { color rgb <0,0,0> }

    }
    rotate <90,0,89+RotAngle>
    translate Location+(Direction*10)+<0,0,-4>+(NRight*-3)
}

'''

def generate(mtg, colorizer = ClassColoring, leaves = False, gc = True, nbsteps = 1, withdate = False):
    import numpy as np
    import pandas as pd
    from math import degrees
    import os
    if not os.path.exists('reconstructionmovie'):
        os.makedirs('reconstructionmovie')
    os.chdir('reconstructionmovie')
    from random import seed
    dates = np.unique(list(mtg.property('DigitDate').values()))
    dates.sort()
    dates = np.unique([pd.Timestamp(d.year,d.month,d.day,23,59) for d in dates])
    for d in dates:
        fname = 'mangodigit'+str(d.year)+'-'+str(d.month)+'-'+str(d.day).zfill(2)
        if leaves: fname += '_leafy'
        fname +='.pov'
        if not os.path.exists(fname):
            seed(0)
            sc = plot(mtg, colorizer=colorizer, leaves=leaves, gc=gc, todate = d, display = False)
            sc.save(fname)
    campos, cpos = bezinterpolatecamerapositions(mtg, dist = 800, nbposition=nbsteps)
    #print len(campos), len(dates)*nbsteps
    assert len(campos) == len(dates)*nbsteps
    v3tostr = lambda v : '<'+str(v.x)+','+str(v.y)+','+str(v.z)+'>'
    for i,cam in enumerate(campos):
        povstream = file('view_'+str(i).zfill(4)+'.pov','w')
        right = direction(cross(direction(cam-cpos),(0,0,1)))*16/9.
        angleSoil = degrees(angle((1,0,0),direction(cam-cpos),(0,0,1)))
        angleRight = degrees(angle((1,0,0),direction(right),(0,0,1)))
        nright = right
        if np.sign(angleSoil) == np.sign(angleRight) == -1 or angleSoil-angleRight < 0:
            print('not change',i)
        else:
            right *= -1
        #print i, cam, right, angleSoil, angleRight, angleSoil-angleRight
        d = dates[i // nbsteps]
        povstream.write(header % (v3tostr(cam),v3tostr(right),v3tostr(nright),v3tostr(cpos), angleSoil, v3tostr(direction(cpos-cam)), str(d.day).zfill(2)+'/'+str(d.month).zfill(2)))
        if withdate : povstream.write(datedef)
        povstream.write('#include "mangodigit'+str(d.year)+'-'+str(d.month)+'-'+str(d.day).zfill(2)+('_leafy' if leaves else '')+'.pov"\n\n\n')
        povstream.close()
    povstream = file('mangoreconstruction.pov','w')
    for i in range(len(campos)):
        povstream.write('#if (clock = %i)\n #include "view_%s.pov"\n#end\n' % (i,str(i).zfill(4)))
    povstream.close()
    cmd = "povray -Imangoreconstruction.pov +H720 +W1280 +FN +A +KFI0 +KI0 +KFF"+str(len(campos)-1)+" +KF"+str(len(campos)-1)+" +L/opt/local/share/povray-3.7/include/ -GA"
    from time import time
    t = time()
    os.system(cmd)
    print('Done in',time()-t, 'sec.')

    #for i in xrange(len(campos)):
    #    cmd = "povray -Iview_"+str(i).zfill(4)+".pov +H720 +W1280 +FN +A -GA"
    #    print cmd
    #    os.system(cmd)
    os.chdir(os.pardir)
