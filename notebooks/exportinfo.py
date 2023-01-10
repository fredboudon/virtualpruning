import openalea.plantgl.all as pgl
from openalea.plantgl.all import *
from openalea.mtg import MTG
import mangoG3 as mg3
import pandas as pd
from math import *

def extend_mtg(mtg, inplace=True, *tableproperties, **properties):
    if not inplace:
        from copy import deepcopy
        mtg = deepcopy(mtg)
    for name, values in properties.items():
        mtg.property(name).update(values)
    for tableproperty in tableproperties:
        for name, values in tableproperty.items():
            mtg.property(name).update(values)
    return mtg

def compute_topooder(mtg):
    root = mtg.roots(scale=mtg.max_scale())[0]
    gurootdepth = mg3.gus_depth_from_root(mtg)
    leafdepth = {}
    rootdepth = {}
    for vid in mg3.get_all_gus(mtg):
        leafdepth[vid] = mg3.get_gu_terminal_min_depth(mtg,vid)
        rootdepth[vid] = gurootdepth[vid]
    return leafdepth, rootdepth

def extend_mtg_with_topooder(mtg):
    leafdepth, rootdepth = compute_topooder(mtg)
    extend_mtg(mtg, TopoLeafDepth=leafdepth, TopoRootDepth=rootdepth)

def compute_leafinfo(mtg, scene):
    pos = {}
    XX = {}
    YY = {}
    ZZ = {}
    El = {}
    Az = {}

    d = pgl.Discretizer()
    shapedict = scene.todict()
    for id in mg3.get_all_leaves(mtg):
        if id in shapedict:
            sh = shapedict[id][0]
            sh.apply(d)
            tr = d.result
            p = tr.pointList.getCenter()
            pos[id] = p
            XX[id] = p.x
            YY[id] = p.y
            ZZ[id] = p.z

            tr.normalPerVertex
            tr.computeNormalList()
            nml = pgl.Vector3.Spherical(sum(tr.normalList,pgl.Vector3(0,0,0)) / len(tr.normalList))
            El[id] = degrees(nml.phi)
            Az[id] = degrees(nml.theta)
    return pos, XX, YY, ZZ, El, Az

def extend_mtg_with_leafinfo(mtg, scene):
    pos, XX, YY, ZZ, El, Az = compute_leafinfo(mtg, scene)
    assert len(pos) > 0
    extend_mtg(mtg,LeafPosition=pos, XX=XX, YY=YY, ZZ=ZZ, El=El, Az=Az)
    return mtg

def spherical_depth(mtg, scale=3):
    if scale == 4:
        pos = mtg.property('LeafPosition')
    else:
        pos = mg3.get_gu_mean_positions(mtg)

    assert len(pos) > 0
    points = pgl.Point3Array(list(pos.values()))
    geomorder = {}
    center = points.getCenter()
    extent = points.getExtent()
    center.z -= extent.z /2
    extent = pgl.norm(extent)
    ids = list(pos.keys())

    def toThetaPhi(v):
        vs = Vector3.Spherical(v) 
        return Vector2(vs.theta, vs.phi)

    anglegrid = pgl.Point2Grid([toThetaPhi(v-center) for v in pos.values()],20)

    for id, opt in pos.items():
        dir = opt-center
        dist = dir.normalize()
        n = dist/extent
        b = center + dir*extent

        pids = anglegrid.query_ball_point(toThetaPhi(dir),pi/8)

        geomorder[id] = 0
        for lid in pids:
            cp, d, u = pgl.closestPointToSegment(points[lid], center, b)
            if d < 20 :
                if u > n:
                    geomorder[id] += 1
    return geomorder

def extend_mtg_with_sphericaldepth(mtg, scale=3):
    geomorder = spherical_depth(mtg, scale=scale)
    extend_mtg(mtg, SphericalDepth=geomorder)

def horizontal_depth(mtg, scale=3):
    if scale == 4:
        pos = mtg.property('LeafPosition')
    else:
        pos = mg3.get_gu_mean_positions(mtg)
    assert len(pos) > 0
    points = pgl.Point3Array(list(pos.values()))
    point2grid = pgl.Point2Grid([Vector2(p.x,p.y) for p in pos.values()],20)
    indices = list(pos.keys())
    geomorder = {}
    zmax = points[points.getZMaxIndex()].z
    for id, opt in pos.items():
        dir = (0,0,1)
        n = opt.z/zmax
        opti = pgl.Vector3(opt.x,opt.y,0)
        optj = pgl.Vector3(opt.x,opt.y,zmax)
        geomorder[id] = 0
        pids = point2grid.query_ball_point((opt.x,opt.y), 20)
        for lid in pids:
            cp, d, u = pgl.closestPointToSegment(points[lid], opti, optj)
            if d < 20 and u > n:
                geomorder[id] += 1
    return geomorder

def extend_mtg_with_horizontalorder(mtg, scale=3):
    geomorder = horizontal_depth(mtg, scale=scale)
    extend_mtg(mtg, HorizontalDepth=geomorder)


def extend_mtg_with_characteristic(mtg, inplace = True):
    if not inplace:
        from copy import deepcopy
        mtg = deepcopy(mtg)
    #import mtgplot as mp
    #scene = mp.representation(mtg,leaves=True, wood=False)
    #extend_mtg_with_leafinfo(mtg, scene)
    extend_mtg_with_topooder(mtg)
    extend_mtg_with_sphericaldepth(mtg)
    extend_mtg_with_horizontalorder(mtg)
    return mtg

def export_as_dataframe(mtg):
    from pandas import DataFrame
    vertices = set(mg3.get_all_gus(mtg))
    proptoexport = {}
    for name,prop in mtg.properties().items():
        if name not in ['edge_type','_line','Position', 'XX', 'YY', 'ZZ', 'AA', 'BB', 'CC']:
            proptoexport[name] = dict([(vid,v) for vid,v in prop.items() if vid in vertices])
    proptoexport['parent'] = dict([(vid, mg3.get_parent(mtg,vid)) for vid in vertices])
    def add_position(proptoexport, name, pos):
        for attname in ['x','y','z']:
            proptoexport[name+attname.upper()] = dict([(vid, getattr(pos[vid],attname)) for vid in vertices])
    add_position(proptoexport, 'TopPos', mg3.get_gu_top_positions(mtg))
    add_position(proptoexport, 'MidPos', mg3.get_gu_mean_positions(mtg))
    add_position(proptoexport, 'BotPos', mg3.get_gu_bottom_positions(mtg))
    pd = DataFrame(proptoexport)
    return pd