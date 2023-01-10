from openalea.mtg import *
from math import pi
from openalea.plantgl.all import Vector3, norm, direction

def get_G3_mtg():
    g = MTG("data/consolidated_mango3d.mtg")

    def pos_prop(mtg):
        xx,yy,zz = mtg.property("XX"),mtg.property("YY"),mtg.property("ZZ")
        return dict ([(node,Vector3(x,-yy[node],-zz[node])) for node,x in xx.items()])

    g.property('Position').update(pos_prop(g))

    def diameter_prop(mtg):
        topdia = lambda x: mtg.property('Diameter').get(x)
        diameters = {}
        for vid in mtg.vertices(scale=3):
            td = topdia(vid)
            if td is None :
                if not mtg.parent(vid) is None:
                    diameters[vid] = diameters[mtg.parent(vid)]
            else:
                diameters[vid] = td
        return diameters

    g.property('Diameter').update(diameter_prop(g))
    g.property('Taille')[6817] = 'A'
    g.property('Taille')[15586] = 'A'
    

    return g

GUScale = 3

def is_gu_point(mtg, vid):
    """ A GU Point is a node in the graph that represent both a digitized point and a GU.
        Conversely, some point that mark the begining of a branch do not represent any GU."""
    return mtg.edge_type(vid) == '<'

def get_gu_bottom_point(mtg,vid):
    return mtg.parent(vid)

def get_gu_top_point(mtg,vid):
    return vid

def is_terminal(mtg,vid):
    return mtg.nb_children(vid) == 0

def get_all_gus(mtg):
    return [vid for vid in mtg.vertices(scale=GUScale) if is_gu_point(mtg, vid)]

def get_ordered_gus(mtg, preorder = True, root = None):
    from openalea.mtg.traversal import pre_order2, post_order2
    orderfunc = pre_order2 if preorder else post_order2
    if root is None:
        root = mtg.roots(scale=GUScale)[0]
    return [root]+[vid for vid in orderfunc(mtg, root) if is_gu_point(mtg, vid)]

def get_all_terminal_gus(mtg):
    return [vid for vid in mtg.vertices(scale=3) if is_gu_point(mtg, vid) and is_terminal(mtg,vid) ]

def get_terminal_gus_from_ancestor(mtg, vid):
    return [vid for vid in mtg.Extremities(vid) if is_gu_point(mtg, vid)]

def get_descendants_gus_from_ancestor(mtg, vid):
    return [vid for vid in mtg.Descendants(vid) if is_gu_point(mtg, vid)]

def get_gu_diameter(mtg, vid):
    """ Diameter is stored in the bottom point.
        :return: diameter in mm  """
    value = mtg.property('Diameter').get(get_gu_bottom_point(mtg,vid))
    if value is None : return value
    return value*10

def set_gu_diameter(mtg, vid, value):
    """ Set diameter in the bottom point. """
    mtg.property('Diameter')[get_gu_bottom_point(mtg,vid)] = value/10.

def get_gu_section(mtg, vid):
    """ Return section value of the GU in mm2 """
    return pi*(get_gu_diameter(mtg, vid)**2)/4

def get_gu_nb_leaf(mtg, vid):
    """ NbLeaf is stored in the top point """
    return mtg.property('NbLeaf').get(get_gu_top_point(mtg,vid),0)

def get_gu_type(mtg, vid):
    """ NbLeaf is stored in the bottom point """
    return mtg.property('UnitType').get(get_gu_bottom_point(mtg,vid))

def get_gu_bottom_position(mtg, vid):
    return mtg.property('Position')[get_gu_bottom_point(mtg,vid)]

def get_gu_top_position(mtg, vid):
    return mtg.property('Position')[get_gu_top_point(mtg,vid)]

def get_gu_mean_position(mtg, vid):
    return (mtg.property('Position')[get_gu_top_point(mtg,vid)]+mtg.property('Position')[get_gu_bottom_point(mtg,vid)])/2

def set_gu_top_position(mtg, vid, value):
    mtg.property('Position')[get_gu_top_point(mtg,vid)] = value

def get_gu_direction(mtg, vid):
    """ Gives direction of the GU. Length of the returned vector represent the lenght of the GU """
    return (get_gu_top_position(mtg, vid) - get_gu_bottom_position(mtg, vid))

def get_gu_normed_direction(mtg, vid):
    """ Gives direction of the GU. Length of the returned vector represent the lenght of the GU """
    return direction(get_gu_direction(mtg, vid))

def get_gu_length(mtg, vid):
    return norm(get_gu_direction(mtg, vid))

def get_gu_depth(mtg, vid1, vid2):
    assert is_gu_point(mtg, vid1) and is_gu_point(mtg, vid2)
    return len([vid for vid in mtg.Path(vid1,vid2) if is_gu_point(mtg, vid)])-1

def get_gu_terminal_min_depth(mtg, vid):
    if is_terminal(mtg, vid) : return 0
    return min([get_gu_depth(mtg, vid, desc) for desc in get_terminal_gus_from_ancestor(mtg, vid)])

def get_gu_property(mtg, vid, propname, toppoint = True):
    return mtg.property(propname)[get_gu_top_point(mtg,vid) if toppoint else get_gu_bottom_point(mtg,vid)]

def set_gu_property(mtg, vid, propname, value, toppoint = True):
    mtg.property(propname)[get_gu_top_point(mtg,vid) if toppoint else get_gu_bottom_point(mtg,vid)] = value

def get_gu_top_positions(mtg):
    return dict([(vid,get_gu_top_position(mtg,vid)) for vid in get_all_gus(mtg)])

def get_gu_bottom_positions(mtg):
    return dict([(vid,get_gu_bottom_position(mtg,vid)) for vid in get_all_gus(mtg)])

def get_gu_mean_positions(mtg):
    return dict([(vid,get_gu_mean_position(mtg,vid)) for vid in get_all_gus(mtg)])

def was_previously_pruned(mtg, vid):
    return ('A' in mtg.property('Taille').get(vid,''))

def get_parent(mtg, vid):
    assert not vid is None
    vid = mtg.parent(vid)
    while not vid is None and not is_gu_point(mtg, vid):
        vid = mtg.parent(vid)
    assert vid is None or is_gu_point(mtg, vid)
    return vid


def get_children(mtg, vid):
    assert not vid is None
    fchildren = []
    for ch in mtg.children(vid):
        if not is_gu_point(mtg, ch):
            fchildren += get_children(mtg, ch)
        else:
            fchildren.append(ch)
    return fchildren

def get_ancestor(mtg, vid, order):
    assert order > 0
    for j in range(order):
        vid = get_parent(mtg, vid)
    return vid

eApical, eLateral = 1,2

def gu_position(mtg, vid):
    return eApical if mtg.edge_type(get_gu_bottom_point(mtg,vid)) == '<' else eLateral


def nbtotalleaves(mtg):
    return gu_recursive_property_from_terminal(mtg, lambda vid, childrenvalues : sum(childrenvalues)+get_gu_nb_leaf(mtg, vid), lambda vid: get_gu_nb_leaf(mtg, vid))

def gus_depth_from_terminal(mtg, agregation = min, terminalorder = 0):
    return gu_recursive_property_from_terminal(mtg, lambda vid, childrenvalues : agregation(childrenvalues)+1, terminalorder)

def gus_depth_from_root(mtg):
    return gu_recursive_property_from_root(mtg)

def gu_recursive_property_from_terminal(mtg, nodeaxiom = lambda vid, childrenvalues : sum(childrenvalues)+1, leafaxiom = lambda vid : 0, root = None):
    from openalea.mtg.traversal import post_order2
    res = {}
    if root is None:
        root = mtg.roots(scale=GUScale)[0]
    for vid in post_order2(mtg, root):
        if is_terminal(mtg, vid): 
            res[vid] = leafaxiom(vid) if callable(leafaxiom) else leafaxiom
        elif is_gu_point(mtg,vid):
            children = get_children(mtg, vid)
            for child in children:
                assert is_gu_point(mtg,child)
            res[vid] = nodeaxiom(vid, [res[child] for child in children])
    return res

def gu_recursive_property_from_root(mtg, nodeaxiom = lambda vid, parentvalue : parentvalue+1, rootaxiom = lambda vid : 0, root = None):
    from openalea.mtg.traversal import pre_order2
    res = {}
    if root is None:
        root = mtg.roots(scale=GUScale)[0]
    for vid in pre_order2(mtg, root):
        if get_parent(mtg, vid) is None: 
            res[vid] = rootaxiom(vid) if callable(rootaxiom) else rootaxiom
        elif is_gu_point(mtg,vid):
            parent = get_parent(mtg, vid)
            res[vid] = nodeaxiom(vid, res[parent])
    return res

def repare_mango_lighted():
    import openalea.plantgl.all as pgl
    s = pgl.Scene('../data/lightedG3.bgeom')
    s2 = pgl.Scene('../data/consolidated_mango3d.bgeom')

    pid2sid = {}
    points = []
    #key = lambda bbx : tuple([round(v,2) for v in bbx.getCenter()])
    key = lambda bbx :  bbx.getCenter()

    for i, (sid, lsh) in enumerate(s2.todict().items()):
        bbx = pgl.BoundingBox(pgl.Scene(lsh))
        c = key(bbx)
        points.append(c)
        pid2sid[i] = sid

    kdtree = pgl.ANNKDTree3(points)

    for i,sh in enumerate(s):
        bbx = pgl.BoundingBox(sh)
        res = kdtree.k_closest_points(key(bbx),1,pgl.norm(bbx.getSize()))
        sh.id = pid2sid[res[0]]

    pgl.Viewer.display(s)

selectionratio = 0.01 # 0.13

def bbox(mtg, selectionratio = selectionratio):
    """
    Compute the bounding box of the tree.
    :return: the dimension in cm3
    """
    import numpy as np
    def minmaxvalue(positions, coord, selectionratio):
        cmax = max(positions[:,coord])
        cmin = min(positions[:,coord])
        clayer = (cmax-cmin) * selectionratio
        downlayer = positions[positions[:,coord] < cmin+clayer,coord]
        mmin = np.mean(downlayer)
        toplayer = positions[positions[:,coord] > cmax-clayer,coord]
        mmax = np.mean(toplayer)
        return mmin, mmax
    positions = np.array([pos for vid, pos in mtg.property('Position').items() if mtg.property('UnitType')[vid] != 'B'])
    xmin, xmax = minmaxvalue(positions, 0, selectionratio)
    ymin, ymax = minmaxvalue(positions, 1, selectionratio)
    zmin, zmax = minmaxvalue(positions, 2, selectionratio)
    return (xmin, xmax), (ymin, ymax), (zmin, zmax)

def volume(mtg, selectionratio = selectionratio):
    """
    Compute the volume tree.
    :return: the volume in m3
    """
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bbox(mtg, selectionratio)
    return (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/ 1000000.

ref_volume = 67.3 # 37.661114634371984

def extend_mtg_with_organs(mtg):
    from copy import deepcopy
    mtg = deepcopy(mtg)
    for vid in get_ordered_gus(mtg):
        nvid = mtg.add_component(vid, label='W', edge_type=mtg.edge_type(vid))
        for lid in range(get_gu_nb_leaf(mtg, vid)):
            mtg.add_child(nvid,label='L', edge_type='+')
    return mtg
