header = '''
/*
 * A povray file generated with GEOM.
 * Example of use : povray -Ifile.pov -Ofile.png +FN +H600 +W800 +A.
 * File Generated with PlantGL 3D Viewer.
 */
#ifndef (__camera_definition__)
#declare __camera_definition__ = true;
camera {{
   fisheye
   angle 150
    location <{xpos},{ypos},{camheight}>
    sky <1,0,0>
    right <1,0,0>
    look_at <{xpos},{ypos},{camheight} + 1>
}}
light_source {{
     <0,0,0>
    color rgb 0
}}
background {{ color rgb <1,1,1> }}
#end // __camera_definition__
'''


import mtgplot as mp
from importlib import reload
reload(mp)

povray_exe = None

def read_povray_exe(fname = 'config.cfg'):
    from sys import platform
    import os
    if os.path.exists(fname):
        return open(fname,'r').read().strip()
    else:
        if platform == "linux" or platform == "linux2" or platform == 'darwin':
            return '/opt/local/bin/povray'
        elif platform == 'win32':
            return 'C:/Program Files/povray/v3.7/bin/pvengine.exe'
        else:
            return 'povray'

def get_povray_exe():
    global povray_exe
    if povray_exe is None:
        povray_exe = read_povray_exe()
    return povray_exe

refcamposition =  (-12, 22, 122) # (25,30,120) # 

def pov_hemispherical_view_from_scene(scene, size = 400, camposition = refcamposition, antialiasing = False, debug = False):
    from sys import platform
    import numpy as np
    import pandas as pd
    from math import degrees
    import os
    import shutil
    xpos, ypos, camheight = camposition
    povray_exe = get_povray_exe()
    if os.path.exists('povray'):
        shutil.rmtree('povray')
    os.makedirs('povray')
    os.chdir('povray')
    if debug:
        print(os.getcwd())
    try:
        tmpfile = 'gapfraction_view'
        if len(scene) > 0:
            scene.save('mangostructure.pov')
        povstream = open(tmpfile+'.pov','w')
        povstream.write(header.format(camheight=camheight, xpos=xpos, ypos=ypos))
        if len(scene) > 0:
            povstream.write('#include "mangostructure.pov"\n\n\n')
        povstream.close()
        cmd = povray_exe
        if platform == "win32":
            cmd +=" /EXIT /RENDER "
        else:
            cmd +=" -I"
        cmd += tmpfile+".pov -O"+tmpfile+".png +H"+str(size)+" +W"+str(size)+" +FN -GA "
        if antialiasing:
            cmd += '+A'
        else:
            cmd += '-A'            
        if platform == "linux" or platform == "linux2" or platform == 'darwin':
            cmd += " &> /dev/null"
        from time import time
        t = time()
        if debug:
            print(cmd)
        os.system(cmd)
        if debug:
            print('Done in',time()-t, 'sec.')
        from PIL import Image
        img = np.asarray(Image.open(tmpfile+".png"))
    except Exception as e:
        os.chdir(os.pardir)
        raise e        
    os.chdir(os.pardir)
    if not debug:
        shutil.rmtree('povray')
    return img

def nbwhite(img):
    import numpy as np
    npimg = np.asarray(img)
    # nb of pix coordinates divide the number of channels
    return len(np.where(npimg>0)[0])//img.shape[2] 

def pov_gap_fraction_from_scene(scene, size = 400, camposition = refcamposition, antialiasing = False, debug = False):
    from openalea.plantgl.all import Scene
    
    img = pov_hemispherical_view_from_scene(scene, size=size, camposition=camposition, antialiasing = antialiasing, debug=debug)
    if debug:
        from matplotlib.pyplot import imshow, show
        imshow(img)
        show()
    nbpix = nbwhite(img)
    if size == 800:
        refwhite = 502652 # 384528780
    elif size == 400:
        refwhite = 125676 #c96142140
    else:
        img = pov_hemispherical_view_from_scene(Scene(), size=size, debug=debug)
        refwhite = nbwhite(img)
    if debug:
        print(nbpix,refwhite)
    return nbpix/refwhite

########################################################

def _setup_pgl_engine(size = 400, camposition = refcamposition, antialiasing = False):
    import openalea.plantgl.all as pgl
    defaultcolor = pgl.Color3(255,255,255)
    engine = pgl.ZBufferEngine(size,size,pgl.eColorBased,defaultcolor)
    engine.setSphericalCamera(150)
    engine.setLight((0,0,0),(0,0,0))
    pos = pgl.Vector3(camposition)
    engine.lookAt(pos,pos+pgl.Vector3(0,0,1),(1,0,0))
    return engine

def pgl_hemispherical_view_from_scene(scene, size = 400, camposition = refcamposition, antialiasing = False, debug = False):
    engine = _setup_pgl_engine(size, camposition, antialiasing)
    engine.process(scene)
    return engine.getImage().to_array()

def pgl_nbwhite(img, size):
    return ((size*size) - dict(img.histogram())[0])

def pgl_gap_fraction_from_scene(scene, size = 400, camposition = refcamposition, antialiasing = False, debug = False):
    engine = _setup_pgl_engine(size, camposition, antialiasing)
    refwhite = pgl_nbwhite(engine.getImage(), size)
    engine.process(scene)
    nbpix = pgl_nbwhite(engine.getImage(), size)
    if debug:
       from matplotlib.pyplot import imshow, show
       imshow(engine.getImage().to_interlaced_array())
       show()
    if debug:
        print(nbpix,refwhite)
    return nbpix/refwhite


def pgl_gapfraction_dynamic(mtg, size = 400, camposition = refcamposition, antialiasing = False, debug = False):
    import numpy as np
    import openalea.plantgl.all as pgl
    dates = sorted(np.unique(list(mtg.property('BurstDate').values())))
    sc = generate_representation(mtg, colorizer=mp.BlackColoring, leaves=True, gc=False, todate=dates[-1])
    scd = sc.todict()
    def vertex_selection(g, date=None):
        return [vid for vid in mtg.vertices(scale=3) if g.property('BurstDate').get(vid, date) == date]
    scene = pgl.Scene(sum([scd[vid] for vid in vertex_selection(mtg) if vid in scd],[]))
    engine = _setup_pgl_engine(size, camposition, antialiasing)
    refwhite = pgl_nbwhite(engine.getImage(), size)
    engine.process(scene)
    nbpix = pgl_nbwhite(engine.getImage(), size)
    result = [(None, nbpix/refwhite)]
    for d in dates:
        scene = pgl.Scene(sum([scd[vid] for vid in vertex_selection(mtg,d) if vid in scd],[]))
        engine.process(scene)
        nbpix = pgl_nbwhite(engine.getImage(), size)
        result += [(d, nbpix/refwhite)]
    return result

########################################################

PglRendering, PovRendering = 1,2

def hemispherical_view_from_scene(scene, size = 400, camposition = refcamposition, antialiasing = False, debug = False, method = PglRendering):
    funcs = { PglRendering : pgl_hemispherical_view_from_scene , PovRendering : pov_hemispherical_view_from_scene }
    return funcs[method](scene, size=size, camposition=camposition, antialiasing = antialiasing, debug=debug)

def generate_representation(mtg, colorizer = mp.BlackColoring, leaves = True, wood = True, gc = True, debug = False, todate=None):
    return mp.representation(mtg, colorizer=colorizer, leaves=leaves, wood = wood, gc=gc, todate=todate)

def hemispherical_view(mtg, colorizer = mp.BlackColoring, size = 400, camposition = refcamposition, antialiasing = False, leaves = True, wood = True, gc = True, debug = False, todate=None, method = PglRendering):
    sc = generate_representation(mtg, colorizer=colorizer, leaves=leaves, wood = wood, gc=gc, todate=todate)
    return hemispherical_view_from_scene(sc, size, camposition=camposition, antialiasing = antialiasing, debug=debug, method=method)

def gap_fraction_from_image(img):
    size = img.shape[0]
    nbpix = nbwhite(img)
    import openalea.plantgl.all as pgl
    defaultcolor = pgl.Color3(255,255,255)
    engine = pgl.ZBufferEngine(size,size,pgl.eColorBased,defaultcolor)
    engine.setSphericalCamera(150)
    refwhite = pgl_nbwhite(engine.getImage(), size)
    return nbpix/refwhite

def gap_fraction_from_scene(scene, size = 400, 
                                   camposition = refcamposition, 
                                   antialiasing = False, 
                                   debug = False,
                                   method = PglRendering):
    funcs = { PglRendering : pgl_gap_fraction_from_scene , PovRendering : pov_gap_fraction_from_scene }
    return funcs[method](scene, size=size, camposition=camposition, antialiasing = antialiasing, debug=debug)


def gap_fraction(mtg, size = 400, 
                      camposition = refcamposition, 
                      antialiasing = False, 
                      debug = False, 
                      todate = None,
                      method = PglRendering):
    
    return gap_fraction_from_scene(generate_representation(mtg, colorizer=mp.BlackColoring, leaves=True, gc=True, todate=todate), size=size, camposition=camposition, antialiasing = antialiasing, debug=debug, method=method)


def gapfraction_dynamic(mtg, size = 400, camposition = refcamposition, antialiasing = False, debug = False, method = PglRendering):
    import numpy as np
    dates = sorted(np.unique(list(mtg.property('BurstDate').values())))

    return [(d,gap_fraction(mtg, size=size, camposition=camposition, antialiasing = antialiasing, debug=debug, todate=d, method=method)) for d in dates]



