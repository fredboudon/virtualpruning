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

def generate_from_representation(scene, size = 400, camheight = 120, xpos = 25, ypos = 30, antialiasing = False, debug = False):
    from sys import platform
    import numpy as np
    import pandas as pd
    from math import degrees
    import os
    import shutil
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

def generate(mtg, colorizer = mp.BlackColoring, size = 400, camheight = 150, xpos = 25, ypos = 30, antialiasing = False, leaves = True, wood = True, gc = True, debug = False, todate=None):
    sc = mp.representation(mtg, colorizer=colorizer, leaves=leaves, wood = wood, gc=gc, todate=todate)
    return generate_from_representation(sc, size, camheight=camheight, xpos = xpos, ypos = ypos, antialiasing = antialiasing, debug=debug)

def nbwhite(img):
    import numpy as np
    npimg = np.asarray(img)
    return npimg[np.where(npimg>0)].sum()

def gap_fraction(mtg, size = 400, camheight = 150, xpos = 25, ypos = 30, antialiasing = False, debug = False, todate=None):
    return gap_fraction_from_scene(mp.representation(mtg, colorizer=mp.BlackColoring, leaves=True, gc=True, todate=todate), size=size, camheight=camheight, xpos = xpos, ypos = ypos, antialiasing = antialiasing, debug=debug)

def gap_fraction_from_scene(scene, size = 400, camheight = 150, xpos = 25, ypos = 30, antialiasing = False, debug = False):
    from openalea.plantgl.all import Scene
    nbpix = nbwhite(generate_from_representation(scene, size=size, camheight=camheight, xpos = xpos, ypos = ypos, antialiasing = antialiasing, debug=debug))
    if size == 800:
        refwhite = 384528780
    elif size == 400:
        refwhite = 96142140
    else:
        img = generate_from_representation(Scene(), size=size, debug=debug)
        if debug:
            from matplotlib.pyplot import imshow, show
            imshow(img)
            show()
        refwhite = nbwhite(img)
    if debug:
        print(nbpix,refwhite)
    return nbpix/refwhite

def gapfraction_dynamic(mtg, size = 400, camheight = 150, xpos = 25, ypos = 30, antialiasing = False, debug = False):
    import numpy as np
    dates = sorted(np.unique(list(mtg.property('BurstDate').values())))
    return [(d,gap_fraction(mtg, size=size, camheight=camheight, xpos = xpos, ypos = ypos, antialiasing = antialiasing, debug=debug, todate=d)) for d in dates]