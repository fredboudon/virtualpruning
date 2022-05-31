header = '''
/*
 * A povray file generated with GEOM.
 * Example of use : povray -Ifile.pov -Ofile.png +FN +H600 +W800 +A.
 * File Generated with PlantGL 3D Viewer.
 */


#ifndef (__camera_definition__)
#declare __camera_definition__ = true;

camera {
   fisheye
   angle 150
    location <25,30,150>
    sky <1,0,0>
    right <1,0,0>
    look_at <25,30,151>
}

light_source {
     <0,0,0>
    color rgb 0
}



background { color rgb <1,1,1> }


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

def generate_from_representation(scene, size = 800, debug = False):
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
        povstream.write(header)
        if len(scene) > 0:
            povstream.write('#include "mangostructure.pov"\n\n\n')
        povstream.close()
        cmd = povray_exe
        if platform == "win32":
            cmd +="/EXIT /RENDER "
        else:
            cmd +="-I"
        cmd += tmpfile+".pov -O"+tmpfile+".png +H"+str(size)+" +W"+str(size)+" +FN -GA -V "
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

def generate(mtg, colorizer = mp.BlackColoring, size = 800, leaves = True, gc = True, debug = False):
    sc = mp.representation(mtg, colorizer=colorizer, leaves=leaves, gc=gc)
    return generate_from_representation(sc, size, debug)

def nbwhite(img):
    import numpy as np
    npimg = np.asarray(img)
    return npimg[np.where(npimg>0)].sum()

def gap_fraction(mtg, size = 800, debug = False):
    return gap_fraction_from_scene(mp.representation(mtg, colorizer=mp.BlackColoring, leaves=True, gc=True), size=size, debug=debug)

def gap_fraction_from_scene(scene, size = 800, debug = False):
    from openalea.plantgl.all import Scene
    return nbwhite(generate_from_representation(scene, size=size, debug=debug))/nbwhite(generate_from_representation(Scene(), size=size, debug=debug))
