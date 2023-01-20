from openalea.plantgl.all import *
import numpy as np
from matplotlib.pyplot import *
from pandas import *
from math import ceil


def dronescan(scene, height = 600, resolution = 2):
    bbx = BoundingBox(scene)
    worldWidth = bbx.getXRange()
    worldheight = bbx.getYRange()
    up = (0,1,0)
    w, h = max(2,int(ceil(worldWidth/resolution))+1), max(2,int(ceil(worldheight/resolution))+1)
    z = ZBufferEngine(w,h,eColorBased)
    z.multithread = False
    minh = height-bbx.upperRightCorner.z
    maxh = height-bbx.lowerLeftCorner.z+minh
    z.setOrthographicCamera(-worldWidth/2., worldWidth/2., -worldheight/2., worldheight/2., minh , maxh)
    eyepos = bbx.getCenter()
    eyepos.z = height
    z.lookAt(eyepos, bbx.getCenter(), up) 
    z.process(scene)
    return z.grabZBufferPoints()

fname0 = "/Users/fboudon/Develop/Acquisition/UAVPhotogrammetry-Emma/nuages_points_drones/av_pruning_%s.txt"
fname1 = "/Users/fboudon/Develop/Acquisition/UAVPhotogrammetry-Emma/nuages_points_drones/ap_pruning_%s.txt"
fname2 = "/Users/fboudon/Develop/Acquisition/UAVPhotogrammetry-Emma/nuages_points_drones/fin_essai_%s.txt"
#s = Scene(fname)
#Viewer.display(s)

def volume_estimate(treepoints, soilheight = 0, unit = 'cm', display = False, soilpoints = None):
    scaling = { 'm' : 1, 'dm' : 1, 'cm' : 100}[unit]
    resolution = 0.2 * scaling
    grid = Point2Grid(resolution, treepoints.sliceCoordinates(0,1))
    def meanheight(grid,i,j):
        idx = grid[i,j]
        midpoint = grid.getVoxelCenter((i,j))
        if len(idx) == 0 : return 0
        mh = treepoints.subset(idx).getCenter().z -soilheight
        #mh = (mpts['Z'][idx].mean() -soilorigin.z)
        #print(mpts['Z'][idx].mean(),soilorigin.z)
        #mh = abs(soilplane.getDistance(Vector3(treepoints[idx[0]].x,treepoints[idx[0]].y,mh)))
        return mh/scaling
    heights = [[meanheight(grid,i,j) for i in range(grid.dimensions()[0])] for j in range(grid.dimensions()[1])]

    tr = Vector3(grid.getVoxelLowerPoint((0,0)),soilheight)
    vxsize = grid.getVoxelSize()/scaling
    if display:
        sc = Scene([Shape(Translated(-tr,PointSet(treepoints)),Material((0,255,0)))
                              ])
        if soilpoints : 
            sc.add(Shape(Translated(-tr,PointSet(soilpoints)),Material((255,0,0))))
        sc.add(Shape(ElevationGrid(heights, vxsize[0], vxsize[1]), Material(transparency=0.5)))
        Viewer.display(sc)

        imshow(heights)
        colorbar()
        show()

    heights = np.array(heights)
    volume = heights.sum()*vxsize[0]*vxsize[1]

    return volume

def volume(fname, simul = False, unit = 'm', display = False):
    mpts = read_csv(fname,sep=' ')
    mpts['X']-=mpts['X'].mean()
    mpts['Y']-=mpts['Y'].mean()
    mpts['Z']-=mpts['Z'].mean()
    pts = Point3Array(mpts[['X','Y','Z']].to_numpy())
    #pts = s[0].geometry.geometry.pointList
    #pts.translate(s[0].geometry.translation)
    zminidx, zmaxidx = pts.getZMinAndMaxIndex()
    zmin = pts[zminidx].z
    zmax = pts[zmaxidx].z
    zdelta = zmax-zmin
    treepoints = pts.filterCoordinates(2,zmin+zdelta*0.4,zmax)
    soilpoints = pts.filterCoordinates(2,zmin, zmin+zdelta*0.2)
    soilorigin, soilplane = Fit.plane(soilpoints)
    return volume_estimate(treepoints, soilheight = soilorigin.z, unit = 'm', display = display, soilpoints=soilpoints)


if __name__ == "__main__":
    volume(fname0 % 'C16', display = True)
#for name in ['B12','C13','C14','C16','D7','D8','E11','F7','F8','F11','H15','H16']:
#    print('volume %s :' % name, volume(fname0 % name),volume(fname1 % name),volume(fname2 % name))



