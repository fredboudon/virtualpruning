localisation={'latitude':-21.32, 'longitude':55.5, 'timezone': 'Indian/Reunion'}

#leaf_prop = { 'Rc' : (0.05, 0.007, 0.078, 0.007), 
#              'Rs' : (0.413, 0.36, 0.455, 0.353),
#              'PAR' : (0.067, 0.025, 0.108, 0.023) }
#wood_prop = { 'Rc' : (0.0001, 0.0001), 'Rs' : (0.0001, 0.0001), 'PAR' : (0.0001, 0.0001)}

leaf_prop = { 'PAR' : (0.067, 0.025, 0.108, 0.023) }
wood_prop = { 'PAR' : (0.0001, 0.0001)}

idshift = 1000

xcenter, ycenter = -21,   -40
xsize,   ysize   = 250+7, 300
xmin, xmax = -xcenter-xsize, -xcenter+xsize
ymin, ymax = -ycenter-ysize,-ycenter+ysize
pattern=(xmin,ymin,xmax,ymax)

import pandas
from alinea.astk.sun_and_sky import sky_sources, sun_sources


def get_light_sources(diffuseratio = 0.3, energy = 100000, date='2021-03-01', starthour = 7, endhour = 18, skydiscretization=46):
    from math import sin, radians
    def todate(hour = 12, date=date):
        return pandas.Timestamp(date+' '+str(hour)+':00', tz=localisation['timezone'])
    hours = pandas.date_range(start=todate(starthour),end=todate(endhour), freq="1H")
    suns = sun_sources(energy*(1-diffuseratio), dates=hours, **localisation)
    skys = sky_sources(sky_type='uoc', irradiance=energy*diffuseratio, turtle_sectors=skydiscretization, **localisation)
    sun_el, sun_az, sun_hei = suns
    sky_el, sky_az , sky_hei = skys
    # every sun position should have the direct energy value
    return ((sun_az, sun_el, [energy*(1-diffuseratio) for hei,el in zip(sun_hei,sun_el)]),  # * sin(radians(el))
            (sky_az, sky_el,sky_hei))

def toCaribuScene(mangoscene, leaf_prop=leaf_prop, wood_prop=wood_prop, idshift=idshift, pattern=pattern, debug = True) :
    from alinea.caribu.CaribuScene import CaribuScene
    import time
    if debug:
        print ('Convert scene for caribu')
        t = time.time()
    geomdict = set([sh.id for sh in mangoscene])
    wavelenghts = list(leaf_prop.keys())
    opt = dict([(k,{}) for k in wavelenghts])
    for vid in geomdict:
        for rv in wavelenghts:
            opt[rv][vid] = (wood_prop if (vid % idshift) == 0 else leaf_prop)[rv]
    cs = CaribuScene(mangoscene, opt=opt, scene_unit='cm', pattern=pattern, debug = False)
    if debug:
        print('done in', time.time() - t)
    return cs


def caribu(scene, sun = None, sky = None, debug = True):
    from alinea.caribu.light import light_sources
    from openalea.plantgl.scenegraph import Scene

    import time
    if debug:
        print('start caribu...')
        t = time.time()
        print('Create light source', end=' ')
    light = []
    if not sun is None:
        sun_az, sun_el, sun_hei = sun 
        light += light_sources(sun_el, sun_az, sun_hei) #, orientation = north) 
    if not sky is None:
        sky_az, sky_el,sky_hei = sky
        light += light_sources(sky_el, sky_az, sky_hei) #, orientation = north)
    if debug:
        print('... ',len(light),' sources.')
    if isinstance(scene, Scene):
        scene = toCaribuScene(scene,debug = debug)
    scene.setLight(light)
    if debug:
        print('Run caribu')
    #raw, agg = scene.run(direct=False, infinite = True, split_face = True, d_sphere = D_SPHERE)
    raw, agg = scene.run(direct=True, infinite = False, split_face = True)
    if debug:
        print('made in', time.time() - t)
    agg = agg['PAR' if 'PAR'in agg else 'default_band']
    agg['irradiance'] = agg['Ei']
    del agg['Ei']
    import pandas as pd
    return pd.DataFrame(agg)

def plantgllight(scene, sun = None, sky = None, debug = True):
    from openalea.plantgl.light.light import scene_irradiance
    import time
    if debug:
        print('start plantgl light...')
        t = time.time()
        print('Create light source', end=' ')
    lights = []
    if sun: 
        lights += list(zip(*sun)) 
    if sky:
        lights += list(zip(*sky))
    if debug:
        print('... ',len(lights),' sources.')
        print('Run plantGL')
    agg = scene_irradiance(scene, lights, horizontal=True, screenresolution= 0.5, scene_unit = 'cm', verbose=False)
    if debug:
        print('made in', time.time() - t)
    return agg

def light(scene, sun = None, sky = None, useplantgl = True, debug = True):
    if useplantgl:
        return plantgllight(scene, sun, sky, debug)
    else:
        return caribu(scene, sun, sky, debug)

def extend_mtg_with_light(mtg, TrPPFD = None, Zeta = None, inplace=True, **properties):
    if not inplace:
        from copy import deepcopy
        mtg = deepcopy(mtg)
    if not TrPPFD is None:
        for propname, propvalues in TrPPFD.items():
            mtg.property('TrPPFD'+propname).update(propvalues)
    if not Zeta is None:
        for propname, propvalues in Zeta.items():
            mtg.property('Zeta'+propname).update(propvalues)
    for proptablename, proptablevalues in properties.items():
        for propname, propvalues in proptablevalues.items():
            mtg.property(proptablename+propname).update(propvalues)


def zeta(TrPPFD):
    from numpy import exp, power
    a = 10.54
    b = 1.596
    d = 1.104
    return d * power((1-exp(-a*TrPPFD))/(1-exp(-a)), 1/b)

def check_scene(scene):
    from openalea.mtg import MTG
    if isinstance(scene, MTG):
        import mtgplot as mp
        return mp.representation(scene, wood = False, leaves=True)
    else:
        return scene

def daily_light_estimation(scene, diffuseratio = 0.3, energy = 100000, date='2021-03-01', starthour = 7, endhour = 18, skydiscretization = 46, debug = True):
    from pandas import concat
    scene = check_scene(scene)
    sun, sky = get_light_sources(diffuseratio=diffuseratio, energy = energy, date=date, starthour = starthour, endhour = endhour, skydiscretization=skydiscretization)
    if sum(sky[2]) > 0:
        skyvalues = light(scene, None, sky, debug = debug)['irradiance']
    else:
        skyvalues = None
    values = {}
    hour = starthour
    for sun_az, sun_el, sun_hei in zip(*sun):
        if sun_hei > 0:
            value = light(scene, ([sun_az], [sun_el], [sun_hei]), None, debug = debug)['irradiance']
            if not skyvalues is None:
                value += skyvalues
        else:
            value = skyvalues
        values[str(hour)+'H'] = value
        hour += 1
    res = concat(values, axis=1)
    res.fillna(0,inplace=True)
    return res 

def daily_light_variables(scene, diffuseratio = 0.3, date='2021-03-01', starthour = 7, endhour = 18, skydiscretization = 46, debug = False):
    energy = 100000
    TrPPFD = daily_light_estimation(scene, 
                            diffuseratio = diffuseratio, 
                            energy = energy, 
                            date=date, 
                            starthour = starthour, 
                            endhour = endhour,
                            skydiscretization=skydiscretization,
                            debug = debug)/energy
    Zeta = zeta(TrPPFD)
    return TrPPFD, Zeta

def light_variables(scene, diffuseratio = 0.3, date='2021-03-01', starthour = 7, endhour = 18, skydiscretization = 46, debug = False):
    # Compute 
    # TrPPFD_min : Mesure ponctuelle min sur la journée
    # Zeta_min : Mesure ponctuelle min sur la journée
    # Zeta_12H : Mesure ponctuelle a 12h
    TrPPFD, Zeta = daily_light_variables(scene, 
                            diffuseratio = diffuseratio, 
                            date=date, 
                            starthour = starthour, 
                            endhour = endhour,
                            skydiscretization=skydiscretization,
                            debug = debug)
    Zeta_mean = Zeta.mean(axis=1)
    return Zeta_mean




