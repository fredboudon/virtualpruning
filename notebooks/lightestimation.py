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


def get_light_sources(diffuseratio = 0.3, energy = 100000, date='2017-08-26', starthour = 7, endhour = 18):
    def todate(hour = 12, date=date):
        return pandas.Timestamp(date+' '+str(hour)+':00', tz=localisation['timezone'])
    hours = pandas.date_range(start=todate(starthour),end=todate(endhour), freq="1H")
    suns = sun_sources(energy*(1-diffuseratio), dates=hours, **localisation)
    skys = sky_sources(sky_type='uoc', irradiance=energy*diffuseratio, **localisation)
    sun_el, sun_az, sun_hei = suns
    sky_el, sky_az , sky_hei = skys
    return (sun_az, sun_el, sun_hei), (sky_az, sky_el,sky_hei)

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
    scene.setLight(light)
    if debug:
        print('Run caribu')
    #raw, agg = scene.run(direct=False, infinite = True, split_face = True, d_sphere = D_SPHERE)
    raw, agg = scene.run(direct=True, infinite = False, split_face = True)
    if debug:
        print('made in', time.time() - t)
    agg = agg['PAR']
    agg['irradiance'] = agg['Ei']/energy
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

def set_light_to_mtg(mtg, lightprop):
    for propname, propvalues in lightprop.items():
        mtg.property(propname).update(propvalues)


def zeta(TrPPFD):
    from numpy import exp, power
    a = 10.54
    b = 1.596
    d = 1.104
    return d * power((1-exp(-a*TrPPFD))/(1-exp(-a)), 1/b)

def daily_light_estimation(scene, diffuseratio = 0.3, energy = 100000, date='2017-08-26', starthour = 7, endhour = 18, debug = True):
    from pandas import concat
    sun, sky = get_light_sources(diffuseratio=diffuseratio, energy = energy, date=date, starthour = starthour, endhour = endhour)
    skyvalues = light(scene, None, sky, debug = debug)['irradiance']
    values = {}
    hour = starthour
    for sun_az, sun_el, sun_hei in zip(*sun):
        values[str(hour)+'H'] = light(scene, ([sun_az], [sun_el], [sun_hei]), None, debug = debug)['irradiance']+skyvalues
        hour += 1
    res = concat(values, axis=1)
    res.fillna(0)
    return res 

def light_variables(scene, diffuseratio = 0.3, date='2017-08-26', starthour = 7, endhour = 18, debug = False):
    # Compute 
    # TrPPFD_min : Mesure ponctuelle min sur la journée
    # Zeta_min : Mesure ponctuelle min sur la journée
    # Zeta_12H : Mesure ponctuelle a 12h
    energy = 100000
    lightestim = daily_light_estimation(scene, 
                            diffuseratio = diffuseratio, 
                            energy = energy, 
                            date=date, 
                            starthour = starthour, 
                            endhour = endhour,
                            debug = debug)/energy
    TrPPFD_min = lightestim.min(axis=1)
    zetavalues = zeta(lightestim)
    Zeta_min = zetavalues.min(axis=1)
    Zeta_12H = zetavalues['12H']
    return TrPPFD_min, Zeta_min, Zeta_12H
