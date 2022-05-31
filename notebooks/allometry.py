import numpy as np
from math import *


def gu_leaf_area(diameter):
    """ Allometric relationship between diameter of a GU and the carried leaf area
     :param diameter: diameter of the GU in mm2
     :return: leaf area in dm2
    """
    return exp(1.093*log(pi*(diameter**2)/4.)-2.146)

def gu_biomass(diameter):
    """ Allometric relationship between diameter of a GU and the carried biomass
     :param diameter: diameter of the GU in mm2
     :return: biomass in g
    """
    return exp(1.2576*log(pi*(diameter**2)/4.)-0.9826)