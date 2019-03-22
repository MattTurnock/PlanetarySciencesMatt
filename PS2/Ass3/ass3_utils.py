import numpy as np
from astropy import units as u
from json_to_dict import constants

def getPropTime(v, d=constants["AU"], outputUnits=u.s):
    """
    Assumes standard SI units if constants is NOne
    :param v:
    :param d:
    :param outputUnits:
    :return:
    """
    t = d/v
    if constants is not None:
        t = t.to(outputUnits)

    return t

def getLorentzFactor(beta):
    gamma = 1/(np.sqrt(1-beta**2))
    return gamma