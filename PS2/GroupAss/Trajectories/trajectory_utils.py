import numpy as np
from astropy import units as u


def C3toDV(C3, muE, rOrbit, outputUnits=None):
    """
    Assumes all units are SI base if not stated
    :param C3:
    :param muE:
    :param rOrbit:
    :param outputUnits:
    :return:
    """
    # print(C3, muE, rOrbit)
    DV = np.sqrt(C3 + (2*muE/rOrbit)) - np.sqrt(muE/rOrbit)

    if outputUnits is not None:
        DV = DV.to(outputUnits)

    return DV