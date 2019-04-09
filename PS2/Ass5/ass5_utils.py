import numpy as np
from astropy import units as u

def getf(Omega, phi, outputUnits=None):
    """
    Gets coriolis parameter
    :param Omega:
    :param phi:
    :param outputUnits:
    :return:
    """
    f = 2 * Omega * np.sin(phi)

    if outputUnits is not None:
        f = f.to(outputUnits)

    return f

def getU(f, Deltah, L, g=9.81*u.m/u.s**2, outputUnits=None):
    """
    Get magnitude of geostrophic current speed from altimetry data
    :param f:
    :param Deltah:
    :param L:
    :param g:
    :param outputUnits:
    :return:
    """
    U = (g/f)*(Deltah/L)

    if outputUnits is not None:
        U = U.to(outputUnits, equivalencies=u.dimensionless_angles())

    return U

def getU2(hB, f, L, rhoA, rhoB, g=9.81*u.m/u.s**2, outputUnits=None):
    """
    Get magnitude of geostrophic current speed from surface data
    :param hB:
    :param f:
    :param L:
    :param rhoA:
    :param rhoB:
    :param g:
    :param outputUnits:
    :return:
    """
    U = ((g*hB)/(f*L))*(1 - rhoB/rhoA)

    if outputUnits is not None:
        U = U.to(outputUnits, equivalencies=u.dimensionless_angles())

    return U

