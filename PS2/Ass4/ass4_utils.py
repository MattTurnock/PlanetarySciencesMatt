import numpy as np
from astropy import units as u

def getn(mu, a, outputUnits=u.rad/u.s):
    n = np.sqrt(mu/a**3)
    if outputUnits is not None:
        n = n.to(outputUnits, equivalencies=u.dimensionless_angles())
    return n

def getp(a, e):
    p = a * (1-e**2)
    return p

def getomegaDot(n, J2, Re, p, i, outputUnits=u.rad/u.s ):
    """
    Rate of change of AoP
    :param n: mean motion
    :param J2: J2 coefficient
    :param Re: Earth radius
    :param p: thta orbital thing a(1-e**2)
    :param i: inlcination
    :param outputUnits:
    :return:
    """
    top = 3 * n * J2 * Re**2
    bottom = 2 * p**2
    side = 2 - (5/2) * (np.sin(i))**2

    omegaDot = (top / bottom) * side

    if outputUnits is not None:
        omegaDot = omegaDot.to(outputUnits)

    return omegaDot