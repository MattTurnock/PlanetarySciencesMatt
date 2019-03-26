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

def getOmegaDot(n, J2, Re, p, i, outputUnits=u.rad/u.s ):
    """
    Rate of change of RAAN
    :param n: mean motion
    :param J2: J2 coefficient
    :param Re: Earth radius
    :param p: thta orbital thing a(1-e**2)
    :param i: inlcination
    :param outputUnits:
    :return:
    """

    OmegaDot = -(3 * n * J2 * Re**2 * np.cos(i)) / (2 * p**2)

    if outputUnits is not None:
        OmegaDot = OmegaDot.to(outputUnits)

    return OmegaDot


def getOmegaDot2(J2, Re, p, i, T, outputUnits=u.rad/u.s ):
    """
    Like OmegaDot function, but uses period T instead of n
    :param J2:
    :param Re:
    :param p:
    :param i:
    :param T:
    :param outputUnits:
    :return:
    """
    OmegaDot = -3*np.pi*J2* (Re/p)**2 *np.cos(i) * (1/T)

    if outputUnits is not None:
        OmegaDot = OmegaDot.to(outputUnits)

    return OmegaDot

def getiFromOmegaDot(OmegaDot, T, J2, p, Re, outputUnits=None):

    cosi = -((T*OmegaDot)/(3*np.pi*J2) * (p/Re)**2)

    i = np.arccos(cosi)
    if outputUnits is not None:
        i = i.to(outputUnits)

    return i


def getfa(fs, f, outputUnits=None):
    """
    Calculates aliasing frequency fa
    :param fs:
    :param f:
    :param outputUnits:
    :return:
    """
    n = int(np.round(f / fs))
    print(n)

    fa = abs(n * fs - f)

    if outputUnits is not None:
        fa = fa.to(outputUnits)

    return fa