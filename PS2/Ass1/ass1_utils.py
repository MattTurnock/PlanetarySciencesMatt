import numpy as np
from astropy import units as u

###################################### Q1 ###############################################

def getBr(Mb, r, theta, outputUnits=u.nT):
    Br = - 2 * (Mb/r**3) * np.cos(theta)
    if outputUnits is not None:
        Br = Br.to(outputUnits)

    return Br

def getBtheta(Mb, r, theta, outputUnits=u.nT):
    Btheta = (Mb/r**3) * np.cos(theta)
    if outputUnits is not None:
        Btheta = Btheta.to(outputUnits)

    return Btheta

###################################### Q2 ###############################################
def getV(E, m, outputUnits=u.m/u.s):
    """
    Gets velocity of a particle with KE E and mass m
    :param E:
    :param m:
    :param outputUnits:
    :return:
    """
    V = np.sqrt(2*E/m)

    if outputUnits is not None:
        V = V.to(outputUnits)

    return V

def getr_lam(m, Vperp, q, B, outputUnits=u.m):
    r_lam = (m*Vperp)/(np.abs(q)*B)
    if outputUnits is not None:
        r_lam = r_lam.to(outputUnits)

    return r_lam


def getXYZ(t, Vperp, Vpar, r_lam):
    omega = Vperp/r_lam
    theta = omega*t

    z = r_lam * np.cos(theta)
    y = r_lam * np.sin(theta)

    x = Vpar * t

    return x,y,z

###################################### Q3 ###############################################

def getRhoS(rho_p, aR, rp, outputUnits=u.Unit("kg/m**3")):
    frac = rho_p**3/aR**3
    rho_s = (16*frac*rp**3)**(1/3)
    # print(rho_s.decompose())
    if outputUnits is not None:
        rho_s = rho_s.to(outputUnits)

    return rho_s

###################################### Q3 ###############################################

def getDV(mu, anom, W, outputUnits=u.m/u.s):
    Vlow = np.sqrt(mu/(anom - 0.5*W))
    Vhigh = np.sqrt(mu/(anom + 0.5*W))
    DV = Vlow - Vhigh
    if outputUnits is not None:
        DV = DV.to(outputUnits)

    return DV

def get_tSpread(anom, DV, outputUnits=u.s):
    tSpread = (2*np.pi*anom)/DV
    if outputUnits is not None:
        tSpread = tSpread.to(outputUnits)

    return tSpread


def getAring(Rout, Rin, outputUnits=u.m**2):
    Aring = np.pi*(Rout - Rin)*(Rout + Rin)
    if outputUnits is not None:
        Aring = Aring.to(outputUnits)

    return Aring

def getAring2(anom, W, outputUnits=u.m**2):
    Rout = anom + 0.5*W
    Rin = anom - 0.5*W
    Aring = getAring(Rout, Rin, outputUnits=outputUnits)

    return Aring


def get_tDiff(A1, A2, Vv, outputUnits=u.s):
    tDiff = (abs(A2-A1))/Vv
    if outputUnits is not None:
        tDiff = tDiff.to(outputUnits)

    return tDiff


