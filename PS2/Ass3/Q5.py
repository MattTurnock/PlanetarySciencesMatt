from PS2.Ass1.ass1_utils import *
from PS2.Ass3.ass3_utils import *
from json_to_dict import constants

from astropy import units as u

def getSigma(Ai, nH, amax, amin, outputUnits=None):

    nSigma = -2 * Ai * nH * np.pi * (amax**-0.5 - amin**-0.5)
    if outputUnits is not None:
        nSigma = nSigma.to(outputUnits)

    return nSigma

def getnd(Ai, nH, amax, amin, outputUnits=None):

    nd = -2/5 * Ai * nH * (amax**-2.5 - amin**-2.5)

    if outputUnits is not None:
        nd = nd.to(outputUnits)

    return nd

def getVH(Tgas, kb=constants["k"], mH=constants["m_p"], outputUnits=None):
    VH = np.sqrt((8*kb*Tgas)/(np.pi * mH))

    if outputUnits is not None:
        VH = VH.to(outputUnits)

    return VH

def getRH2(nH, VH, Sigma, nd, SH, epsilon, outputUnits=None):

    RH2 = 0.5 * nH * VH * Sigma * nd * SH * epsilon

    if outputUnits is not None:
        RH2 = RH2.to(outputUnits)

    return RH2



Ai = 7E-26 * u.cm**2.5
nH = 100*u.cm**-3
Tgas = 100*u.K
SH = 1
epsilon = 1
VH = getVH(Tgas, outputUnits=u.m/u.s)
print("Hydrogen velocity: %s\n" %VH)

amin_all = [50*u.angstrom, 100*u.angstrom]
amax_all = [2500*u.angstrom, 2500*u.angstrom]

Sigmas = []
sigmas = []
nds = []

for i in range(len(amin_all)):
    amin = amin_all[i]
    amax = amax_all[i]
    Sigma = getSigma(Ai, nH, amax, amin, outputUnits=u.cm**-1)
    Sigmas.append(Sigma)
    print("Integrated area Sigma for %s to %s  : %s" %(amin, amax, Sigma))

    nd = getnd(Ai, nH, amax, amin, outputUnits=u.cm**-3)
    nds.append(nd)
    print("Dust density : %s" %nd)

    sigma = Sigma/nd
    sigmas.append(sigma)
    print("little sigma (for RH2 equation) : %s" %sigma)

    RH2 = getRH2(nH, VH, sigma, nd, SH, epsilon, outputUnits=u.cm**-3/u.s)
    print(RH2 / nH ** 2)
    print("H2 formation rate : %s\n" %RH2)

print("Irradiated formation rate = %s" %(10E-16 * nH**2 * u.cm**3/u.s))












