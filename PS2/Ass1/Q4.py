import numpy as np
from astropy import units as u
from json_to_dict import constants

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

muSaturn = constants["muSaturn"]

anoms = [80000*u.km, 80000*u.km, 120000*u.km]
Ws = [1*u.km, 100*u.km, 1*u.km]
DVs = []
tSpreads = []

for i in range(len(anoms)):
    anom = anoms[i]
    W = Ws[i]
    DV = getDV(muSaturn, anom, W)
    DVs.append(DV)
    tSpread = get_tSpread(anom, DV, outputUnits=u.year)
    tSpreads.append(tSpread)
    print("Radius %s; Width %s; Gives DV %s; tSpread %s" %(anom, W, DV, tSpread))

#######################################################################################################################

print("\n")

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

Vv = 0.01*u.m**2/u.s
tDiffs = []

for i in range(len(anoms)):
    anom = anoms[i]
    W1 = Ws[i]
    W2 = W1 * 2

    A1 = getAring2(anom, W1)
    A2 = getAring2(anom, W2)
    tDiff = get_tDiff(A1, A2, Vv, outputUnits=u.year)
    tDiffs.append(tDiff)
    print("Radius %s; Nominal width %s; tDiff %s" %(anom, W1, tDiff))


# anom = 80000*u.km
# W1 = 1*u.km
# W2 = W1*2
#
#
# A1 = getAring2(anom, W1)
# A2 = getAring2(anom, W2)
# tDiff = get_tDiff(A1, A2, Vv, outputUnits=u.year)
# print(tDiff)
