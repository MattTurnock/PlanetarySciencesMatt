import numpy as np
from astropy import units as u
from json_to_dict import constants
from PS2.Ass1.ass1_utils import *



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
