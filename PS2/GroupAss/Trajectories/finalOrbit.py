import numpy as np
from astropy import units as u

from json_to_dict import constants



def geta(mu, T, outputUnits=None):

    a = np.cbrt( (mu * T**2)/(4*np.pi**2) )#^(1/3)

    if outputUnits is not None:
        a = a.to(outputUnits)

    return a

def getV(mu, r, a, outputUnits=None):

    V = np.sqrt(mu* (2/r - 1/a))

    if outputUnits is not None:
        V = V.to(outputUnits)

    return V

def getDV(mu, Ap1, Pe1, Ap2, Pe2, manPoint="Ap", outputUnits=None):

    a1 = 0.5 * (Ap1 + Pe1)
    a2 = 0.5 * (Ap2 + Pe2)

    if manPoint=="Ap":
        V1 = getV(mu, Ap1, a1)
        V2 = getV(mu, Ap2, a2)

    elif manPoint=="Pe":
        V1 = getV(mu, Pe1, a1)
        V2 = getV(mu, Pe2, a2)

    DV = abs(V2 - V1)

    if outputUnits is not None:
        DV = DV.to(outputUnits)

    return DV


muV = constants["muVenus"]
Rv = constants["RVenus"]
T = 0.1*u.day

a = geta(muV, T, outputUnits=u.km)
print("Semi-major axis for %s orbit : %s" %(T, a))
print("Circular orbit height : %s" %(a-Rv))

Ap1 = (Rv + 500*u.km)
Ap2 = Ap1
Pe1 = (Rv + 50*u.km)
Pe2 = Ap1

DV = getDV(muV, Ap1, Pe1, Ap2, Pe2, manPoint="Ap", outputUnits=u.m/u.s)
print("Delta V for circularistion of 500km orbit : %s" %DV)

#############################################################################
# Calcs for insertion
print("\n==================INSERTION CALCS===================")

# From parabolic to 500x250km (DV1)
Ap1_1 =  999999999999999999*u.km #Rv + 500*u.km
Ap2_1 = (Rv + 500*u.km)  #Rough height arbitrary
Pe1_1 = (Rv + 250*u.km)  # Height of atmosphere as given by Wikipedia
Pe2_1 = (Rv + 250*u.km)

DV1 = getDV(muV, Ap1_1, Pe1_1, Ap2_1, Pe2_1, manPoint="Pe", outputUnits=u.km/u.s)
print("DV1 from parabolic to 500x250 orbit : %s" %DV1)


# From 500x250 to 500x50, for balloon re-entry insertion
Ap1_2 = (Rv + 500*u.km)
Ap2_2 = (Rv + 500*u.km)
Pe1_2 = (Rv + 250*u.km)
Pe2_2 = (Rv + 50*u.km)

DV2 = getDV(muV, Ap1_2, Pe1_2, Ap2_2, Pe2_2, manPoint="Ap", outputUnits=u.km/u.s)
print("DV2 from 500x250 to 500x50 orbit : %s" %DV2)


# From 500x50 to 500x500, for final nominal orbiter
Ap1_3 = (Rv + 500*u.km)
Ap2_3 = (Rv + 500*u.km)
Pe1_3 = (Rv + 50*u.km)
Pe2_3 = (Rv + 500*u.km)

DV3 = getDV(muV, Ap1_3, Pe1_3, Ap2_3, Pe2_3, manPoint="Ap", outputUnits=u.km/u.s)
print("DV3 from 500x250 to 500x500 orbit : %s" %DV3)






