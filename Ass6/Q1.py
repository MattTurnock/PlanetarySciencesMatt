import numpy as np
from astropy import units as u

from json_to_dict import constants

pi = np.pi
#################################################################################################################
# Part c
print("\n========QUESTION C==========")
def getT(a, mu, outputUnits=None):

    T = 2*pi * np.sqrt(a**3/mu)

    if outputUnits is not None:
        T = T.to(outputUnits)

    return T

def geta(Ap, Pe, outputUnits=None):

    a = 0.5*(Ap + Pe)

    if outputUnits is not None:
        a = a.to(outputUnits)

    return a


muSun = constants["muSun"]
Ap = 16*u.AU
Pe = 0.5*u.AU
a = geta(Ap, Pe, u.AU)
T = getT(a, muSun, outputUnits=u.year)
print("Comet 2019/ABC orbit period with a = %s : %s" %(a, T))


#################################################################################################################
# Part d
print("\n========QUESTION D==========")

def getTP(aP, a, e, i, outputUnits=None):

    TP = aP/a + 2*np.sqrt((1-e**2) * a/aP) * np.cos(i)

    if outputUnits is not None:
        TP = TP.to(outputUnits)

    return TP

aP = constants["rJupiter"]
a1 = 4.76*u.AU
e1 = 0.781
i1 = 0*u.deg
a2 = 4.37*u.AU
e2 = 0.8004
i2 = 0*u.deg

TP1 = getTP(aP, a1, e1, i1, outputUnits=u.dimensionless_unscaled)
TP2 = getTP(aP, a2, e2, i2, outputUnits=u.dimensionless_unscaled)
DTP = abs(TP1 - TP2)
percDiffTP = 100*DTP / TP1


print("TP of orbit 1 : %s" %TP1)
print("TP of orbit 2 : %s" %TP2)
print("Percentage difference between values : %s" %percDiffTP)




