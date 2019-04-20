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
e1 = (Ap1 - Pe1)/(Ap1 + Pe1)
e2 = (Ap2 - Pe2)/(Ap2 + Pe2)
a1 = 0.5*(Ap1 + Pe1)
a2 = 0.5*(Ap2 + Pe2)


V1 = getV(muV, Ap1, a1, outputUnits=u.m/u.s)
V2 = getV(muV, Ap2, a2, outputUnits=u.m/u.s)
DV = V2 - V1
print("Delta V for circularistion of 500km orbit : %s" %DV)






