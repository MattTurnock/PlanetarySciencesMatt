import numpy as np
from astropy import units as u

from PS2.Ass5.ass5_utils import *
############################################################################################
# Part a) - magnitude

g = 9.8*u.Unit("m/s**2")
Omega = 7.29E-5*u.Unit("rad/s")
rho_av = 1027.5*u.Unit("kg/m**3")
phi = 36*u.deg
L = 200*u.km
Deltah = -1.1*u.m

f = getf(Omega, phi, outputUnits=u.rad/u.s)
print("Coriolis parameter f : %s" %f)

U = getU(f, Deltah, L, g=g, outputUnits=u.m/u.s)
print("Geostrophic speed u (from altimetry) : %s" %U)

rhoA = 1027.1*u.Unit("kg/m**3")
rhoB = 1027.9*u.Unit("kg/m**3")

P_level = 2000 * 0.1*u.bar
hA = (P_level/(rhoA * g)).to(u.m)
print("hA : %s" %hA)
hB = (P_level/(rhoB * g)).to(u.m)
print("hB : %s" %hB)
print("Delta h : %s" %(hA-hB))

L2 = 250*u.km
slope_AB = ((hA - hB)/L2).to(u.mm/u.km)
print("Slope from A to B : %s" %slope_AB)

U2 = getU2(hB, f, L2, rhoA, rhoB, outputUnits=u.m/u.s)
print("Geostrophic speed u (from oceanography) : %s" %U2)






