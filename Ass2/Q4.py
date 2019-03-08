from astropy import units as u
from json_to_dict import constants
from Q2 import get_g_r
import numpy as np
pi = np.pi

######################################################
#part a)
def get_Htot(R, H):
    Htot = 4*pi*R**2*H
    return Htot.to(u.W)

H = 4E-2*(u.J/(u.m**2*u.s))
Htot = get_Htot(constants["RMars"], H)

Q_year = Htot.to(u.J/u.year)
print("Total heat lost by Mars each year due to heat flux: Q_year = %s\n" %Q_year)

#####################################################
#part b)

def get_T(Htot, dt, c, M):
    T = Htot*dt/(c*M)
    return T.decompose()

cp = 1.2E3 *u.J/(u.kg*u.K)
M = constants["muMars"]/constants["G"]
tss = 4.571E9 *u.year

T_mars = get_T(Htot, tss, cp, M)
print("Temperature of Mars after the age of the solar system: T_mars = %s" %T_mars)

