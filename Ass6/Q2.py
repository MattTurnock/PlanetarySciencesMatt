import numpy as np
from astropy import units as u

from json_to_dict import constants

pi = np.pi
#################################################################################################################
# Part a
print("\n========QUESTION A==========")

def getTemp(L, Ab, d, epsilon_ir, sigma=constants["sigma_SB"], outputUnits=None):

    print((L*(1-Ab)).decompose())
    print((16*d**2*epsilon_ir*sigma).decompose())
    T = (  (L*(1-Ab))/(16*pi*d**2*epsilon_ir*sigma)   )**(1/4)

    if outputUnits is not None:
        T = T.to(outputUnits)

    return T

Ab = 0.1
epsilon_ir = 0.9
d = 15*u.AU
L = 4E26*u.W

T = getTemp(L, Ab, d, epsilon_ir, outputUnits=u.K)
print("Comet temperature before outgassing : %s" %T)



#################################################################################################################
# Part c
print("\n========QUESTION C==========")

NA = constants["N_A"]
M_H2O = (2*1.008 + 15.999)*u.g/u.mol
m_H2O = (M_H2O/NA).to(u.kg)
print("H2O molar mass %s, mass %s" %(M_H2O, m_H2O))
constii = 1.2E22 * m_H2O
print("mdot expression constant : %s" %constii)

#################################################################################################################
# Part d
print("\n========QUESTION D==========")

def getDR(M, NA, t, d, rho, outputUnits=None):
    # d = (d.to(u.AU)).value
    mdot_ish = (1.2E22*M/NA)/u.s
    DR = mdot_ish * (t/(4*d**2*rho))

    if outputUnits is not None:
        DR = DR.to(outputUnits)

    return DR

T = 75*u.year
t = 0.1 * T
d_av = 1.5*u.AU
rho = 600*u.kg*u.m**-3

DR = getDR(M_H2O, NA, t, d_av, rho)
print("Change in ice layer thickness : %s" %DR.decompose())





