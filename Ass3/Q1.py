from astropy import units as u
from json_to_dict import constants
import numpy as np
pi = np.pi

##############################################################
#Part (b)

def get_L(R, T, sigma_SB=constants["sigma_SB"], astropy_units = "W"):
    L = 4 * pi * R**2 * sigma_SB * T**4
    return L.to(astropy_units)

L_PC = get_L(constants["R_PC"], constants["T_PC"], astropy_units="W")
print('Luminosity of Proxima Centauri: %s (= %s) \n' %(L_PC, L_PC.to(u.solLum)))

##############################################################
#Part (c)

def get_F(L, r, astropy_units="W/(m**2)"):
    F = L/(4 * pi * r**2)
    return F.to(astropy_units)

F_PC = get_F(L_PC, constants["AU"])
#F_E = get_F(1*u.solLum, constants["AU"])
print('Flux at a distance of 1 AU: %s \n' %F_PC)

##############################################################
#Part (d)

def get_r_forF(L, F, astropy_units = "km"):
    r = np.sqrt(L/(4 * pi * F))
    return r.to(astropy_units)

r = get_r_forF(L_PC, constants["F_Earth"])
print('Orbital radius for equivalent flux to Earth: %s (=%s) \n' %(r, r.to(u.AU)))

##############################################################
#Part (e)

def get_lambda_max(T, wie_lam=2.9E-3*(u.m*u.K), units="micron"):
    lambda_max = wie_lam/T
    return lambda_max.to(units)

lambda_max = get_lambda_max(constants["T_PC"])
print('Proxima centauri peak flux wavelength: %s' %lambda_max)
