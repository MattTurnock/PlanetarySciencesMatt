from astropy import units as u
from json_to_dict import constants
import numpy as np
pi = np.pi

##############################################################
#Part (a)
def get_rho_av(mu, R, G):
    rho_av = (mu / G) / ((4 / 3) * pi * R ** 3)
    return rho_av

def get_g_r(mu, r):
    g_r = mu/r**2
    return g_r.to(u.m/(u.s**2))

#must input as astropy unit types
def get_Pc(mu, R, G):
    rho_av = get_rho_av(mu, R, G)
    g_R = get_g_r(mu,R)
    return (rho_av*g_R*R).to(u.GPa)

#print((constants['muEarth']/constants['RE']**2).to(u.m/(u.s**2)))
Pc_E = get_Pc(constants['muEarth'], constants['RE'], constants['G'])
Pc_Moon = get_Pc(constants['muMoon'], constants['RMoon'], constants['G'])
Pc_Saturn = get_Pc(constants['muSaturn'], constants['RSaturn'], constants['G'])

prnt=False
if prnt:
    print("Pressure at the centre of the Earth (const g and rho): Pc_E = %s" %Pc_E)
    print("Pressure at the centre of the Moon (const g and rho): Pc_moon =  %s" %Pc_Moon)
    print("Pressure at the centre of Saturn (const g and rho): %s \n" %Pc_Saturn)

##############################################################
#Part (b)

#Needs unit variables
def get_P_r(mu,R,G,r):
    rho_av = get_rho_av(mu,R,G)
    P_r = 2/3 * rho_av**2 * np.pi * G * (R**2 - r**2)
    return P_r.to(u.GPa)

# Note that the following estimates are about half of what they should be - changing gravity is accounted for
# but not changing density, and gravity is 0 at core, therefore average g is 0.5*g_R
Pc_E_b = get_P_r(constants['muEarth'], constants['RE'], constants['G'],0*u.km)
Pc_Moon_b = get_P_r(constants['muMoon'], constants['RMoon'], constants['G'],0*u.km)
Pc_Saturn_b = get_P_r(constants['muSaturn'], constants['RSaturn'], constants['G'],0*u.km)

if prnt:
    print("Pressure at the centre of the Earth (changing g const rho): Pc_E = %s" %Pc_E_b)
    print("Pressure at the centre of the Moon (changing g const rho): Pc_moon =  %s" %Pc_Moon_b)
    print("Pressure at the centre of Saturn (changing g const rho): %s \n" %Pc_Saturn_b)

    print(get_rho_av(constants['muSaturn'], constants['RSaturn'], constants['G']).to(u.kg/(u.m)**3))