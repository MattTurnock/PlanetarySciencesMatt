from astropy import units as u
from json_to_dict import constants
from Q2 import get_g_r
import numpy as np
pi = np.pi

############################################
#Part a)

def get_omega(T):
    omega = 2*pi/T
    return omega

def get_g_eff(mu, R, T, ret_g_c=False):
    omega = get_omega(T)
    g_p = get_g_r(mu, R)
    g_c = omega**2*R
    g_eff = g_p - g_c
    if ret_g_c:
        return g_eff.to(u.m/(u.s**2)), g_c.to(u.m/(u.s**2))
    else:
        return g_eff.to(u.m/(u.s**2))

def do_and_print(mu, R, T, body):
    g_p = get_g_r(mu, R)
    g_eff, g_c = get_g_eff(mu, R, T, ret_g_c=True)
    R_g = g_c / g_p
    first = body[0]

    print("Regular gravitational acceleration at surface of %s: g_p_%s = %s" %(body, first, g_p))
    print("Effective gravitational acceleration at surface of %s: g_eff_%s = %s" %(body, first, g_eff))
    print("Ratio between effective and regular acceleration at surface of %s: R_g_%s = %s %% \n" %(body, first, R_g*100))


do_and_print(constants["muEarth"], constants["RE"], constants["TEarth"], "Earth")
do_and_print(constants["muVenus"], constants["RVenus"], constants["TVenus"], "Venus")
do_and_print(constants["muJupiter"], constants["RJupiter"], constants["TJupiter"], "Jupiter")
