from astropy import units as u
from json_to_dict import constants
import numpy as np
pi = np.pi

##############################################################
#Part (a)

def get_H(T, g, mu_a, k=constants["k"], m_amu=constants["m_amu"], units=u.km):
    H = (k*T)/(g * mu_a * m_amu)
    return H.to(units)

def get_mu_a(ns, Ms):
    product = np.multiply(ns,Ms)
    mu_a = sum(product)
    return mu_a

T0 = (288*u.K)
g0 = 9.81*(u.m/u.s**2)

ns = [0.78, 0.22]
Ms = [2*14, 2*16]
mu_a = get_mu_a(ns, Ms)

H0 = get_H(T0, g0, mu_a)
print('Scale height of Earth atmosphere at surface: %s \n' %H0)

##############################################################
#Part (b)

def get_g_z(z, g0=9.81*(u.m/u.s**2), R=constants["RE"], units=u.m/u.s**2):
    r0 = R
    r1 = R + z
    g1 = g0*(r0/r1)**2
    return g1.to(units)

T1 = (215*u.K)
g0 = 9.81*(u.m/u.s**2)
z = 10*u.km
g1 = get_g_z(z)

ns = [0.78, 0.22]
Ms = [2*14, 2*16]
mu_a = get_mu_a(ns, Ms)

H1 = get_H(T1, g1, mu_a)
print('Scale height of Earth atmosphere at 10km: %s \n' %H1)

##############################################################
#Part (c)

def get_Pz(z, P0=1*u.bar, H=7.6*u.km, units="Pa", do_PzP0=False):
    Pz_P0 = np.exp(-z/H)
    if do_PzP0:
        return Pz_P0
    elif do_PzP0==False:
        Pz = Pz_P0 * P0
        return Pz.to(units)

z=100*u.km
Pz = get_Pz(z)
print('Atmospheric pressure according to exponential atmosphere at 100km: %s \n' %Pz)

##############################################################
#Part (d)

T_tit = (95*u.K)
g_tit = 1.35*(u.m/u.s**2)

ns_tit = [1.0]
Ms_tit = [2*14]
mu_a_tit = get_mu_a(ns, Ms)#

H_tit = get_H(T_tit, g_tit, mu_a_tit)
print('Scale height of Titan atmosphere at surface: %s \n' %H_tit)

##############################################################
#Part (e)

def get_z(H, Pz_p0, units=u.km):
    z = -H * np.log(Pz_p0)
    return z.to(units)

running=True
T_tit_av = 100*u.K
R_tit = 2576*u.km
z=0*u.km
dz = 100.0*u.m
while running:
    g_z = get_g_z(z, g0=g_tit, R=R_tit)
    H_z = get_H(T_tit_av, g_z, mu_a_tit)
    Pz_P0 = get_Pz(z, H=H_z, do_PzP0=True)
    ratio = 0.0002/100
    if Pz_P0 < ratio:
        g_z_final = g_z
        H_z_final = H_z
        Pz_P0_final = Pz_P0
        running = False

    else:
        z += dz

print('For the analysis of the height of Titans atmosphere where pressure is 0.0002 percent surface:')
print('Final pressure ratio: %s' %Pz_P0_final)
print('Final gravitational acceleration: %s' %g_z_final)
print('Final scale height: %s' %H_z_final)
print('Final height: %s \n' %z)

z = get_z(H_tit, 2E-6)
print('Compare this to the height if H was constant: %s' %z)