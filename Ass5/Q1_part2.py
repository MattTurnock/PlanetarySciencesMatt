from astropy import units as u
from json_to_dict import constants
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab
pi = np.pi

######################################################################################

######################################################################################
#Part F, G, H

def get_tmt0(N_t, N_t0, tau12, t0=0):
    log=np.log
    tmt0 = -tau12/log(2) * log(N_t/N_t0)

    return tmt0

tau12_Pd = 17.0 # days
N_t0_Pd = 1.00  # g
N_t_Pd = 0.125  # g
t_Pd = get_tmt0(N_t_Pd, N_t0_Pd, tau12_Pd)
t_Pd = t_Pd * u.d

tau12_Nb = 20000 # years
N_t0_Nb = 2.00   # g
N_t_Nb =  0.0625 # g
t_Nb = get_tmt0(N_t_Nb, N_t0_Nb, tau12_Nb)
t_Nb = t_Nb * u.year

print("For an isotope of Pd-103 with mass %s g, to decay to %s g, the time is: t = %s" %(N_t0_Pd, N_t_Pd, t_Pd))
print("For an isotope of Nb-103 with mass %s g, to decay to %s g, the time is: t = %s" %(N_t0_Nb, N_t_Nb, t_Nb))

def get_tau12(t, N_t, N_t0, t0=0):
    log=np.log
    tau12 = -(t-t0) * log(2) / log(N_t/N_t0)

    return tau12

t_is = 75.0  # years
N_t_is = 1.25 # g
N_t0_is = 10.0 # g
tau12 = get_tau12(t_is, N_t_is, N_t0_is)
tau12 = tau12 * u.year
print("\nFor an isotope with mass %s g, to decay to %s g in 75 years, the half-life is: tau12 = %s" %(N_t0_is, N_t_is, tau12))


######################################################################################
#Part I

def get_Vt(D, rho_R, CD, rho_air=1.225*u.Unit("kg/m^3"), g=9.81*u.Unit("m/s^2")):
    Vt = np.sqrt( (4*D*rho_R*g)/(3*CD*rho_air) )
    Vt = Vt.decompose()
    return Vt

D = 15*u.m
rho_R = 2000*u.Unit("kg/m^3")
CD = 1.0
Vt = get_Vt(D, rho_R, CD)
print(Vt)

######################################################################################
#Part J
# CHECK THESE NUMBERINOS

def get_q(V, rho_air=1.225*u.Unit("kg/m^3")):
    q = 0.5 * rho_air * V**2
    q = q.to(u.Unit("bar"))
    return q

q = get_q(Vt)
print("\nCHECK Dynamic pressure of air on rock (or vice versa) is: q = %s (= %s)" %(q, q.to(u.kPa)))
P = q + 101325*u.Pa
print("CHECK Therefore total pressure is P = %s (= %s)" %(P, P.to(u.kPa)))

######################################################################################
#Part k

