from astropy import units as u
from json_to_dict import constants
import numpy as np
pi = np.pi

######################################################################################

def get_K(Ms, Mp, P, i, e=0*u.one, G=constants["G"]):
    top = (2*pi*G)**(1/3) * Mp * np.sin(i)
    bottom = P**(1/3) * (Ms + Mp)**(2/3) * np.sqrt(1-e**2)
    K = top/bottom
    K = K.decompose()
    return K

def get_Mp(Ms, K, P, e=0*u.one, G=constants["G"], i=0*u.deg, units=True):
    if units:
        Ms = (Ms.to(u.kg)).value
        K = (K.to(u.m/u.s)).value
        P = (P.to(u.s)).value
        e = e.value
        G = (G.to(u.Unit("m**3/(kg*s*s)"))).value
        i = (i.to(u.rad)).value

    A = ( K * P**(1/3)/(2*pi*G)**(1/3))**1.5 * (np.sqrt(1-e**2))**1.5
    Mp = A*Ms/(np.sin(i)-A)

    if units:
        Mp = Mp*u.kg

    return Mp

def get_Mp_2(Ms, K_true, P, minMp=0, maxMp=3, iters=100, i=90*u.deg):

    Mps = np.linspace(minMp,maxMp,iters) * constants["MJupiter"]
    Kdiffs = []
    for Mp in Mps:
        K = get_K(Ms, Mp, P, i)
        Kdiffs.append((K-K_true).value)
    Kdiffs = np.abs(np.array(Kdiffs))
    mindex = np.argmin((Kdiffs))

    return Kdiffs[mindex], Mps[mindex]/ constants["MJupiter"]

def get_a(mu, P, units=True, unit_return=u.AU):
    inner = mu*(P/(2*pi))**2
    a = inner**(1/3)
    if units:
        a = a.decompose()
        a = a.to(unit_return)

    return a

Ms = constants["MSun"]
K1 = 32*u.m/u.s
P1 = 150 * constants["TEarth"]
K2 = 88*u.m/u.s
P2 = 1200 * constants["TEarth"]
i=90*u.deg

#######################################################################################################################
# Part D
print("=====================================Part D============================================")
calc1=True
if calc1:
    Kdiff1, Mp1 = get_Mp_2(Ms, K1, P1, minMp=0, maxMp=1, iters=10000, i=i)
    print("(Minimum) mass of planet 1 (in Jupiter masses), Mp1: %s \nWith K error: %s m/s \n" %(Mp1, Kdiff1))

    Kdiff2, Mp2 = get_Mp_2(Ms, K2, P2, minMp=4, maxMp=5, iters=10000, i=i)
    print("(Minimum) mass of planet 2 (in Jupiter masses), Mp2: %s \nWith K error: %s m/s \n" %(Mp2, Kdiff2))


########################################################################################################################
# Part G
print("=====================================Part G============================================")
iG = 30*u.deg
calc2=True
if calc2:
    Kdiff1G, Mp1G = get_Mp_2(Ms, K1, P1, minMp=1, maxMp=2, iters=10000, i=iG)
    print("Mass of planet 1 (in Jupiter masses), for i=30 deg  Mp1G: %s \nWith K error: %s m/s \n" %(Mp1G, Kdiff1G))

    Kdiff2G, Mp2G = get_Mp_2(Ms, K2, P2, minMp=0, maxMp=100, iters=10000, i=iG)
    print("Mass of planet 2 (in Jupiter masses), for i=30 deg  Mp2G: %s \nWith K error: %s m/s \n" %(Mp2G, Kdiff2G))

########################################################################################################################
# Part H
print("=====================================Part H============================================")
mu = constants["G"]*constants["MSun"]
a1 = get_a(mu, P1)
print("Semi-major axis of planet 1, a1 (ie D1): %s \n" %a1)

a2 = get_a(mu, P2)
print("Semi-major axis of planet 2, a2 (ie D2): %s \n" %a2)

########################################################################################################################
# Part I
print("=====================================Part I============================================")

Mi = 0.8*Ms
mu_i = Mi*constants["G"]

calc3=True
if calc3:
    Kdiff1_I, Mp1_I = get_Mp_2(Mi, K1, P1, minMp=0, maxMp=1, iters=10000, i=i)
    print("(Minimum) mass of planet 1 (in Jupiter masses), for 0.8 solar mass star, Mp1_I: %s \nWith K error: %s m/s \n" %(Mp1_I, Kdiff1_I))

    Kdiff2_I, Mp2_I = get_Mp_2(Mi, K2, P2, minMp=3, maxMp=5, iters=10000, i=i)
    print("(Minimum) mass of planet 2 (in Jupiter masses), for 0.8 solar mass star, Mp2_I: %s \nWith K error: %s m/s \n" %(Mp2_I, Kdiff2_I))


a1_I = get_a(mu_i, P1)
print("Semi-major axis of planet 1, for a 0.8 solar mass star, a1_I (ie D1_I): %s \n" %a1_I)

a2_I = get_a(mu_i, P2)
print("Semi-major axis of planet 2, for a 0.8 solar mass star, a2_I (ie D2_I): %s \n" %a2_I)