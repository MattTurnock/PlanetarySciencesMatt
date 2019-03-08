import numpy as np
from scipy import special as sp

import json
import os
AE4878_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(AE4878_path, 'constants_matt.json')) as  handle:
    course_constants = json.loads(handle.read())

from matplotlib import pyplot as plt

def get_Ua_n(mu_2, r_1, r_12, n, psi):
    term1 = mu_2/r_12
    term2 = (r_1/r_12)**n
    Pn = sp.legendre(n)
    Ua_n = term1*term2*Pn(np.cos(psi))
    return Ua_n

#define legendre polynomials
P2 = sp.legendre(2)
P3 = sp.legendre(3)

mu_M = 4.9048695E12
r_E = course_constants['RE']['val']*1000
r_EM = 60*r_E
g = 9.81
psi=2*3.14159/4#np.deg2rad(90)
#print(mu_M, r_E)

N = range(2,6)
for n in N:
    Ua_n = get_Ua_n(mu_M, r_E, r_EM, n, psi)
    #print((Ua_n/g)*100)

def P2(x):
    return 0.5*(3*x**2 - 1)

pi = np.pi
plot=False
if plot:
    xCenter = 0
    yCenter = 0
    xRadius = r_E + 36.3/100
    yRadius = r_E + 9.07/100
    theta = np.linspace(0, pi*2, 1000)
    x = xRadius * np.cos(theta) + xCenter
    y = yRadius * np.sin(theta) + yCenter
    plt.plot(x, y)
    plt.grid()
    plt.show()



#Define all final constants
bodies = ['Sun','Moon', 'Jupiter']
r_E = course_constants['RE']['val']
psi_primary=0
psi_secondary = np.deg2rad(90)
Ns = [2,3,4]
g = 9.80665/1000
AU = course_constants['AU']['val']
precision = 100
tides = {}
for body in bodies:
    primary_bulges = []
    secondary_bulges = []
    mu_2 = course_constants['mu'+body]['val']
    if body == 'Sun':
        r_12 = course_constants['AU']['val']
    elif body == 'Jupiter':
        r_12_max = course_constants['r'+body]['val']*AU + AU
        r_12_min = course_constants['r'+body]['val']*AU - AU
        r_12s = [r_12_max, r_12_min]
        #print(r_12s)
    else:
        r_12 = course_constants['r'+body]['val']

    if body == 'Jupiter':
        for r_12 in r_12s:
            primary_bulges = []
            secondary_bulges = []
            for n in Ns:
                primary_bulges.append(round(1000*100*get_Ua_n(mu_2, r_E, r_12, n, psi_primary)/g, precision))
                secondary_bulges.append(round(1000*100*get_Ua_n(mu_2, r_E, r_12, n, psi_secondary)/g, precision))
            if r_12 == r_12_max:
                tides[body + '_rmax_primary'] = primary_bulges
                tides[body + '_rmax_secondary'] = secondary_bulges
            elif r_12 == r_12_min:
                tides[body + '_rmin_primary'] = primary_bulges
                tides[body + '_rmin_secondary'] = secondary_bulges
    else:
        for n in Ns:
            primary_bulges.append(round(1000*100*get_Ua_n(mu_2, r_E, r_12, n, psi_primary)/g, precision))
            secondary_bulges.append(round(1000*100*get_Ua_n(mu_2, r_E, r_12, n, psi_secondary)/g, precision))
        tides[body + '_primary'] = primary_bulges
        tides[body + '_secondary'] = secondary_bulges

#Values in cm
for tide in tides:
    a=4
    print(tide, tides[tide])

