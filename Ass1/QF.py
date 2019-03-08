import numpy as np
import math
import json
import os
AE4878_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(AE4878_path, 'constants_matt.json')) as  handle:
    course_constants = json.loads(handle.read())

def get_D22(M2, W22sig, psi22):
    bottom = W22sig*np.sin(psi22)
    return M2/bottom

def get_W22sig(R, rho, k2, h2, sigma, A22):
    pi = np.pi

    return 4 * pi * R**2 * rho * (1 + k2 - h2) * sigma * A22

def get_A22(g, H, f22):
    return g * H * f22

def get_fnm(n,m):
    pi=np.pi
    Nnm = (2/(2*n+1))*(math.factorial(n+m))/(math.factorial(n-m))
    return (2*pi*Nnm)**(-0.5)*(-1)**m

def get_sigma():
    pi = np.pi
    omega_E = 2*pi/(24*60*60)
    omega_M = 2*pi/(30.4167*24*60*60)
    return 2*(omega_E -omega_M)

g = 9.81 #m/s^2
H = 0.63194
RE = 6378136.0 #m
rho = 1025 #kg/m^3
k2 = 0.303
h2 = 0.612
sigma = get_sigma() #rad/s
M2 = 2.416E12
psi22 = np.deg2rad(60)

f22 = get_fnm(2,2)
A22 = get_A22(g, H, f22)
W22sig = get_W22sig(RE, rho, k2, h2, sigma, A22)

D22 = get_D22(M2, W22sig, psi22)
print('D22 = %s cm' %(D22*100))