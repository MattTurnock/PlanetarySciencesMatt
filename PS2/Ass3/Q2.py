from PS2.Ass1.ass1_utils import *
from PS2.Ass3.ass3_utils import *

from astropy import units as u

def getMI(Z, rho, A, useAstroUnits=True):
    """
    Gets material index
    :param Z: atomic no
    :param rho: density
    :param A: mass number
    :return:
    """
    if useAstroUnits is True:
        rho = rho.to(u.g/u.cm**3)
    MI = Z * rho**-1 * A**(-2/3)

    return MI

def getAZ(As, Zs):
    """
    gets effective A and Z given lists of the mass numbers (As) and atomic numbers (Zs). Note that mass and atomic
    numbers must be for each atom in the molecule
    :param Ws:
    :param As:
    :param Zs:
    :return:
    """
    Asum = sum(As)
    Aeff = 0
    ZeffTemp = 0
    for i in range(len(As)):
        Ai = As[i]
        Zi = Zs[i]

        Wi = Ai/Asum
        Aeff += Wi*Ai

        ZeffTemp += Wi * Zi/Ai

    Zeff = Aeff * ZeffTemp

    return Aeff, Zeff

# Define mass number and atomic number for some elements
A_H = 1.008
Z_H = 1

A_O = 15.999
Z_O = 8

A_C = 12.011
Z_C = 6

A_Al = 26.982
Z_Al = 13

A_Pb = 207.2
Z_Pb = 82

materials = ["Liquid Hydrogen",
             "Water",
             "Polyethylene",
             "Aluminium",
             "Lead"]

As_all = [[A_H, A_H],
          [A_H, A_H, A_O],
          [A_C, A_C, A_H, A_H, A_H, A_H],
          [A_Al],
          [A_Pb]]


Zs_all = [[Z_H, Z_H],
          [Z_H, Z_H, Z_O],
          [Z_C, Z_C, Z_H, Z_H, Z_H, Z_H],
          [Z_Al],
          [Z_Pb]]

rhos = [70.99*u.g/u.L,
        0.997*u.g/u.mL,
        0.933*u.g/u.cm**3,
        2.70*u.g/u.cm**3,
        11.34*u.g/u.cm**3]
MIs=[]
Aeffs=[]
Zeffs=[]

for i in range(len(materials)):
    material = materials[i]
    As = As_all[i]
    Zs = Zs_all[i]
    # print(As, Zs)
    Aeff, Zeff = getAZ(As, Zs)
    Aeffs.append(Aeff)
    Zeffs.append(Zeff)

    rho = rhos[i]
    MI = getMI(Zeff, rho, Aeff, useAstroUnits=True)
    MIs.append(MI)
    print("For %s, Aeff = %s, Zeff = %s, rho = %s, MI = %s" %(material, Aeff, Zeff, rho.to(u.g/u.cm**3), MI.value))


# print(getAZ([1,1,16], [1,1,8]))
# print(getMI(7.963, rhos[1], 14.333))






