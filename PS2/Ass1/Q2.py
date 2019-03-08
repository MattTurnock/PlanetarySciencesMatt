import numpy as np
from astropy import units as u
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

from json_to_dict import constants

def getV(E, m, outputUnits=u.m/u.s):
    V = np.sqrt(2*E/m)

    if outputUnits is not None:
        V = V.to(outputUnits)

    return V

def getr_lam(m, Vperp, q, B, outputUnits=u.m):
    r_lam = (m*Vperp)/(np.abs(q)*B)
    if outputUnits is not None:
        r_lam = r_lam.to(outputUnits)

    return r_lam

E = 80*u.keV
m = constants["m_e"]
q = constants["q_e"]
B =45000*u.nT

V = getV(E, m)
r_lam = getr_lam(m, V, q, B)

print("Larmor radius of a %s electron a magnetic field of %s (the average for the Earth): %s\n"
      "The above assumes that the whole velocity of the electron is perpendicular to the magnetic field" %(E, B, r_lam))
print("Gives a velocity of %s" %V)



#######################################################################################################################
# Part B

def getXYZ(t, Vperp, Vpar, r_lam):
    omega = Vperp/r_lam
    theta = omega*t

    z = r_lam * np.cos(theta)
    y = r_lam * np.sin(theta)

    x = Vpar * t

    return x,y,z

# V = V.value
alpha = np.deg2rad(25)

Vperp = V * np.sin(alpha)
Vpar = V * np.cos(alpha)
r_lam = getr_lam(m, Vperp, q, B)
r_lam = r_lam.value
Vperp = Vperp.value
Vpar = Vpar.value

times = np.linspace(0,0.00001, 1001)
outArray = np.zeros((len(times), 4))
outArray[:, 0] = times

for i in range(len(times)):
    t = times[i]
    outArray[i, 1:] = getXYZ(t, Vperp, Vpar, r_lam)

print(outArray)



fig = plt.figure()
ax = plt.axes(projection='3d')

xline, yline, zline = outArray[:, 1], outArray[:, 2], outArray[:, 3]
ax.plot3D(xline, yline, zline)
plt.xlabel("x [m]")
plt.ylabel("y [m]")

ax.set_zlabel("z [m]")

ax.view_init(elev=40, azim=210)
plt.savefig("plots/electronPlot.pdf")

# for ii in range(0, 360, 10):
#     ax.view_init(elev=40., azim=ii)
#     plt.savefig("plots/movie%d.png" % ii)

plt.show()
