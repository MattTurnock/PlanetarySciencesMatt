import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from json_to_dict import constants

def getBr(Mb, r, theta, outputUnits=u.nT):
    Br = - 2 * (Mb/r**3) * np.cos(theta)
    if outputUnits is not None:
        Br = Br.to(outputUnits)

    return Br

def getBtheta(Mb, r, theta, outputUnits=u.nT):
    Btheta = (Mb/r**3) * np.cos(theta)
    if outputUnits is not None:
        Btheta = Btheta.to(outputUnits)

    return Btheta

# def getBr(M, r, theta, mu0=constants["mu0"], outputUnits=u.nT):
#     Br = (mu0 * M * np.sin(theta))/(2*np.pi*r**3)
#     if outputUnits is not None:
#         # Br = Br.to(outputUnits)
#         Br = Br.decompose()
#     return Br

Mb = constants["Mb_Jupiter"]
r = constants["rEuropa"]
print(Mb)
print(r)
# print(getBr(Mb, r, (9.6)*u.deg))
print(getBr(Mb, r, 90*u.deg))

thetas = np.linspace(0*u.deg, 360*u.deg, 1000)
thetaVals = []
Brs = []
BrVals = []
Bthetas = []
BthetaVals = []
for i in range(len(thetas)):
    theta = thetas[i]
    thetaVals.append(theta.value)

    Br = getBr(Mb, r, theta)
    Brs.append(Br)
    BrVals.append(Br.value)

    Btheta = getBtheta(Mb, r, theta)
    Bthetas.append(Btheta)
    BthetaVals.append(Btheta.value)

thetaVals = np.array(thetaVals)
thetaCoVals = thetaVals

plt.figure()
plt.plot(thetaVals, BrVals)
plt.grid()
plt.xlabel(r"$\phi$ [deg]")
plt.ylabel(r"$B_r$ [nT]")
plt.savefig("plots/Q1_Br.pdf")

plt.figure()
plt.plot(thetaVals, BthetaVals)
plt.grid()
plt.xlabel(r"$\phi$ [deg]")
plt.ylabel(r"$B_{\theta}$ [nT]")
plt.savefig("plots/Q1_Btheta.pdf")


# plt.show()
