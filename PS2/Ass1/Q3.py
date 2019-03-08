import numpy as np
from astropy import units as u
from json_to_dict import constants

def getRhoS(rho_p, aR, rp, outputUnits=u.Unit("kg/m**3")):
    frac = rho_p**3/aR**3
    rho_s = (16*frac*rp**3)**(1/3)
    # print(rho_s.decompose())
    if outputUnits is not None:
        rho_s = rho_s.to(outputUnits)

    return rho_s

M_Neptune = constants["M_Neptune"]
R_Neptune = constants["R_Neptune"]
rhoNeptune = (M_Neptune/(4/3 * np.pi * R_Neptune**3)).to(u.Unit("kg/m**3"))

rNaiad = constants["rNaiad"]
rThallassa = constants["rThalassa"]
rDespina = constants["rDespina"]
rGalatea = constants["rGalatea"]
rLarissa = constants["rLarissa"]
rProteus = constants["rProteus"]

satellites_r = [rNaiad, rThallassa, rDespina, rGalatea, rLarissa, rProteus]
satellites_names = ["Naiad", "Thallassa", "Despina", "Galatea", "Larissa", "Proteus"]
satellites_rho = []

for i in range(len(satellites_names)):
    name = satellites_names[i]
    r = satellites_r[i]
    rho = getRhoS(rhoNeptune, r, R_Neptune)
    print("For the moon %s, the required density is %s at a radius of %s" %(name, rho, r))

