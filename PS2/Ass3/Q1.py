from PS2.Ass1.ass1_utils import *
from PS2.Ass3.ass3_utils import *

from astropy import units as u



###################################################################################################################
AU = constants["AU"]

# Q1a - propagation from sun to earth of a photon
v_photon = constants["c"]
t_SE_photon = getPropTime(v_photon, d=AU, outputUnits=u.min)
print("Propagation time of X-ray photon: %s\n" %t_SE_photon)

# Q1b - propagation from sun to earth of a 10Mev proton
E_proton = 10*u.MeV
m_proton = constants["m_p"]
v_proton = getV(E_proton, m_proton, outputUnits=u.m/u.s)
v_proton_rel = v_proton/getLorentzFactor(v_proton/v_photon)

t_SE_proton = getPropTime(v_proton, d=AU, outputUnits=u.min)
t_SE_proton_rel = getPropTime(v_proton_rel, d=AU, outputUnits=u.min)
print("Proton velocity: %s (= %s c)" %(v_proton, v_proton/v_photon))
print("Proton relativistic velocity: %s (= %s c)" %(v_proton_rel, v_proton_rel/v_photon))
print("Propagation time of proton: %s" %t_SE_proton)
print("Propagation time of proton with relativity: %s\n" %t_SE_proton_rel)

# Q1c - propagation from sun to earth of a 10Mev proton
v_CMEs = [200*u.km/u.s, 450*u.km/u.s, 1000*u.km/u.s]
t_SE_CMEs = []
for i in range(len(v_CMEs)):
    v_CME = v_CMEs[i]
    t_CME = getPropTime(v_CME, d=AU, outputUnits=u.day)
    t_SE_CMEs.append(t_CME)
    print("Propagation time of CME with speed %s: %s" %(v_CME, t_CME))

