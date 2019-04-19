import numpy as np
from astropy import units as u
from json_to_dict import constants


######################################################################################
####### Set adjustable parameters ##########

# Fractions of time payloads are on. CCD, Spectrum analyser, and everything else
Pccd_av_frac = 0.5 # Assumes only active at night
Psa_av_frac = 1.0 # Assumes constantly active
Pelse_av_frac = 0.5 # Arbitrary activity fraction


# Mission time
t_total = 9 * u.day

######################################################################################
# Set known values for instrument power
Pasi_max = 3.9*u.W
Pneph_max = 1.56*u.W
Ptls_max = 13*u.W
Prad_max = 6.5*u.W
Pmag_max = 4.25*u.W
Pccd_max = 11*u.W
Psa_max = 0.5*u.W

# Set known values for power distribution
Ppl_frac = 29.3/100
Pstruct_frac = 1.3/100
Ptherm_frac = 1.3/100
Ppow_frac = 13.3/100
Pttc_frac = 24/100
Pproc_frac = 14.7/100
Padcs_frac = 16/100

# Known values for orbit
Rv = constants["RVenus"]
h = 500*u.km

######################################################################################
# Calculate derived parameters / budgets

# Power of most instruments
Pelse_max = Pasi_max + Pneph_max + Ptls_max + Prad_max + Pmag_max

# Calculate average power of instruments (and total average power)
Pccd_av = Pccd_av_frac * Pccd_max
Psa_av = Psa_av_frac * Psa_max
Pelse_av = Pelse_av_frac * Pelse_max
Ppl_tot_av = Pccd_av + Psa_av + Pelse_av

# Calculate maximum possible payload power usage in a given instant
Ppl_max = Pccd_max + Psa_max + Pelse_max

# Calculate power of different subsystems (and the total), based on average power usage
Ptot_av = Ppl_tot_av / Ppl_frac
Pstruct = Pstruct_frac * Ptot_av
Ptherm = Ptherm_frac * Ptot_av
Ppow = Ppow_frac * Ptot_av
Pttc = Pttc_frac * Ptot_av
Pproc = Pproc_frac * Ptot_av
Padcs = Padcs_frac * Ptot_av
Psubsystems = Pstruct + Ptherm + Ppow + Pttc + Pproc + Padcs

# Find fraction on-time
print(type(Rv), type(h))
alpha = np.arccos(Rv / (Rv + h)) # Alpha calculated from orbit geometry
C = alpha / (np.pi *u.rad) # Refers to fraction of time ttc is on
print(C)
# C = 0.5 # Refers to fraction of time ttc is on

# Calculate maximum possible spacecraft power usage (assumes payload fully operating, but others only sized based on average PL)
Ptot_max = Psubsystems + Ppl_max


#########################################################################################################
# Power budget calcs and printing for different options

print("Maximum power usage of s/c : %s" %Ptot_max)
print("Average power usage of s/c : %s" %Ptot_av)
print("Operating time : %s" %t_total)

#C refers to fraction of time ttc operates for
Ptot_av_energy = Ppl_tot_av + Pstruct + Ptherm + Ppow + C*Pttc + Pproc + Padcs
Etot = (Ptot_av_energy * t_total).to(u.W*u.hour)
print("Total energy usage Etot (ttc on %s of the time) : %s" %(C, Etot))

#####################################################
# Primary batteries only
Espec = 800  *u.Unit("W*hour/kg")
print("\n====PRIMARY BATTERIES==== \nSpecific energy : %s" %Espec)
m_batt = Etot/Espec
print("Total battery mass (assuming constant Ptot_av use) : %s" %m_batt)

#####################################################
# Fuel Cells only
Pspec_FC = 275 * u.W/u.kg
m_FC = Ptot_max / Pspec_FC
print("\n====FUEL CELLS==== \nWith specific power %s, mass of fuel cell is : %s" %(Pspec_FC, m_FC))
print("^^ Calculated using Ptot_max")

mdot_fuel = 0.36 *u.Unit("kg * kW**-1 * hour**-1")
m_fuel = (mdot_fuel * Etot).to(u.kg)
print("With fuel flow %s, fuel mass is : %s" %(mdot_fuel, m_fuel))
m_FC_total = m_FC + m_fuel
print("Total system mass is : %s" %m_FC_total)






