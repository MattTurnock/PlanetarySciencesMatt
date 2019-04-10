import numpy as np
from astropy import units as u

Pmax_PL = 44.71*u.W
PL_power_frac = 29.3/100

P_ccd = 11*u.W
Pmax_PL = Pmax_PL - 0.5*P_ccd

Pmax = Pmax_PL/PL_power_frac
P_ttc = (24/100) * Pmax
Pmax = Pmax - 0.5*P_ttc

print("Maximum power usage (sum of all) : %s" %Pmax)
t_total = 10*u.day
print("Operating time : %s" %t_total)
Emax = (Pmax * t_total).to(u.W*u.hour)
print("Maximum energy usage Emax : %s" %Emax)

#####################################################
# Primary batteries only
Espec = 800 *u.Unit("W*hour/kg")
print("\n====PRIMARY BATTERIES==== \nSpecific energy : %s" %Espec)
m_batt = Emax/Espec
print("Total battery (and power system) mass : %s" %m_batt)

#####################################################
# Fuel Cells only
Pspec_FC = 275 * u.W/u.kg
m_FC = Pmax / Pspec_FC
print("\n====FUEL CELLS==== \nWith specific power %s, mass of fuel cell is : %s" %(Pspec_FC, m_FC))

mdot_fuel = 0.36 *u.Unit("kg * kW**-1 * hour**-1")
m_fuel = (mdot_fuel * Emax).to(u.kg)
print("With fuel flow %s, fuel mass is : %s" %(mdot_fuel, m_fuel))
m_FC_total = m_FC + m_fuel
print("Total system mass is : %s" %m_FC_total)






