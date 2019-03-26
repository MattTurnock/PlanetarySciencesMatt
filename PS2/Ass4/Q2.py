from PS2.Ass4.ass4_utils import *
from json_to_dict import constants

mu = constants["muEarth"]
J2 = constants["J2"]
Re = constants["RE"]

a_1 = 1330*u.km + constants["RE"]
print(a_1)
i_1 = 66*u.deg
e_1 = 0
Nr_1 = 127
Nd_1 = -10

p_1 = getp(a_1, e_1)

thetaDot = ((2*np.pi*u.rad)/(u.sday)).to(u.rad/u.s)
print("Earth rotation rate thetaDot : %s" %thetaDot)

n_1 = getn(mu, a_1)
print("mean motion n (==f_dot) : %s" %n_1)
T_1 = ((2*np.pi*u.rad/n_1)*127).to(u.sday)
print(T_1)

omegaDot_1 = getomegaDot(n_1, J2, Re, p_1, i_1, outputUnits=u.rad/u.s)
print("omegaDot_1 : %s" %omegaDot_1)

OmegaDot_1 = getOmegaDot(n_1, J2, Re, p_1, i_1, outputUnits=u.rad/u.s)
print("OmegaDot_1 : %s" %OmegaDot_1)

omegaEDot_1 = omegaDot_1/-12.7
print("omegaDotE_1 : %s\n" %omegaEDot_1)

print("Nr/Nd : %s" %(Nr_1/Nd_1))
print("omegaDot0_1/omegaDotE_1 : %s" %(omegaDot_1/omegaEDot_1))

#################################################################################################################
# Part b
print("\n")
a_2 = Re + 780*u.km
e_2=0

OmegaDotReq = (2*np.pi*u.rad/(1*u.year)).to(u.rad/u.s)
print("OmegaDot for sync : %s" %OmegaDotReq)

n_2 = getn(mu, a_2)
print("n_2 : %s" %n_2)

p_2 = getp(a_2, e_2)
i_GT = 98.12*u.deg

OmegaDotFound = getOmegaDot(n_2, J2, Re, p_2, i_GT)
print("OmegaDot for ground track : %s" %OmegaDotFound)


#################################################################################################################
# Part c
print("\n")

freqUnit = 1/u.sday
Dstar = constants["DS"]
Nd_2 = 3

f_M2S2 = 2*freqUnit
print("M2 frequency : %s" %f_M2S2)

fs_1 = abs((1/(Nd_1 * Dstar)).to(freqUnit))
print("Sampling frequency 1 : %s" %fs_1)

fs_2 = abs((1/(Nd_2 * Dstar)).to(freqUnit))
print("Sampling frequency 2 : %s\n" %fs_2)

fa_1 = getfa(fs_1, f_M2S2, outputUnits=freqUnit)
print("Aliasing frequency 1 : %s" %fa_1)

fa_2 = getfa(fs_2, f_M2S2, outputUnits=freqUnit)
print("Aliasing frequency 2 : %s" %fa_2)









