from PS2.Ass4.ass4_utils import *
from json_to_dict import constants

mu = constants["muEarth"]
J2 = constants["J2"]
Re = constants["RE"]

a_1 = 1330*u.km + constants["RE"]
i_1 = 66*u.deg
e_1 = 0

n_1 = getn(mu, a_1)
p_1 = getp(a_1, e_1)
omegaDot_1 = getomegaDot(n_1, J2, Re, p_1, i_1, outputUnits=u.rad/u.s)

omegaEDot_1 = omegaDot_1/-12.7
print(omegaEDot_1)












