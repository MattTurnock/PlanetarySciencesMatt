
from json_to_dict import constants

from PS2.Ass2.utils import *

hs = [constants["h_everest"],
      constants["h_olympusMons"],
      constants["h_hellasBasin"],
      constants["h_theiaMons"]]

rhoCrusts = [constants["rhoCrust_Earth"],
             constants["rhoCrust_Mars"],
             constants["rhoCrust_Mars"],
             constants["rhoCrust_Venus"]]
rhoMantles = [constants["rhoMantle_Earth"],
              constants["rhoMantle_Mars"],
              constants["rhoMantle_Mars"],
              constants["rhoMantle_Venus"]]

structures = ["Himalaya", "Olympus Mons", "Hellas Basin", "Theia Mons"]
rConts=[]

for i in range(len(hs)):
    h = hs[i]
    rhoCrust = rhoCrusts[i]
    rhoMantle = rhoMantles[i]
    structure = structures[i]

    rCont = getRCont(h, rhoCrust, rhoMantle, outputUnits="m")
    rConts.append(rCont)

    print("r_continent for %s: %s" %(structure, rCont))


