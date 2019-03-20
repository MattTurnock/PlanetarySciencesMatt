import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pareto
from astropy import units as u
from json_to_dict import constants

from PS2.GroupAss.Trajectories.trajectory_utils import *

# print(C3toDV(2*u.Unit("km**2/s**2"), constants["muEarth"], (6378 + 100)*u.km, outputUnits=u.km/u.s))

trajDataLoc = "trajectory.csv"

printing = False
plotting = False

trajStringData = pd.read_csv(trajDataLoc, sep=',', usecols=[0, 1, 2, 9], dtype=str)
# print(trajStringData.dtypes)
trajNumericData = pd.read_csv(trajDataLoc, sep=",", usecols=[3, 4, 5, 6, 7, 8], dtype=np.float)
if printing: print(trajNumericData.dtypes)


# print(trajNumericData)





# Calculate DV from C3 assuming 150km orbit alt
muE = constants["muEarth"]
rOrbit = constants["RE"] + 150*u.km
C3Units = u.Unit("km**2/s**2")
outputUnits = u.km/u.s
trajNumericData["Total DV with C3 (km/s)"] = np.zeros(len(trajNumericData["Total DV (km/s)"]))


for i in range(len(trajNumericData["Total DV (km/s)"])):
    DV = trajNumericData["Total DV (km/s)"].loc[i]
    C3 = trajNumericData["C3 (km2/s2)"].loc[i]
    C3DV = C3toDV(C3*C3Units, muE, rOrbit, outputUnits=u.km/u.s).value
    DV_all = DV + C3DV
    trajNumericData["Total DV with C3 (km/s)"].loc[i] = DV_all

DVs = trajNumericData["Total DV (km/s)"]
DVS_C3s = trajNumericData["Total DV with C3 (km/s)"]
transferTimes = trajNumericData["Duration (days)"]









# trajNumericData["Total DV with C3 (km/s)"] = trajNumericData["Total DV (km/s)"] + C3toDV( trajNumericData["C3 (km2/s2)"] * C3Units, muE, rOrbit, outputUnits=u.km/u.s).value
# print(trajNumericData)

# Pareto's duration against total dV
paretoIndex1 = 0
paretoIndex2 = 6

# # Pareto's Post injection dV against total dV
# paretoIndex1 = 4
# paretoIndex2 = 6

paretoPointsDataframe = pareto.eps_sort(trajNumericData, objectives=[paretoIndex1, paretoIndex2])
paretoPoints = np.array(paretoPointsDataframe)

paretoDVs = paretoPoints[:, paretoIndex2]
paretoDurations = paretoPoints[:, paretoIndex1]

indices=[]
for i in range(len(paretoPoints)):
    for j in range(len(trajNumericData)):
        k = np.array(trajNumericData.loc[j])
        if list(paretoPoints[i]) == list(trajNumericData.loc[j]):
            indices.append(j)

for i in range(len(indices)):
    strData = trajStringData.loc[indices[i]]
    numData = trajNumericData.loc[indices[i]]
    if printing: print("\nOption %s, with string data:\n%s" %(i, strData))
    if printing: print("And numeric data:\n%s\n" %numData)

plt.figure()
plt.plot(paretoPoints[:, paretoIndex1], paretoPoints[:, paretoIndex2], color="red")
plt.scatter(transferTimes, DVS_C3s)
if paretoIndex1 == 0: xlabel = "Duration (days)"
elif paretoIndex1 == 4: xlabel = "Post-injection dV (km/s)"
plt.xlabel(xlabel)
plt.ylabel("Total DV (km/s)")
plt.grid(which="both")
ax = plt.gca()
# ax.set_yscale('log')
ax.set_xscale('log')

for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_formatter(ScalarFormatter())
plt.savefig("paretoFront.pdf")
if plotting: plt.show()
