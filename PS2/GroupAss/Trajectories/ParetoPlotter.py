import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pareto


trajDataLoc = "trajectory.csv"


trajStringData = pd.read_csv(trajDataLoc, sep=',', usecols=[0, 1, 2, 9], dtype=str)
print(trajStringData.dtypes)
trajNumericData = pd.read_csv(trajDataLoc, sep=",", usecols=[3, 4, 5, 6, 7, 8], dtype=np.float)
print(trajNumericData.dtypes)

DVs = trajNumericData["Total DV (km/s)"]
transferTimes = trajNumericData["Duration (days)"]


print(trajNumericData.columns)
paretoPoints = np.array(pareto.eps_sort(trajNumericData, objectives=[0, 5]))

paretoDVs = paretoPoints[:, 5]
paretoDurations = paretoPoints[:, 0]
print(paretoDurations)

plt.plot(paretoDurations, paretoDVs, color="red")
plt.scatter(transferTimes, DVs)
plt.xlabel("Duration (days)")
plt.ylabel("Total DV (km/s)")
plt.grid(which="both")
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.savefig("paretoFront.pdf")
plt.show()