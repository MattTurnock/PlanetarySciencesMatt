import numpy as np
from os.path import join
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.cm import coolwarm
import matplotlib

from json_to_dict import constants
from PS2.Ass2.utils import *

titling = False

def getBouguer(rho, h, G=constants["G"], outputUnits=u.Unit("m/s^2"), returnValue=False):
    dgb = 2*np.pi*G*rho*h
    if outputUnits is not None:
        dgb = dgb.to(outputUnits)

    if returnValue is True:
        dgb = dgb.value

    return dgb

def getBouguerSimple(rho, h, G=constants["G"].value):
    dgb = 2 * np.pi * G * rho * h

    return dgb

MOLA_dataDir = "MOLA_data"

topoCsvLoc = join(MOLA_dataDir, "MolaPEDR_-90N90N_0E360E_20190311T134316682_topo_csv.csv")

colNames = ["LON", "LAT", "h"]
topoDataFull = pd.read_csv(topoCsvLoc, names=colNames, dtype=float, header=0, usecols=[0,1,2])
topoData = topoDataFull.iloc[::100, :]
# print(sorted(topoData.loc[:, "Bouguer"]))
topoData.loc[:, "Bouguer"] = getBouguerSimple((constants["rhoCrust_Mars"].to(u.Unit("kg/m^3"))).value, topoData.loc[:, "h"])


fig, ax = plt.subplots()
if titling: plt.title("Bouguer Plot")

s = ax.scatter(topoData.loc[:, "LON"], topoData.loc[:, "LAT"], c=topoData.loc[:, "Bouguer"], cmap=coolwarm)

plt.grid()
plt.ylim(-90,90)
plt.xlim(0,360)
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.savefig("BouguerPlot.png")
# plt.show()

fig, ax = plt.subplots()
if titling: plt.title("Topography Plot")
s = ax.scatter(topoData.loc[:, "LON"], topoData.loc[:, "LAT"], c=topoData.loc[:, "h"], cmap=coolwarm)


plt.grid()
plt.ylim(-90,90)
plt.xlim(0,360)
plt.xlabel("Longitude [deg]")
plt.ylabel("Latitude [deg]")
plt.savefig("TopographyPlot.png")
plt.show()
