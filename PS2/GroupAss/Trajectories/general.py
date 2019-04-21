from PS2.GroupAss.Trajectories.trajectory_utils import *
from PS2.GroupAss.Trajectories.ParetoPlotter import *
u.imperial.enable()

print(trajNumericData.dtypes)

optionIndices = []
for i in range(len(trajNumericData["Post-injection DV (km/s)"])):
    piDV = trajNumericData["Post-injection DV (km/s)"].loc[i]
    totDV = trajNumericData["Total DV with C3 (km/s)"].loc[i]
    if (piDV <= 0.52) and (totDV <= 7.65):
        optionIndices.append(i)

print(trajNumericData.loc[optionIndices, "Post-injection DV (km/s)":"Total DV with C3 (km/s)"])
print(trajNumericData.loc[optionIndices[1]])
print(trajStringData.loc[optionIndices[1]])

def doTsiolkovsky(Isp, g, mInit, mFinal, outputUnits=None):

    DV = Isp * g * np.log(mInit/mFinal)

    if outputUnits is not None:
        DV = DV.to(outputUnits)

    return DV

g = 9.81*u.Unit("m/s**2")
lbUnit = u.imperial.lb
# Isp = 270*u.s

mPays = [100*u.kg, 200*u.kg, 400*u.kg]
mKicks = [2578.8*lbUnit, 2530.8*lbUnit, 3072*lbUnit, 1924*lbUnit, 1485.7*lbUnit, 1196.7*lbUnit, 810.9*lbUnit, 796.2*lbUnit, 579*lbUnit, 662.3*lbUnit]
mProps = [2345.3*lbUnit, 2350.1*lbUnit, 2835*lbUnit, 1698*lbUnit, 1392*lbUnit, 1113.6*lbUnit, 744.8*lbUnit, 735.6*lbUnit, 511.4*lbUnit, 601.6*lbUnit]
Isps = [293.7*u.s, 289.8*u.s, 293.5*u.s, 286.97*u.s, 290.4*u.s, 292.3*u.s, 291.4*u.s, 287.9*u.s, 272.1*u.s, 286.5*u.s]
stageNames = ["STAR 37FMV", "STAR 37FM", "STAR 31", "Orion 38", "STAR 30E", "STAR 30BP", "STAR 27H", "STAR 27", "STAR 26C", "STAR 20"]
DVs = []
for j in range(len(mPays)):
    mPay = mPays[j]
    print("\n")
    for i in range(len(mKicks)):

        mKick = mKicks[i]
        mProp = mProps[i]
        Isp = Isps[i]

        mInit = mPay + mKick
        mFinal = mInit - mProp

        DV = doTsiolkovsky(Isp, g, mInit, mFinal, outputUnits=u.km/u.s)
        DVs.append(DV)
        print("%s has mass %s (%s), and provides %s for a %s bus" %(stageNames[i], mKick, mKick.to(u.kg),  DV, mPay))






