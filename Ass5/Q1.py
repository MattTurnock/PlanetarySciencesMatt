from astropy import units as u
from json_to_dict import constants
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab
pi = np.pi

######################################################################################

######################################################################################
#PartB

#leftcol is daughter
leftcol = [0.12193,
           0.12260,
           0.12322,
           0.12371,
           0.12570,
           0.13283,
           0.13362,
           0.13656,
           0.13744,
           0.13286,
           0.13354,
           0.11953,
           0.14143]

#rightcol is parent
rightcol = [0.3360,
            0.3446,
            0.3521,
            0.3588,
            0.3835,
            0.4747,
            0.4839,
            0.5232,
            0.5330,
            0.4746,
            0.4828,
            0.3058,
            0.5841]


# calc the trendline
z = np.polyfit(rightcol, leftcol, 1)
p = np.poly1d(z)

trendXs = np.linspace(0, 0.8, 100)
Nd_t0 = p(0)
print("Initial daughter to stable nuclide ratio Nd_t0 (Os-187 / Os-188) occurs at y-intercept and is: %s" %Nd_t0)

show = False
save = False
do_Gas=True


plt.figure()
legend_plots = ['Linear Trendline']
legend_scatters = ['Plotted Points']
legend = legend_plots + legend_scatters
plt.scatter(rightcol, leftcol, marker='o', s=30)
pylab.plot(trendXs,p(trendXs), color='red', linewidth=1)
plt.grid()
# plt.xlim([0,0.15])
plt.xlabel("Np(t) ratio (Re-187 / Os-188)")
plt.ylabel("Nd(t) ratio (Os-187 / Os-188)")
if show: plt.legend(legend)
if save: plt.savefig("isochron_1.pdf")
if show: plt.show()


####################################################################################################################
# Part D

def doGaplot(t, t12, y_int, trendXs, t0=0):
    lmbda = np.log(2)/t12
    slope = np.exp(lmbda*(t-t0)) - 1
    def y(x, m, c):
        y = m*x + c
        return y

    trendXs = trendXs
    trendYs = y(trendXs, slope, y_int)
    # print(trendYs)

    plt.plot(trendXs, trendYs, linestyle='--')

    return None

t12 = 41.2e9


if do_Gas:
    Gas = [0,1,2,3,4]
    for Ga in Gas:
        t = Ga*10**9
        doGaplot(t, t12, Nd_t0, trendXs, t0=0)
        legend_plots.append("Theoretical isochron at t=%s Ga" %Ga)
    legend = legend_plots + legend_scatters
    plt.legend(legend)
    if save: plt.savefig("isochron_2.pdf")
    plt.show()


####################################################################################################################
# Part E THIS IS CURRENTLY WRONG K?

def get_tmt0(tau12, Np_t, Nd_t, Nd_t0):
    term1 = tau12/np.log(2)
    frac = Np_t/(Np_t + Nd_t - Nd_t0)
    tmt0 = -term1*np.log2(frac)

    return tmt0

tau12 = t12
Np_t = np.mean(rightcol)
Nd_t = np.mean(leftcol)
Nd_t0 = Nd_t0

print("Half life = %s Gyrs" %(tau12*10**-9))
print("Np(t) = %s \nNd(t) = %s \nNd(t0) = %s \n" %(Np_t, Nd_t, Nd_t0))

age = get_tmt0(tau12, Np_t, Nd_t, Nd_t0)
print("CHECK Meteorite age = %s Gyrs" %(age*10**-9))

