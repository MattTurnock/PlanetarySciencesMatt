from json_to_dict import constants
import itertools
from PS2.Ass4.MGOD_ass2_scripts.repeating_orbit_utils import do_a_iterations
import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt

def list_perms(list1, list2):
    perms_temp = list(itertools.product(list1, list2))
    perms=[]
    for perm in perms_temp:
        perms.append(list(perm))
    return perms
pi = np.pi

k=3
js = np.array([43])
jks = js/k
i_all = np.linspace(0,180,181)

#Function definition to do all calculations
def docalc_plot(i_all, js, k, approach='1', prnt=True, plot=True, plotshow=True, calc=True, extra='', interval=30, linewidth=1.5):
    if calc:
        input_data = list_perms(js,i_all)
        if prnt: print('data in: \n',input_data, '\n')
        out_data = []
        for combo in input_data:

            jk = combo[0]/k
            i = combo[1]*u.deg
            a = do_a_iterations(jk, i, approach=approach)
            h = a - constants["RE"]
            if 200*u.km < h < 1200*u.km:
                combo.append(h.value)
                out_data.append(combo)
        out_data = np.array(out_data)
        np.save('out_data%s' %extra, out_data)

    out_data = np.load('out_data%s.npy' %extra)
    if prnt: print('data out: \n',out_data,'\n')

    if plot:
        split_data = np.split(out_data, np.where(np.diff(out_data[:,0]))[0]+1)
        plt.figure()
        for data in split_data:
            xs = data[:,2]
            ys = data[:,1]
            plt.plot(xs, ys, linewidth=linewidth)
        legend=[]
        for j in js:
            string = "(j,k) = (%s, 3)" %j
            legend.append(string)
        plt.legend(legend, markerscale=5., bbox_to_anchor=(1.05, 1))
        plt.grid()
        plt.ylabel('Inclination, i [deg]')
        plt.xlabel('Orbital altitude, h [km]')
        plt.savefig('repeaters%s.pdf' %extra, bbox_inches="tight")
        if plotshow: plt.show()

    tabulated = out_data[~(out_data[:,1]%interval!=0.0), :]
    if prnt: print('tabulated data: \n', tabulated)
    np.savetxt('out_data%s.txt' %extra, tabulated, delimiter=' & ', newline=' \\\\\n', fmt='%i & %i & %1.2f')

#Do for approach 1 and approach 2
docalc_plot(i_all, js, k, approach='2', prnt=True, plot=True, plotshow=False, calc=True, extra='_2')
docalc_plot(i_all, js, k, approach='1', prnt=True, plot=True, plotshow=False, calc=True, extra='_1')


#Make a combined plot
out_data_1 = np.load('out_data_1.npy')
out_data_2 = np.load('out_data_2.npy')
split_data_1 = np.split(out_data_1, np.where(np.diff(out_data_1[:, 0]))[0] + 1)
split_data_2 = np.split(out_data_2, np.where(np.diff(out_data_2[:, 0]))[0] + 1)

# NOTE: split data 2 is for the approach 2 (solid line)
print("SPLIT DATA")
print(split_data_2[0])



linewidth=1.0
plot=True
plotshow=True
if plot:
    plt.figure()
    for data in split_data_1:
        xs = data[:, 2]
        ys = data[:, 1]
        plt.plot(xs, ys, linestyle='dashed', linewidth=linewidth )
    for data in split_data_2:
        xs = data[:, 2]
        ys = data[:, 1]
        plt.plot(xs, ys, linewidth=linewidth)
    legend = []
    for j in js:
        string = "(j,k) = (%s, 3)" % j
        legend.append(string)
    plt.legend(legend, markerscale=5.,bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.grid()
    plt.ylabel('Inclination, i [deg]')
    plt.xlabel('Orbital altitude, h [km]')
    plt.savefig('repeaters_both.pdf',bbox_inches="tight")
    if plotshow: plt.show()

