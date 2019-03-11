import numpy as np


from json_to_dict import constants

from PS2.Ass2.utils import *

densUnit = u.Unit("g/cm^3")

rhoIces = [0.937*densUnit,
           0.7*densUnit]

rhoWater = 1.0*densUnit

hs = [-2*u.km,
      -500*u.m,
      500*u.m,
      2*u.km,
      5*u.km]

rConts=np.zeros( ( len(hs), len(rhoIces) ), dtype=object )

for i in range(len(hs)):
    for j in range(len(rhoIces)):
        h = hs[i]
        rhoIce = rhoIces[j]

        rCont = getRCont(h, rhoIce, rhoWater, outputUnits=u.km)
        rConts[i, j] = rCont

        print("r_continent for h = %s and rhoIce = %s: %s" %(h, rhoIce, rCont))

print(rConts)
# npArray2LatexTable(rConts, "test.txt", header="cock")
header ="""
\\begin{table}[H]
    \centering
    \caption{Caption}
    \label{tab:my_label}
    \\begin{tabular}{}
"""
footer = """
    \end{tabular}
\end{table}
"""
npArray2LatexTable(rConts, "Q1c.txt")
# \begin{table}[H]
#     \centering
#     \caption{Caption}
#     \label{tab:my_label}
#     \begin{tabular}{}
#          &  \\
#          &
#     \end{tabular}
# \end{table}