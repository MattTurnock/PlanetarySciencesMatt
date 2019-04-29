import numpy as np

def getThings(P_pl, P_pl_frac=0.22, P_else_fracs=[0.01, 0.15, 0.10, 0.18, 0.11, 0.12, 0.11]):

    P_tot = P_pl/P_pl_frac
    print("Total power : %s" %P_tot)

    P_else = []
    printTemplate = "Power of %s subsystem : %s"
    for i in range(len(P_else_fracs)):
        P_else_frac = P_else_fracs[i]
        P_else_temp = P_tot * P_else_frac

        if i==0:
            else_string = "structure"
        elif i==1:
            else_string = "thermal"
        elif i==2:
            else_string = "power"
        elif i==3:
            else_string = "TTC"
        elif i==4:
            else_string = "on-board processing"
        elif i==5:
            else_string = "ADCS"
        elif i==6:
            else_string = "propulsion"
        else:
            else_string = "u dun fuked up"

        P_else.append(P_else_temp)
        print("Power of %s subsystem : %s" %(else_string, P_else_temp))


getThings(38.48)
print("\n")
getThings(71.42)