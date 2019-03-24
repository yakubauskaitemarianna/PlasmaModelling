from data_getters import Data
from methods import Runge_Kutta_IV_3var_1step
from math import pi
import scipy
from scipy import integrate
from scipy import arange
from prettytable import PrettyTable

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_To_m(I):
    table_I_T0 = Data.get_table_I_T0()
    len_I_T0 = len(table_I_T0)

    In, Ik = 0, 1

    if I > 1200:
        In, Ik = len_I_T0 - 2,  len_I_T0 - 1

    for i in range(len_I_T0 - 1):
        if table_I_T0[i][0] < I < table_I_T0[i + 1][0]:
            In = i
            Ik = i + 1

    X0 = table_I_T0[In][0]
    X1 = table_I_T0[Ik][0]
    Y0 = table_I_T0[In][1]
    Y1 = table_I_T0[Ik][1]
    To = Y0 + (Y1 - Y0) * (I - X0) / (X1 - X0)
    Y0 = table_I_T0[In][2]
    Y1 = table_I_T0[Ik][2]
    m = Y0 + (Y1 - Y0) * (I - X0) / (X1 - X0)
    return To, m


def get_T(I, r):
    To, m = get_To_m(I)
    if m < 0: 
        m = 1
    return (Data.T_w - To) * (r / Data.Rad) ** m + To


def get_sigma(T):
    sigma_table = Data.get_table_sigma()
    len_sigma = len(sigma_table)

    Tn, Tk = 0, 1
    if T > 14000:
        Tn, Tk = len_sigma - 2, len_sigma - 1

    for i in range(len_sigma - 1):
        if sigma_table[i][0] < T < sigma_table[i + 1][0]:
            Tn = i
            Tk = i + 1

    X0 = sigma_table[Tn][0]
    X1 = sigma_table[Tk][0]
    Y0 = sigma_table[Tn][1]
    Y1 = sigma_table[Tk][1]
    sigma = Y0 + (Y1 - Y0) * (T - X0) / (X1 - X0)
    return sigma


def sigma_integrand(r, I):
    sigma = get_sigma(get_T(I, r))
    # print("R:", r, "\n")
    return sigma * r


def get_Rp(r, I):
    x = []
    r0 = 0
    step = 1e3
    for i in range(0, 1000):
        x.append(r0)
        r0 += Data.Rad/step
    y = [sigma_integrand(a, I) for a in x]
    return Data.l_e / (2 * pi * scipy.integrate.simps(y, x))


def dI_dt(t, Uc, I):
    Rp = get_Rp(t, I)
    #Rp=0
    return (Uc - (Data.R_k + Rp) * I) / Data.L_k


def dUc_dt(t, Uc, I):
    return -1 / Data.C_k * I


if __name__ == '__main__':
    ts = list(arange(Data.t_0, (Data.t_n + Data.tau), Data.tau))
    I, Uc, Rp = [Data.I_0], [Data.U_0], [get_Rp(Data.t_0, Data.I_0)]

    Ii = Data.I_0
    Ui = Data.U_0
    Rpi = get_Rp(Data.t_0, Ii)
    Rpi = 0
    result = PrettyTable()
    result.field_names = ["T", "I", "Uc", "Rp"]
   


    for t in ts[1:]:
        Uc_next, I_next = Runge_Kutta_IV_3var_1step(dI_dt, dUc_dt, t, Ui, Ii, Data.tau)
        Rpi = get_Rp(t, Ii)
        #Rpi = 0

        I.append(I_next)
        Uc.append(Uc_next)
        Rp.append(Rpi)

        Ii = I_next
        Ui = Uc_next

        result.add_row(["{:4.2e}".format(t), "{:4.5f}".format(Ii), "{:4.5f}".format(Ui), "{:4.5f}".format(Rpi)])

    result.add_row(["T_max = "+"{:4.2e}".format(ts[-1]), "I_max = "+"{:4.5f}".format(max(I)), "U_max "+"{:4.5f}".format(max(Uc)), "Rp_max "+"{:4.5f}".format(max(Rp))])
    result.add_row(["T_min = "+"{:4.2e}".format(ts[0]), "I_min = "+"{:4.5f}".format(min(I)), "U_min = "+"{:4.5f}".format(min(Uc)), "Rp_min = "+"{:4.5f}".format(min(Rp))])
    print(result)
    
    # График
    df = pd.DataFrame({'t' : ts, 'I' : I, 'Uc' : Uc, 'Rp' : Rp})
    plt.style.use('seaborn-darkgrid')
    palette = plt.get_cmap('Set1')
    num = 0
    for column in df.drop('t', axis=1):
        num += 1
        plt.plot(df['t'], df[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)
    plt.legend(loc=2, ncol=1)
    plt.title("I(t), Uc(t), Rp(t)")
    plt.xlabel('time, ms')
    plt.show()
    

    


