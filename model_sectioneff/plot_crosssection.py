import math
import matplotlib.pyplot as plt
import numpy as np

Mn = 0.985
Mv = 0.14


def cross_section_holo_2(t , Epho = 1.84, D_0 = -1.52, m_A = 1.612, m_D = 1.206, expo=3, A_0=0.430, N=2.032): #Formula in PHYSICAL REVIEW D 106, 086004 (2022)
    
    N_e = N #in nb GeV-2
    #N_e = 7.768 * 0.389 * 2/3
       
    s = (Mn**2 + 2*Mn*Epho)
    
    tilde_F = ((s-Mn**2)/2)**4
    
    eta = Mv**2 / ( 2*(s-Mn**2) - Mv**2 + t )

    #A_0 = 0.414
    A_t = A_0 / ( 1- (t/(m_A**2)) )**expo
    D_t = D_0 / ( 1- (t/(m_D**2)) )**expo

    
    sigma = N_e**2 * ( 1. / (64*math.pi*(s-Mn**2)**2) ) * ((A_t+(eta**2)*D_t)**2 / (A_0**2)) * tilde_F * 8
    return sigma

def Phi_dsigm_dt(tM, Eg = 3):
    Mp = 0.9383  # GeV
    M_Phi = 1.019  # GeV
    Pi = 3.14159
    Fermi = 1.e-13  # cm

    r = 1. * Fermi
    N_2g = 0.01


    t_slope = 1.13  

    F_2g = np.exp(t_slope * tM)
    s = Mp * Mp + 2 * Mp * Eg
    x = (2 * Mp * M_Phi + M_Phi * M_Phi) / (s - Mp * Mp)

    nue = 1 / (16 * Pi * (s - Mp * Mp) * (s - Mp * Mp))

    dSigm_dt = ((N_2g * nue * ((1 - x) ** 2)) / ((r ** 2) * (M_Phi) ** 2)) * F_2g * ((s - Mp * Mp) ** 2)
    print(nue)
    return dSigm_dt


t = [-i/100 for i in range(2000)]
Cs = []
for i in range(2000):
    Cs.append(cross_section_holo_2(t[i]))


plt.yscale("log")
plt.xlim(-20, 0)
plt.ylim(0.0000001, 2)
plt.xlabel("t")
plt.ylabel("dsigma/dt")
plt.plot(t, Cs)
plt.show()
