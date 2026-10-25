import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
from uncertainties import ufloat, nominal_value, std_dev
from scipy.optimize import curve_fit

# Global constants


e = const.e # elementary charge
e_0 = const.epsilon_0
N_A = 6.022e23 # Avogadro's constant
m_a = 6.645e-27 # alpha mass in kg
E_a = 5.638*1e6*e # alpha energy in Joule
v_a = np.sqrt(2*E_a/m_a)
print(f"Alpha particle initial velocity = {v_a:.2f}")

# Atomic numbers
Z_Au = 79
Z_Bi = 83
Z_Al = 13

# Atomic densities #/m^3
N_Au = 19.32/197*N_A*1e6 # gcm^-3/gmol^-1*#mol^-1*1e6cm^3
N_Bi = 9.78/209*N_A*1e6
N_Al = 2.7/27*N_A*1e6

# Foil thickness in m 
D_Au = 2e-6 
D_Bi = 1e-6
D_Al = 3e-6

def Rutherford(theta,Z,D,N,A,t):
    return A*(Z**2*e**4*D*N)/((4*np.pi*e_0)**2*m_a**2*v_a**4)/(np.sin((theta-t)/2)**4)

# Calibration measurement to find the center of the beam

data = np.genfromtxt('calibration.csv',delimiter=',',unpack=True)
calibration = {"theta": data[0], "N": data[1], "p":data[2]}

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(calibration["theta"],calibration["N"]/90,xerr=0.1,yerr=np.sqrt(calibration["N"])/90,fmt=".",capsize=2,label=r"Calibration data")
ax.set_xlabel(r"$\theta$ [°]")
ax.set_ylabel(r"$N$ [#/s]")
ax.vlines(1,-5,30,ls="--",color="red",label="Maximum")
ax.set_xticks(np.arange(-5,6))
ax.set_ylim(-1,29)
ax.grid()
ax.legend()
fig.savefig("calibration.pdf")

# Rutherford-Experiment for Au-foil 2um

data = np.genfromtxt('Au2.csv',delimiter=",",unpack=True)
Au2 = {"theta":np.deg2rad(data[0]), "N":data[1], "t": data[2], "p": data[3]}

N = []
for i in range(len(Au2["N"])):
    num = Au2["N"][i]
    time = Au2["t"][i]
    N.append(ufloat(num,np.sqrt(num)) / ufloat(time,1))

fig, ax = plt.subplots(layout="constrained")

ax.errorbar(np.rad2deg(Au2["theta"][2:]), [i.n for i in N[2:]], xerr=0.1, yerr=[i.s for i in N[2:]], fmt="x", capsize=2, label=r"fitted data")
ax.errorbar(np.rad2deg(Au2["theta"][:2]), [i.n for i in N[:2]], xerr=0.1, yerr=[i.s for i in N[:2]], fmt="x", color="grey", capsize=2, label=r"Not fitted data")
popt, pcov = curve_fit(Rutherford, Au2["theta"][2:], [i.n for i in N[2:]], sigma=[i.s for i in N[2:]], absolute_sigma=True, p0=[Z_Au, D_Au, N_Au, 1, 0])
print(popt)
x = np.linspace(1e-6,np.pi/2,1000)
ax.plot(np.rad2deg(x),Rutherford(x,*popt), label=r"Fit ~ $sin(\theta/2)^{-4}$")
ax.set_ylim(-1,15)
ax.set_xlim(-1,20)
ax.set_xticks(np.arange(0,19,2))
ax.set_xlabel(r"$\theta$ [°]")
ax.set_ylabel(r"$N$ [#/s]")
ax.grid()
ax.legend()

fig.savefig("Au2.pdf")


# Rutherford-Experiment for Bi-foil 1um

data = np.genfromtxt('Bi1.csv',delimiter=",",unpack=True)
Bi1 = {"theta":np.deg2rad(data[0]), "N":data[1], "t": data[2], "p": data[3]}

N = []
for i in range(len(Bi1["N"])):
    num = Bi1["N"][i]
    time = Bi1["t"][i]
    N.append(ufloat(num,np.sqrt(num)) / ufloat(time,1))

fig, ax = plt.subplots(layout="constrained")

ax.errorbar(np.rad2deg(Bi1["theta"][2:]), [i.n for i in N[2:]], xerr=0.1, yerr=[i.s for i in N[2:]], fmt="x", capsize=2, label=r"Fitted data")
ax.errorbar(np.rad2deg(Bi1["theta"][:2]), [i.n for i in N[:2]], xerr=0.1, yerr=[i.s for i in N[:2]], fmt="x", capsize=2, color="grey", label=r"Not fitted data")
popt, pcov = curve_fit(Rutherford, Bi1["theta"][2:], [i.n for i in N[2:]], sigma=[i.s for i in N[2:]], absolute_sigma=True, p0=[Z_Bi, D_Bi, N_Bi, 1, 0])
print(popt)
x = np.linspace(1e-6,np.pi/2,1000)
ax.plot(np.rad2deg(x),Rutherford(x,*popt), label=r"Fit ~ $sin(\theta/2)^{-4}$")
ax.set_ylim(-1,21)
ax.set_xlim(-1,20)
ax.set_xticks(np.arange(0,19,2))
ax.set_yticks(np.arange(0,21,2))
ax.set_xlabel(r"$\theta$ [°]")
ax.set_ylabel(r"$N$ [#/s]")
ax.grid()
ax.legend()
fig.savefig("Bi1.pdf")

# Rutherford-Experiment for Al-foil 3um

data = np.genfromtxt('Al3.csv',delimiter=",",unpack=True)
Al3 = {"theta":np.deg2rad(data[0]), "N":data[1], "t": data[2], "p": data[3]}

N = []
for i in range(len(Al3["N"])):
    num = Al3["N"][i]
    time = Al3["t"][i]
    N.append(ufloat(num,np.sqrt(num)) / ufloat(time,1))

fig, ax = plt.subplots(layout="constrained")

ax.errorbar(np.rad2deg(Al3["theta"]),[i.n for i in N],xerr=0.1,yerr=[i.s for i in N],fmt="x",capsize=2, label=r"Fitted data")
popt, pcov = curve_fit(Rutherford,Al3["theta"][:],[i.n for i in N[:]],sigma=[i.s for i in N[:]],absolute_sigma=True, p0=[Z_Al, D_Al, N_Al, 1, 0])
print(popt)
x = np.linspace(1e-6,np.pi/2,1000)
ax.plot(np.rad2deg(x),Rutherford(x,*popt), label=r"Fit ~ $sin(\theta/2)^{-4}$")
ax.set_ylim(-1,17)
ax.set_xlim(-1,20)
ax.set_xticks(np.arange(0,19,2))
ax.set_xlabel(r"$\theta$ [°]")
ax.set_ylabel(r"$N$ [#/s]")
ax.legend()
ax.grid()
fig.savefig("Al3.pdf")