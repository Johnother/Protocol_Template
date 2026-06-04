import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit, fsolve
from uncertainties import ufloat
import uncertainties.unumpy as unp

### define custom functions here ###



def gauss(x, A, x0, sigma, C):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + C

def gaushermite(x, A, x0, w, C):
    return A * (x - x0)**2 * np.exp(-2 * (x - x0)**2 / w**2) + C
### end of functions ###


tem00 = np.genfromtxt("TEM00.txt")
tem01 = np.genfromtxt("TEM01.txt")


#fit für tem00 mode

x = tem00[:,0]
I = tem00[:,1]
I_err = tem00[:,2]


# Startwerte
p0 = [np.max(I), x[np.argmax(I)], 2.0, np.min(I)]

# gewichteter Fit mit Fehlern
popt, pcov = curve_fit(gauss,x,I,p0=p0,sigma=I_err,absolute_sigma=True)

A, x0, sigma, C = popt
A_err, x0_err, sigma_err, C_err = np.sqrt(np.diag(pcov))

print(f"A     = {A:.3f} ± {A_err:.3f}")
print(f"x0    = {x0:.3f} ± {x0_err:.3f} mm")
print(f"sigma = {sigma:.3f} ± {sigma_err:.3f} mm")
print(f"C     = {C:.3f} ± {C_err:.3f}")

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(x,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

x_fit = np.linspace(np.min(x), np.max(x), 1000)

ax.plot(x_fit,gauss(x_fit, *popt),label="Gauß-Fit")

ax.set_xlabel(r"x [mm]")
ax.set_ylabel(r"I [\mu A]")

ax.set_xlim(-16,16)
ax.set_ylim(0,8)

ax.grid()
ax.legend()

fig.savefig("../content/images/tem00.pdf")


#Fit für TEM01 mode

x = tem01[:,0]
I = tem01[:,1]
I_err = tem01[:,2]


p0 = [np.max(I),0.5,5.0,np.min(I)]

popt, pcov = curve_fit(gaushermite,x,I,p0=p0,sigma=I_err,absolute_sigma=True)

I0, x0, w, C = popt
I0_err, x0_err, w_err, C_err = np.sqrt(np.diag(pcov))

print("TEM01-Fit:")
print(f"I0 = {I0:.3f} ± {I0_err:.3f}")
print(f"x0 = {x0:.3f} ± {x0_err:.3f} mm")
print(f"w  = {w:.3f} ± {w_err:.3f} mm")
print(f"C  = {C:.3f} ± {C_err:.3f}")


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(x,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

x_fit = np.linspace(np.min(x), np.max(x), 1000)

ax.plot(x_fit,gaushermite(x_fit, *popt),label="TEM01-Fit")

ax.set_xlabel(r"x [mm]")
ax.set_ylabel(r"I [\mu A]")

ax.set_xlim(-16,16)
ax.set_ylim(0,8)

ax.grid()
ax.legend()

fig.savefig("../content/images/tem01.pdf")


#stabilitätsbedingung konvex konvex

konkon = np.genfromtxt("konkon.txt")


L = konkon[:,0]
L_err = konkon[:,1]
I = konkon[:,2]
I_err = konkon[:,3]


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(L,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")


ax.set_xlabel(r"L [m]")
ax.set_ylabel(r"I [mW]")

ax.set_xlim(np.min(L)-0.1,np.max(L)+0.1)
ax.set_ylim(0,7)

ax.grid()
ax.legend()

fig.savefig("../content/images/konkon.pdf")

#stabilitätsbedingung konvex konvex

plankon = np.genfromtxt("plankon.txt")


L = plankon[:,0]
L_err = plankon[:,1]
I = plankon[:,2]
I_err = plankon[:,3]


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(L,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")


ax.set_xlabel(r"L [m]")
ax.set_ylabel(r"I [mW]")

ax.set_xlim(np.min(L)-0.05,np.max(L)+0.05)
ax.set_ylim(0,5)

ax.grid()
ax.legend()

fig.savefig("../content/images/plankon.pdf")