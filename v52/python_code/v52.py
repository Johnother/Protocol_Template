import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit, fsolve
from uncertainties import ufloat
import uncertainties.unumpy as unp

### functions ###






### Frequenzabhängigkeit ###


imp = np.genfromtxt("impedanz.txt")

f = imp[:,0]*10**3
L = unp.uarray(imp[:,1]*10**(-6),imp[:,2]*10**(-6))
R = unp.uarray(imp[:,3],imp[:,4])
C = unp.uarray(imp[:,5]*10**(-12),imp[:,6]*10**(-12))

omega = 2*np.pi*f

Z = unp.sqrt(R**2 + (omega**2)*L**2)

phi = unp.arctan((omega*L)/(R))* 180 / np.pi


with open("imp_phi.txt", "w") as g:
    for i in range(len(f)):
        g.write(f"{f[i]} & {Z[i].nominal_value:.5f} & {Z[i].std_dev:.5f} & {phi[i].nominal_value:.5f} & {phi[i].std_dev:.5f} \\\\\n")



fig, ax = plt.subplots(layout="constrained")
ax.errorbar(omega,unp.nominal_values(L),yerr=unp.std_devs(L),fmt="x",color="darkorange",capsize=3,label=r"$L\left(\omega\right)$")

ax.set_xlabel(r"$\omega$ [rad/s]")
ax.set_ylabel(r"$L$ [H]")
ax.set_title(r"Induktivitätsbelag $L\left(\omega\right)$")
ax.grid()
ax.legend()

fig.savefig("../content/images/ind.pdf")





fig, ax = plt.subplots(layout="constrained")

ax.errorbar(omega,unp.nominal_values(R),yerr=unp.std_devs(R),fmt="x",color="yellowgreen",capsize=3,label=r"$R\left(\omega\right)$")
ax.set_title(r"Ohmscher Belag $R\left(\omega\right)$")
ax.set_xlabel(r"$\omega$ [rad/s]")
ax.set_ylabel(r"$R$ [$\Omega$]")

ax.grid()
ax.legend()

fig.savefig("../content/images/ohm.pdf")




fig, ax = plt.subplots(layout="constrained")

ax.errorbar(omega,unp.nominal_values(Z),yerr=unp.std_devs(Z),fmt="x",color="steelblue",capsize=3,label=r"$\left|Z\right| \left(\omega\right)$")
ax.set_title(r"Impedanzbetrag $\left|Z\right| \left(\omega\right)$")
ax.set_xlabel(r"$\omega$ [rad/s]")
ax.set_ylabel(r"$Z$ [$\Omega$]")

ax.grid()
ax.legend()

fig.savefig("../content/images/imp.pdf")




fig, ax = plt.subplots(layout="constrained")

ax.errorbar(omega,unp.nominal_values(C),yerr=unp.std_devs(C),fmt="x",color="crimson",capsize=3,label=r"$C \left(\omega\right)$")
ax.set_title(r"Kapazitätsbelag $C \left(\omega\right)$")
ax.set_xlabel(r"$\omega$ [rad/s]")
ax.set_ylabel(r"$C$ [F]")

ax.grid()
ax.legend()

fig.savefig("../content/images/kap.pdf")





#### Bestimmung der Laenge ###

ld = np.genfromtxt("len_daem.txt")



index = ld[:,0]
delt = unp.uarray(ld[:,1]*10**(-9),ld[:,2]*10**(-9))
A0 = unp.uarray(ld[:,3],ld[:,5])
A1 = unp.uarray(ld[:,4],ld[:,5])

epsr=2.25
vphas=const.c/np.sqrt(epsr)
print("Phasengeschwindigkeit:")
print(vphas)
#vphas = 2e8

L = vphas*delt/2



### Bestimmung der Daempfung ###

print(A0)
print(A1)

Alpha =  (1/(2*L))*unp.log(A0/A1)

daem = ufloat(np.mean(unp.nominal_values(Alpha)),np.std(unp.nominal_values(Alpha)))
print(Alpha)
print(f"alpha    = {daem.nominal_value:.8f}    \u00B1 {daem.std_dev:.8f} ")

### Dielektrizität ###

eps = ((const.c*delt)/(2*L))**2



### Für Tabelle ###

with open("laenge_daempfung.txt", "w") as f:
    for i in range(len(L)):
        f.write(f"{index[i]} & {L[i].nominal_value:.3f} & {L[i].std_dev:.3f} & {Alpha[i].nominal_value:.4f} & {Alpha[i].std_dev:.4f}  \\\\\n")


