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
ax.errorbar(omega,unp.nominal_values(L),yerr=unp.std_devs(L),fmt="x",color="darkorange",capsize=3,label="Induktivitätsbelag")


#x_fit=np.linspace(np.min(x),np.max(x),1000)

#ax.plot(x_fit,gauss(x_fit, *popt),label="Gaußförmiger-Fit")

ax.set_xlabel(r"$\omega$ [Hz]")
ax.set_ylabel(r"$L$ [H]")

#ax.set_xlim(-16,16)
#ax.set_ylim(0,8)

ax.grid()
ax.legend()

fig.savefig("../content/images/ind.pdf")

fig, ax = plt.subplots(layout="constrained")

ax.errorbar(omega,unp.nominal_values(R),yerr=unp.std_devs(R),fmt="x",color="yellowgreen",capsize=3,label="Ohmscher Belag")


#x_fit=np.linspace(np.min(x),np.max(x),1000)

#ax.plot(x_fit,gauss(x_fit, *popt),label="Gaußförmiger-Fit")

ax.set_xlabel(r"$\omega$ [Hz]")
ax.set_ylabel(r"$R$ [$\Omega$]")

#ax.set_xlim(-16,16)
#ax.set_ylim(0,8)

ax.grid()
ax.legend()

fig.savefig("../content/images/ohm.pdf")

fig, ax = plt.subplots(layout="constrained")

ax.errorbar(omega,unp.nominal_values(Z),yerr=unp.std_devs(Z),fmt="x",color="c",capsize=3,label="Impedanzbetrag")

#x_fit=np.linspace(np.min(x),np.max(x),1000)

#ax.plot(x_fit,gauss(x_fit, *popt),label="Gaußförmiger-Fit")

ax.set_xlabel(r"$\omega$ [Hz]")
ax.set_ylabel(r"$Z$ [$\Omega$]")

#ax.set_xlim(-16,16)
#ax.set_ylim(0,8)

ax.grid()
ax.legend()

fig.savefig("../content/images/imp.pdf")


fig, ax = plt.subplots(layout="constrained")

ax.errorbar(omega,unp.nominal_values(C),yerr=unp.std_devs(C),fmt="x",color="crimson",capsize=3,label="Kapatzitätsbelag")

#x_fit=np.linspace(np.min(x),np.max(x),1000)

#ax.plot(x_fit,gauss(x_fit, *popt),label="Gaußförmiger-Fit")

ax.set_xlabel(r"$\omega$ [Hz]")
ax.set_ylabel(r"$C$ [F]")

#ax.set_xlim(-16,16)
#ax.set_ylim(0,8)

ax.grid()
ax.legend()

fig.savefig("../content/images/kap.pdf")

#### Bestimmung der Laenge ###

ld = np.genfromtxt("len_daem.txt")



index = ld[:,0]
delt = unp.uarray(ld[:,1]*10**(-9),ld[:,2]*10**(-9))
A0 = unp.uarray(ld[:,3],ld[:,5])
A1 = unp.uarray(ld[:,4],ld[:,5])

vphas = 2e8

L = vphas*delt/2



### Bestimmung der Daempfung ###

print(A0)
print(A1)

Alpha =  (1/L)*unp.log(A0/A1)

daem = ufloat(np.mean(unp.nominal_values(Alpha)),np.std(unp.nominal_values(Alpha)))
print(Alpha)
print(f"alpha    = {daem.nominal_value:.8f}    \u00B1 {daem.std_dev:.8f} ")

### Dielektrizität ###

eps = ((const.c*delt)/(2*L))**2



### Für Tabelle ###

with open("laenge_daempfung.txt", "w") as f:
    for i in range(len(L)):
        f.write(f"{index[i]} & {L[i].nominal_value} & {L[i].std_dev} & {Alpha[i].nominal_value:.5f} & {Alpha[i].std_dev:.5f} & {eps[i].nominal_value:.10f} & {eps[i].std_dev:.10f} \\\\\n")






















########

#p0 = [np.max(I), x[np.argmax(I)], 2.0]
#
#
#popt, pcov = curve_fit(gauss,x,I,p0=p0,sigma=I_err,absolute_sigma=True)
#
#A, x0, w= popt
#A_err, x0_err, w_err = np.sqrt(np.diag(pcov))
#
#
#print("TEM00-Fit:")
#print(f"A     = {A:.3f} ± {A_err:.3f}")
#print(f"x0    = {x0:.3f} ± {x0_err:.3f} mm")
#print(f"w = {w:.3f} ± {w_err:.3f} mm")
##print(f"C     = {C:.3f} ± {C_err:.3f}")
#
#fig, ax = plt.subplots(layout="constrained")
#ax.errorbar(x,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")
#
#x_fit=np.linspace(np.min(x),np.max(x),1000)
#
#ax.plot(x_fit,gauss(x_fit, *popt),label="Gaußförmiger-Fit")
#
#ax.set_xlabel(r"$x$ [mm]")
#ax.set_ylabel(r"$I$ [\mu A]")
#
#ax.set_xlim(-16,16)
#ax.set_ylim(0,8)
#
#ax.grid()
#ax.legend()
#
#fig.savefig("../content/images/tem00.pdf")


