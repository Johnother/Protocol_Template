import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit, fsolve
from uncertainties import ufloat
import uncertainties.unumpy as unp

### define custom functions here ###

def polynome(x,a,b):
    return a*x**b

def f_G(A0,A,b):
    return (A0/(A*np.sqrt(2)))**(1/b)



### end of functions ###

# Resistors used for the linear amplifier
R1 = np.array([10000,15000,47000,68000])
R2 = 100000
A = R2/R1
print(A)
# read in data for the four linear amplifier setups
lin1 = np.genfromtxt("linear_1.txt")
lin2 = np.genfromtxt("linear_2.txt")
lin3 = np.genfromtxt("linear_3.txt")
lin4 = np.genfromtxt("linear_4.txt")

lin1[:,1] = 2*lin1[:,1]
lin2[:,1] = 2*lin2[:,1]
lin3[:,1] = 2*lin3[:,1]
lin4[:,1] = 2*lin4[:,1]

# Resistor and capacitor for integrator and differentiator
R = np.array([100000,1000])
C = np.array([100e-9,1e-6])

# read in data for the inegrator and differentiator
int = np.genfromtxt("integrator.txt")
diff = np.genfromtxt("differenzierer.txt")

# read in data for the Schmitt-Trigger

schmitt = np.genfromtxt("schmitt_trigger.txt")

# visualize the lin-amp data in a loglog diagram

x = np.linspace(1,1e6,100000)
print(x)

fig, ax = plt.subplots()
ax.scatter(lin1[10:,0],lin1[10:,1],marker='x',label=r"Messwerte")
ax.scatter(lin1[:10,0],lin1[:10,1],marker='x',color="grey",label=r"Nicht verwendet")
popt1, pcov1 = curve_fit(polynome,lin1[10:36,0],lin1[10:36,1])
perr1 = np.sqrt(np.diag(pcov1))
A0_1 = ufloat(popt1[0],perr1[0])
# print("A0 =",A0_1,"\t","b =",popt1[1],"+/-",perr1[0])
ax.plot(x,polynome(x,*popt1),color="green",label=r"Fit: Plateau")
print("Plateau 1:\nA =",A0_1,'\t','b =',popt1[1],"+/-",perr1[0],'\n')
popt2, pcov2 = curve_fit(polynome,lin1[37:54,0],lin1[37:54,1])
perr2 = np.sqrt(np.diag(pcov2))
A_1 = ufloat(popt2[0],perr2[0])
b_1 = ufloat(popt2[1],perr2[1])
# print("A_lin =",popt2[0],"+/-",perr2[0],"\t","b_lin =",popt2[1],"+/-",perr2[1])
ax.plot(x[400:],polynome(x[400:],*popt2),color="orange",label=r"Fit: abfallender Bereich")
print("Abfall 1:\nA =",A_1,'\t','b =',b_1,'\n')
f1 = f_G(A0_1,A_1,b_1)
print("Grenzfrequenz 1: ",f1)
ax.vlines(f1.n,0,100,ls="--",color="red",label=r"$f_\mathrm{G}=22\pm 6\,\mathrm{kHz}$")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"f [Hz]")
ax.set_ylabel(r"Verstärkung V")
ax.set_ylim(0,100)
ax.grid()
ax.legend()
fig.savefig("lin1.pdf")

fig, ax = plt.subplots()
ax.scatter(lin2[13:,0],lin2[13:,1],marker='x',label=r"Messwerte")
ax.scatter(lin2[:13,0],lin2[:13,1],marker='x',color="grey",label=r"Nicht verwendet")
popt1, pcov1 = curve_fit(polynome,lin2[13:18,0],lin2[13:18,1])
perr1 = np.sqrt(np.diag(pcov1))
A0_2 = ufloat(popt1[0],perr1[0])
# print("A0 =",A0_2,"\t","b =",popt1[1],"+/-",perr1[0])
ax.plot(x,polynome(x,*popt1),color="green",label=r"Fit: Plateau")
print("Plateau 2:\nA =",A0_2,'\t','b =',popt1[1],"+/-",perr1[0],'\n')
popt2, pcov2 = curve_fit(polynome,lin2[19:24,0],lin2[19:24,1])
perr2 = np.sqrt(np.diag(pcov2))
A_2 = ufloat(popt2[0],perr2[0])
b_2 = ufloat(popt2[1],perr2[1])
print("Abfall 2:\nA =",A_2,'\t','b =',b_2,'\n')
ax.plot(x[900:],polynome(x[900:],*popt2),color="orange",label=r"Fit: abfallender Bereich")
# print("Abfall 2:\nA =",popt2[0],'\t','b =',popt2[1],'\n')
f2 = f_G(A0_2,A_2,b_2)
print("Grenzfrequenz 2: ",f2)
ax.vlines(f2.n,0,100,ls="--",color="red",label=r"$f_\mathrm{G}=37\pm 15\,\mathrm{kHz}$")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"f [Hz]")
ax.set_ylabel(r"Verstärkung V")
ax.set_ylim(0,60)
ax.grid()
ax.legend()
fig.savefig("lin2.pdf")

fig, ax = plt.subplots()
ax.scatter(lin3[12:,0],lin3[12:,1],marker='x',label=r"Messwerte")
ax.scatter(lin3[:12,0],lin3[:12,1],marker='x',color="grey",label=r"Nicht verwendet")
popt1, pcov1 = curve_fit(polynome,lin3[12:16,0],lin3[12:16,1])
perr1 = np.sqrt(np.diag(pcov1))
A0_3 = ufloat(popt1[0],perr1[0])
ax.plot(x,polynome(x,*popt1),color="green",label=r"Fit: Plateau")
print("Plateau 3:\nA =",A0_3,'\t','b =',popt1[1],"+/-",perr1[0],'\n')
popt2, pcov2 = curve_fit(polynome,lin3[17:22,0],lin3[17:22,1])
perr2 = np.sqrt(np.diag(pcov2))
A_3 = ufloat(popt2[0],perr2[0])
b_3 = ufloat(popt2[1],perr2[1])
ax.plot(x[2000:],polynome(x[2000:],*popt2),color="orange",label=r"Fit: abfallender Bereich")
print("Abfall 3:\nA =",A_3,'\t','b =',b_3,'\n')
f3 = f_G(A0_3,A_3,b_3)
print("Grenzfrequenz 3: ",f3)
ax.vlines(f3.n,0,100,ls="--",color="red",label=r"$f_\mathrm{G}=120\pm 40\,\mathrm{kHz}$")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"f [Hz]")
ax.set_ylabel(r"Verstärkung V")
ax.set_ylim(0,20)
ax.grid()
ax.legend()
fig.savefig("lin3.pdf")

fig, ax = plt.subplots()
ax.scatter(lin4[11:,0],lin4[11:,1],marker='x',label=r"Messwerte")
ax.scatter(lin4[:11,0],lin4[:11,1],marker='x',color="grey",label=r"Nicht verwendet")
popt1, pcov1 = curve_fit(polynome,lin4[11:15,0],lin4[11:15,1])
perr1 = np.sqrt(np.diag(pcov1))
A0_4 = ufloat(popt1[0],perr1[0])
ax.plot(x,polynome(x,*popt1),color="green",label=r"Fit: Plateau")
print("Plateau 4:\nA =",A0_4,'\t','b =',popt1[1],"+/-",perr1[0],'\n')
popt2, pcov2 = curve_fit(polynome,lin4[15:19,0],lin4[15:19,1])
perr2 = np.sqrt(np.diag(pcov2))
A_4 = ufloat(popt2[0],perr2[0])
b_4 = ufloat(popt2[1],perr2[1])
ax.plot(x[6000:],polynome(x[6000:],*popt2),color="orange",label=r"Fit: abfallender Bereich")
print("Abfall 4:\nA =",A_4,'\t','b =',b_4,'\n')
f4 = f_G(A0_4,A_4,b_4)
print("Grenzfrequenz 4: ",f4)
ax.vlines(f4.n,0,100,ls="--",color="red",label=r"$f_\mathrm{G}=170\pm 50\,\mathrm{kHz}$")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r"f [Hz]")
ax.set_ylabel(r"Verstärkung V")
ax.set_ylim(0,10)
ax.grid()
ax.legend()
fig.savefig("lin4.pdf")



# visualize the phase-shift for the lin-amp
lin1[:,2] = lin1[:,2]*lin1[:,0]*360/1000
lin2[:,2] = lin2[:,2]*lin2[:,0]*360/1000
lin3[:,2] = lin3[:,2]*lin3[:,0]*360/1000
lin4[:,2] = lin4[:,2]*lin4[:,0]*360/1000

fig, ax = plt.subplots(2,2,layout="constrained")
ax[0,0].scatter(lin1[:,0],lin1[:,2],marker='x',label=r"$A_1=10$")
ax[0,1].scatter(lin2[:,0],lin2[:,2],marker='x',color="orange",label=r"$A_2=6.66$")
ax[1,0].scatter(lin3[:,0],lin3[:,2],marker='x',color="green",label=r"$A_3=2.13$")
ax[1,1].scatter(lin4[:,0],lin4[:,2],marker='x',color="red",label=r"$A_4=1.47$")

ax[0,0].set_xscale('log')
ax[0,1].set_xscale('log')
ax[1,0].set_xscale('log')
ax[1,1].set_xscale('log')

ax[0,0].set_xlabel(r"f [Hz]")
ax[0,1].set_xlabel(r"f [Hz]")
ax[1,0].set_xlabel(r"f [Hz]")
ax[1,1].set_xlabel(r"f [Hz]")

ax[0,0].set_ylabel(r"\phi [°]")
ax[0,1].set_ylabel(r"\phi [°]")
ax[1,0].set_ylabel(r"\phi [°]")
ax[1,1].set_ylabel(r"\phi [°]")

ax[0,0].grid()
ax[0,1].grid()
ax[1,0].grid()
ax[1,1].grid()

ax[0,0].legend()
ax[0,1].legend()
ax[1,0].legend()
ax[1,1].legend()


fig.savefig("lin_phase.pdf")

