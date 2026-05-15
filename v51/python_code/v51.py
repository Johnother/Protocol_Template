import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit

### define custom functions here ###

def polynome(x,a,b):
    return a*x**b


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
ax.scatter(lin1[:,0],lin1[:,1],marker='x',label=r"$A_1=10$")
popt, pcov = curve_fit(polynome,lin1[10:36,0],lin1[10:36,1])
ax.plot(x,polynome(x,*popt))
print("Plateau 1:\nA =",popt[0],'\t','b =',popt[1],'\n')
popt, pcov = curve_fit(polynome,lin1[37:54,0],lin1[37:54,1])
ax.plot(x[500:],polynome(x[500:],*popt))
print("Abfall 1:\nA =",popt[0],'\t','b =',popt[1],'\n')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.legend()
fig.savefig("lin1.pdf")

fig, ax = plt.subplots()
ax.scatter(lin2[:,0],lin2[:,1],color="red",marker='x',label=r"$A_2=6.66$")
popt, pcov = curve_fit(polynome,lin2[13:18,0],lin2[13:18,1])
ax.plot(x,polynome(x,*popt))
print("Plateau 2:\nA =",popt[0],'\t','b =',popt[1],'\n')
popt, pcov = curve_fit(polynome,lin2[19:24,0],lin2[19:24,1])
ax.plot(x[1000:],polynome(x[1000:],*popt))
print("Abfall 2:\nA =",popt[0],'\t','b =',popt[1],'\n')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.legend()
fig.savefig("lin2.pdf")

fig, ax = plt.subplots()
ax.scatter(lin3[:,0],lin3[:,1],color="green",marker='x',label=r"$A_3=2.13$")
popt, pcov = curve_fit(polynome,lin3[12:16,0],lin3[12:16,1])
ax.plot(x,polynome(x,*popt))
print("Plateau 3:\nA =",popt[0],'\t','b =',popt[1],'\n')
popt, pcov = curve_fit(polynome,lin3[17:22,0],lin3[17:22,1])
ax.plot(x[2000:],polynome(x[2000:],*popt))
print("Abfall 3:\nA =",popt[0],'\t','b =',popt[1],'\n')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.legend()
fig.savefig("lin3.pdf")

fig, ax = plt.subplots()
ax.scatter(lin4[:,0],lin4[:,1],color="orange",marker='x',label=r"$A_4=1.47$")
popt, pcov = curve_fit(polynome,lin4[11:15,0],lin4[11:15,1])
ax.plot(x,polynome(x,*popt))
print("Plateau 4:\nA =",popt[0],'\t','b =',popt[1],'\n')
popt, pcov = curve_fit(polynome,lin4[15:19,0],lin4[15:19,1])
ax.plot(x[7000:],polynome(x[7000:],*popt))
print("Abfall 4:\nA =",popt[0],'\t','b =',popt[1],'\n')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.legend()
fig.savefig("lin4.pdf")



# visualize the phase-shift for the lin-amp

fig, ax = plt.subplots()
ax.scatter(lin1[:,0],lin1[:,2],marker='x',label=r"$A_1=10$")
ax.scatter(lin2[:,0],lin2[:,2],marker='x',label=r"$A_2=6.66$")
ax.scatter(lin3[:,0],lin3[:,2],marker='x',label=r"$A_3=2.13$")
ax.scatter(lin4[:,0],lin4[:,2],marker='x',label=r"$A_4=1.47$")

ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.legend()
fig.savefig("lin_phase.pdf")

