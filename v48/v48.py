import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
from uncertainties import ufloat, unumpy as unp
from scipy.optimize import curve_fit

# Global constants




#functions

def bmean(T):
    return np.mean(np.diff(T))

def expo(x,a,b):
    return a*np.exp(b*x)


# Measurement for b=2

data2 = np.genfromtxt('b2.txt')
t2 = data2[:,0]
T2 = unp.uarray(data2[:,1],data2[:,2])
I2 = unp.uarray(-data2[:,3],data2[:,4])

bm2 = bmean(T2)

print("bmean for b2:")
print(bm2)

T2 = T2 + 273.15


mask = np.ones(len(T2), dtype=bool)

mask[:5] = False
mask[43:] = False
mask[11:14] = False
mask[16] = False
mask[17] = False
mask[21:39] = False


popt, pcov = curve_fit(expo,unp.nominal_values(T2[mask]),unp.nominal_values(I2[mask]),sigma=unp.std_devs(I2[mask]),absolute_sigma=True,p0=[0.01,0.1])
perr = np.sqrt(np.diag(pcov))
a = ufloat(popt[0], perr[0])
b = ufloat(popt[1], perr[1])


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(unp.nominal_values(T2[mask]),unp.nominal_values(I2[mask]),yerr=unp.std_devs(I2[mask]),fmt="x",color="darkorange",capsize=3,label="Hintergrund")
ax.errorbar(unp.nominal_values(T2[~mask]),unp.nominal_values(I2[~mask]),yerr=unp.std_devs(I2[~mask]),fmt="x",color="grey",capsize=3,label="Messwerte")

x_lin = np.linspace(min(unp.nominal_values(T2[mask])),max(unp.nominal_values(T2[mask])),500)
ax.plot(x_lin, expo(x_lin, *popt), label="Exponentieller Fit")

ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b2.pdf")







# Measurement for b=1.5

data15 = np.genfromtxt('b15.txt')
t15 = data15[:,0]
T15 = unp.uarray(data15[:,1],data15[:,2])
I15 = unp.uarray(-data15[:,3],data15[:,4])

bm15 = bmean(T15)

print("b_mean for b15:")
print(bm15)

T15 = T15 + 273.15


mask = np.ones(len(T15), dtype=bool)

mask[:8] = False
mask[52:] = False
mask[10] = False
mask[12] = False
mask[13] = False
mask[15:18] = False
mask[23:45] = False


popt, pcov = curve_fit(expo,unp.nominal_values(T15[mask]),unp.nominal_values(I15[mask]),sigma=unp.std_devs(I15[mask]),absolute_sigma=True,p0=[0.01,0.1])
perr = np.sqrt(np.diag(pcov))
a = ufloat(popt[0], perr[0])
b = ufloat(popt[1], perr[1])


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(unp.nominal_values(T15[mask]),unp.nominal_values(I15[mask]),yerr=unp.std_devs(I15[mask]),fmt="x",color="darkorange",capsize=3,label="Hintergrund")
ax.errorbar(unp.nominal_values(T15[~mask]),unp.nominal_values(I15[~mask]),yerr=unp.std_devs(I15[~mask]),fmt="x",color="grey",capsize=3,label="Messwerte")

x_lin = np.linspace(min(unp.nominal_values(T15[mask])),max(unp.nominal_values(T15[mask])),500)
ax.plot(x_lin, expo(x_lin, *popt), label="Exponentieller Fit")

ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b15.pdf")
