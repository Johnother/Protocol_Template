import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
from uncertainties import ufloat, unumpy as unp
from scipy.optimize import curve_fit

# Global constants




#functions

def bmean(T):
    return np.mean(np.diff(T))


# Measurement for b=2

data2 = np.genfromtxt('b2.txt')
t2 = data2[:,0]
T2 = unp.uarray(data2[:,1],data2[:,2])
I2 = unp.uarray(data2[:,3],data2[:,4])

bm2 = bmean(T2)

print("bmean for b2:")
print(bm2)

T2 = T2 + 273.15

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(unp.nominal_values(T2),unp.nominal_values(I2),yerr=unp.std_devs(I2),fmt="x",capsize=3,label="Messwerte")


ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b2.pdf")







# Measurement for b=1.5

data15 = np.genfromtxt('b15.txt')
t15 = data15[:,0]
T15 = unp.uarray(data15[:,1],data15[:,2])
I15 = unp.uarray(data15[:,3],data15[:,4])

bm15 = bmean(T15)

print("b_mean for b15:")
print(bm15)

T15 = T15 + 273.15

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(unp.nominal_values(T15),unp.nominal_values(I15),yerr=unp.std_devs(I15),fmt="x",capsize=3,label="Messwerte")


ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b15.pdf")
