import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
from uncertainties import ufloat, umath, unumpy as unp
from scipy.optimize import curve_fit

# Global constants

k_b = const.k
e = const.electron_volt

#functions

def bmean(T):
    return np.mean(np.diff(T))

def expo(x,a,b):
    return a*np.exp(b*x)

def expom(x,a,b):
    return a*np.exp(-b/x)

def expou(x,a,b):
    return a*umath.exp(b*x)

def lin(x,a,b):
    return a*x+b

# Measurement for b=2

data2 = np.genfromtxt('b2.txt')
t2 = data2[:,0]
T2 = unp.uarray(data2[:,1],data2[:,2])
I2 = unp.uarray(-data2[:,3],data2[:,4])

bm2 = bmean(T2)

print("bmean for b2:")
print(bm2)

T2 = T2 + 273.15


#find background

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
ax.plot(x_lin, expo(x_lin, *popt),color="darkorchid", label="Exponentieller Fit")

ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$-I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b2_background.pdf")

print("fit parameters:")
print(a)
print(b)




#methode one

mask[:] = True

for i in range(len(T2)):
    I2[i]= (I2[i]-expou(T2[i], a, b))

mask[:19] = False
mask[41:] = False

mask2 = np.ones(len(T2), dtype=bool)

mask2[:21]= False
mask2[29:] = False



#popt, pcov = curve_fit(expo,unp.nominal_values(T2[mask2]),unp.nominal_values(I2[mask2]),sigma=unp.std_devs(I2[mask2]),absolute_sigma=True,p0=[0.01,0.1])
#perr = np.sqrt(np.diag(pcov))
#c = ufloat(popt[0], perr[0])
#d = ufloat(popt[1], perr[1])


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(unp.nominal_values(T2[mask]),unp.nominal_values(I2[mask]),yerr=unp.std_devs(I2[mask]),fmt="x",color="orange",capsize=3,label="cleansed data")
#ax.errorbar(unp.nominal_values(T2[mask2]),unp.nominal_values(I2[mask2]),yerr=unp.std_devs(I2[mask2]),fmt="x",color="darkorange",capsize=3,label="Fit Werte")
#ax.errorbar(unp.nominal_values(T2[~mask]),unp.nominal_values(I2[~mask]),yerr=unp.std_devs(I2[~mask]),fmt="x",color="grey",capsize=3,label="Messwerte")

#x_lin = np.linspace(min(unp.nominal_values(T2[mask2])),max(unp.nominal_values(T2[mask2])),500)
#ax.plot(x_lin, expo(x_lin, *popt),color="darkorchid", label="Exponentieller Fit")


ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$-I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b2.pdf")

#print("fit parameters:")
#print(c)
#print(d)


I2_log = unp.log(abs(I2))


popt, pcov = curve_fit(lin,1/(unp.nominal_values(T2[mask2])),unp.nominal_values(I2_log[mask2]),sigma=unp.std_devs(I2_log[mask2]),absolute_sigma=True,p0=[0.01,0.1])
perr = np.sqrt(np.diag(pcov))
c = ufloat(popt[0], perr[0])
d = ufloat(popt[1], perr[1])


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(1/(unp.nominal_values(T2[mask])),unp.nominal_values(I2_log[mask]),yerr=unp.std_devs(I2_log[mask]),fmt="x",color="grey",capsize=3,label="Hintergrund")
ax.errorbar(1/(unp.nominal_values(T2[mask2])),unp.nominal_values(I2_log[mask2]),yerr=unp.std_devs(I2_log[mask2]),fmt="x",color="darkorange",capsize=3,label="Fit Werte")
#ax.errorbar(unp.nominal_values(T2[~mask]),unp.nominal_values(I2[~mask]),yerr=unp.std_devs(I2[~mask]),fmt="x",color="grey",capsize=3,label="Messwerte")

x_lin = np.linspace(min(unp.nominal_values(T2[mask2])),max(unp.nominal_values(T2[mask2])),500)
ax.plot(1/(x_lin), lin(1/(x_lin), *popt),color="darkorchid", label="Lin Fit")

ax.set_xlabel(r"$1/T$ [1/K]")
ax.set_ylabel(r"$-I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b2_log.pdf")

print("fit parameters:")
print(c)
print(d)
print("activation energy:")
w=-(c*k_b)/e
print(w)


#####################################################################
print("\n")
# Measurement for b=1.5

data15 = np.genfromtxt('b15.txt')
t15 = data15[:,0]
T15 = unp.uarray(data15[:,1],data15[:,2])
I15 = unp.uarray(-data15[:,3],data15[:,4])

bm15 = bmean(T15)

print("b_mean for b15:")
print(bm15)

T15 = T15 + 273.15


#find background

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
ax.plot(x_lin, expo(x_lin, *popt),color="darkorchid", label="Exponentieller Fit")

ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$-I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b15_background.pdf")

print("fit parameters:")
print(a)
print(b)



#methode one

mask[:] = True

for i in range(len(T15)):
    I15[i]= (I15[i]-expou(T15[i], a, b))

mask[:21] = False
mask[45:] = False

mask2 = np.ones(len(T15), dtype=bool)

mask2[:23]= False
mask2[32:] = False



#popt, pcov = curve_fit(expo,unp.nominal_values(T15[mask2]),unp.nominal_values(I15[mask2]),sigma=unp.std_devs(I15[mask2]),absolute_sigma=True,p0=[0.01,0.1])
#perr = np.sqrt(np.diag(pcov))
#c = ufloat(popt[0], perr[0])
#d = ufloat(popt[1], perr[1])


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(unp.nominal_values(T15[mask]),unp.nominal_values(I15[mask]),yerr=unp.std_devs(I15[mask]),fmt="x",color="orange",capsize=3,label="cleansed data")
#ax.errorbar(unp.nominal_values(T15[mask2]),unp.nominal_values(I15[mask2]),yerr=unp.std_devs(I15[mask2]),fmt="x",color="darkorange",capsize=3,label="Fit Werte")
#ax.errorbar(unp.nominal_values(T2[~mask]),unp.nominal_values(I2[~mask]),yerr=unp.std_devs(I2[~mask]),fmt="x",color="grey",capsize=3,label="Messwerte")

#x_lin = np.linspace(min(unp.nominal_values(T15[mask2])),max(unp.nominal_values(T15[mask2])),500)
#ax.plot(x_lin, expo(x_lin, *popt),color="darkorchid", label="Exponentieller Fit")

ax.set_xlabel(r"$T$ [K]")
ax.set_ylabel(r"$-I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b15.pdf")

#print("fit parameters:")
#print(c)
#print(d)

I15_log = unp.log(abs(I15))


popt, pcov = curve_fit(lin,1/(unp.nominal_values(T15[mask2])),unp.nominal_values(I15_log[mask2]),sigma=unp.std_devs(I15_log[mask2]),absolute_sigma=True,p0=[0.01,0.1])
perr = np.sqrt(np.diag(pcov))
c = ufloat(popt[0], perr[0])
d = ufloat(popt[1], perr[1])


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(1/(unp.nominal_values(T15[mask])),unp.nominal_values(I15_log[mask]),yerr=unp.std_devs(I15_log[mask]),fmt="x",color="grey",capsize=3,label="Hintergrund")
ax.errorbar(1/(unp.nominal_values(T15[mask2])),unp.nominal_values(I15_log[mask2]),yerr=unp.std_devs(I15_log[mask2]),fmt="x",color="darkorange",capsize=3,label="Fit Werte")
#ax.errorbar(unp.nominal_values(T2[~mask]),unp.nominal_values(I2[~mask]),yerr=unp.std_devs(I2[~mask]),fmt="x",color="grey",capsize=3,label="Messwerte")

x_lin = np.linspace(min(unp.nominal_values(T15[mask2])),max(unp.nominal_values(T15[mask2])),500)
ax.plot(1/(x_lin), lin(1/(x_lin), *popt),color="darkorchid", label="Lin Fit")

ax.set_xlabel(r"$1/T$ [1/K]")
ax.set_ylabel(r"$-I$ [pA]")


ax.grid()
ax.legend()

fig.savefig("images/b15_log.pdf")

print("fit parameters:")
print(c)
print(d)

print("activation energy:")
w=-(c*k_b)/e
print(w)