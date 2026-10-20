import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat, nominal_value, std_dev
from scipy.optimize import curve_fit

def Rutherford(theta,A,B):
    return (A**2/B**2)*1/(np.sin((theta)/2)**4)

# Calibration measurement to find the center of the beam

data = np.genfromtxt('calibration.csv',delimiter=',',unpack=True)
calibration = {"theta": data[0], "N": data[1], "p":data[2]}

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(calibration["theta"],calibration["N"],xerr=0.1,yerr=np.sqrt(calibration["N"]),fmt=".",capsize=2,label=r"Calibration measurement")
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

ax.errorbar(Au2["theta"],[i.n for i in N],xerr=np.deg2rad(0.1),yerr=[i.s for i in N],fmt="x",capsize=2)
popt, pcov = curve_fit(Rutherford,Au2["theta"][2:],[i.n for i in N[2:]],p0=[79,1e-10],sigma=[i.s for i in N[2:]],absolute_sigma=True)
print(popt)
x = np.linspace(1e-6,np.pi/2,1000)
ax.plot(x,Rutherford(x,*popt))
ax.set_ylim(0,15)
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

ax.errorbar(Bi1["theta"],[i.n for i in N],xerr=np.deg2rad(0.1),yerr=[i.s for i in N],fmt="x",capsize=2)
popt, pcov = curve_fit(Rutherford,Bi1["theta"][:],[i.n for i in N[:]],p0=[79,1e-10],sigma=[i.s for i in N[:]],absolute_sigma=True)
print(popt)
x = np.linspace(1e-6,np.pi/2,1000)
ax.plot(x,Rutherford(x,*popt))
ax.set_ylim(0,15)
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

ax.errorbar(Al3["theta"],[i.n for i in N],xerr=np.deg2rad(0.1),yerr=[i.s for i in N],fmt="x",capsize=2)
popt, pcov = curve_fit(Rutherford,Al3["theta"][:],[i.n for i in N[:]],p0=[79,1e-10],sigma=[i.s for i in N[:]],absolute_sigma=True)
print(popt)
x = np.linspace(1e-6,np.pi/2,1000)
ax.plot(x,Rutherford(x,*popt))
ax.set_ylim(0,15)
fig.savefig("Al3.pdf")