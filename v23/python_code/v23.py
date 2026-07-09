import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit, fsolve
from scipy.signal import find_peaks
from scipy.special import sph_harm_y
from uncertainties import ufloat
import uncertainties.unumpy as unp

### functions ###

def theta(alpha):
    return unp.arccos(0.5*unp.cos(alpha) - 0.5)

def lin(x,a,b):
    return a*x+b

### frequenzabhängigkeit ###

### Wasserstoff-Atom ###
H180 = np.genfromtxt("H_180_spektrum.dat",delimiter=" ")
H180[:,0] /= 1000.0

peaks, props = find_peaks(H180[:,1],height=5)
print(peaks)
peaks = np.delete(peaks,(2,8,12))
print(peaks)
fig, ax = plt.subplots(layout="constrained")
ax.plot(H180[:,0],H180[:,1],label=r"Messwerte")
ax.vlines(H180[peaks,0],ymin=0,ymax=26,ls="--",color="red",label=r"Resonanzfrequenzen")
secax = ax.secondary_xaxis("top")
secax.set_xticks(H180[peaks,0])
secax.set_xticklabels(H180[peaks,0],rotation=45,ha="center")
secax.set_xlabel(r"Resonanzfrequenzen [kHz]")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,26)

ax.grid()
ax.legend()

fig.savefig("H_180_Spektrum.pdf")



alpha = np.arange(5,181,5)
amplitude = np.genfromtxt("Winkelmessung/0.dat",delimiter=" ")

for a in alpha:
    file = f"Winkelmessung/{a}.dat"
    data = np.genfromtxt(file,delimiter=" ")
    amplitude = np.c_[amplitude,[round(data[i,1],3) for i in range(len(data[:,1]))]]

resonance = [peaks[0],peaks[1],peaks[2],peaks[5]]
f_res = H180[resonance,0]
print(" Verwendete Resonanzfrequenzen:")
for i in range(len(f_res)):
    print(f"{f_res[i]:.3f}")

resonance = [61,340,597,1089]
print(resonance)
alpha = np.deg2rad(np.arange(0,181,5))
# alpha= np.arange(0,181,5)
theta = theta(alpha)

theta_theory = np.deg2rad(np.linspace(0,361,1000))
phi = np.zeros_like(theta_theory)

print(len(amplitude[resonance[0],:]))
print(amplitude)

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(1,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[0],1:])
ax.plot(theta,amplitude[resonance[0],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^1_0$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_2.31kHz.pdf")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(2,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[1],1:])
ax.plot(theta,amplitude[resonance[1],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^1_0$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_3.705kHz.pdf")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(3,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[2],1:])
ax.plot(theta,amplitude[resonance[2],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^1_0$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_4.99kHz.pdf")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(5,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[3],1:])
ax.plot(theta,amplitude[resonance[3],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^1_0$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_7.45kHz.pdf")



amp_3 = np.genfromtxt("3mm_Aufspaltung.dat",delimiter=" ")
amp_3[:,0] /= 1000.0
amp_9 = np.genfromtxt("9mm_Aufspaltung.dat",delimiter=" ")
amp_9[:,0] /= 1000.0
amp_12 = np.genfromtxt("12mm_Aufspaltung.dat",delimiter=" ")
amp_12[:,0] /= 1000.0

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_3[:,1],height=4)
ax.plot(amp_3[:,0],amp_3[:,1],label=r"Aufspaltung 3mm Ring")
secax = ax.secondary_xaxis("top")
ax.vlines(amp_3[peaks,0],ymin=0,ymax=5.9,ls="--",color="red",label=r"Resonanzfrequenzen")
secax.set_xlabel(r"$f$ [kHz]")
secax.set_xticks(amp_3[peaks,0])
secax.set_xticklabels(amp_3[peaks,0],rotation=45,ha="center")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(1.8,2.6)
ax.set_ylim(0,5.9)
ax.grid()
ax.legend()
fig.savefig("3mm_Aufspaltung.pdf")



fig, ax = plt.subplots(layout="constrained")

peaks9, props = find_peaks(amp_9[:,1],height=5)

ax.plot(amp_9[:,0],amp_9[:,1],label=r"Aufspaltung 9mm Ring")
secax = ax.secondary_xaxis("top")
ax.vlines(amp_9[peaks9,0],ymin=0,ymax=6.5,ls="--",color="red",label=r"Resonanzfrequenzen")
secax.set_xlabel(r"$f$ [kHz]")
secax.set_xticks(amp_9[peaks9,0])
secax.set_xticklabels(amp_9[peaks9,0],rotation=45,ha="center")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(1.8,2.6)
ax.set_ylim(0,6.5)
ax.grid()
ax.legend()
fig.savefig("9mm_Aufspaltung.pdf")






fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_12[:,1],height=3)
ax.plot(amp_12[:,0],amp_12[:,1],label=r"Aufspaltung 12mm Ring")
ax.vlines(amp_12[peaks,0],ymax=6.9,ymin=0,ls="--",color="red",label=r"Resonanzfrequenzen")
secax = ax.secondary_xaxis("top")
secax.set_xlabel(r"$f$ [kHz]")
secax.set_xticks(amp_12[peaks,0])
secax.set_xticklabels(amp_12[peaks,0],rotation=45,ha="center")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(1.8,2.6)
ax.set_ylim(0,6.9)
ax.grid()
ax.legend()
fig.savefig("12mm_Aufspaltung.pdf")



alpha = np.arange(5,181,5)
amplitude = np.genfromtxt("9mm_Winkel/0.dat",delimiter=" ")
for a in alpha:
    file = f"9mm_Winkel/{a}.dat"
    data = np.genfromtxt(file,delimiter=" ")
    amplitude = np.c_[amplitude,data[:,1]]

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
alpha = np.deg2rad(np.arange(0,181,5))
Y = np.abs(sph_harm_y(0,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[peaks9[0],1:])

ax.plot(alpha,amplitude[peaks9[0],1:],"x",label=r"Messwewrte")
ax.plot(theta_theory,Y_N,label=r"$Y^0_0$")
ax.grid(True)
ax.legend()
fig.savefig("Polar_9_1.pdf")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(1,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[peaks9[1],1:])


ax.plot(alpha,amplitude[peaks9[1],1:],"x",label=r"Messwewrte")
ax.plot(theta_theory,Y_N,label=r"$Y^1_0$")
ax.grid(True)
ax.legend()
fig.savefig("Polar_9_2.pdf")



### Zylinder-Messungen ###

amp_50 = np.genfromtxt("50mm_Zylinder/1.dat",delimiter=" ")
for i in range(2,13):
    file = f"50mm_Zylinder/{i}.dat"
    data = np.genfromtxt(file,delimiter=" ")
    amp_50 = np.c_[amp_50,data[:,1]]

amp_50[:,0] /= 1000.0

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_50[:,1],height=5)
peaks = np.delete(peaks,(0,2))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF = []
DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_50mm = {DeltaF[0]} kHz")

ax.plot(amp_50[:,0],amp_50[:,1])
ax.plot(amp_50[peaks,0],amp_50[peaks,1],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_1.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,2],height=5)
# peaks = np.delete(peaks,(0,2))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_100mm = {DeltaF[1]} kHz")

ax.plot(amp_50[:,0],amp_50[:,2])
ax.plot(amp_50[peaks,0],amp_50[peaks,2],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_2.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,3],height=2)
# peaks = np.delete(peaks,(0,2))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_150mm = {DeltaF[2]} kHz")

ax.plot(amp_50[:,0],amp_50[:,3])
ax.plot(amp_50[peaks,0],amp_50[peaks,3],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_3.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,4],height=2.5)
# peaks = np.delete(peaks,(0,2))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_200mm = {DeltaF[3]} kHz")

ax.plot(amp_50[:,0],amp_50[:,4])
ax.plot(amp_50[peaks,0],amp_50[peaks,4],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_4.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,5],height=2.7)
# peaks = np.delete(peaks,(0,2))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_250mm = {DeltaF[4]} kHz")

ax.plot(amp_50[:,0],amp_50[:,5])
ax.plot(amp_50[peaks,0],amp_50[peaks,5],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_5.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,6],height=2)
# peaks = np.delete(peaks,(0,2))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_300mm = {DeltaF[5]} kHz")

ax.plot(amp_50[:,0],amp_50[:,6])
ax.plot(amp_50[peaks,0],amp_50[peaks,6],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_6.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,7],height=2)
peaks = np.delete(peaks,(0))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_350mm = {DeltaF[6]} kHz")

ax.plot(amp_50[:,0],amp_50[:,7])
ax.plot(amp_50[peaks,0],amp_50[peaks,7],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_7.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,8],height=2)
# peaks = np.delete(peaks,(0,2))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_400mm = {DeltaF[7]} kHz")

ax.plot(amp_50[:,0],amp_50[:,8])
ax.plot(amp_50[peaks,0],amp_50[peaks,8],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_8.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,9],height=2)
peaks = np.delete(peaks,(1))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_450mm = {DeltaF[8]} kHz")

ax.plot(amp_50[:,0],amp_50[:,9])
ax.plot(amp_50[peaks,0],amp_50[peaks,9],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_9.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,10],height=1.9)
peaks = np.delete(peaks,(1,2,3))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_500mm = {DeltaF[9]} kHz")

ax.plot(amp_50[:,0],amp_50[:,10])
ax.plot(amp_50[peaks,0],amp_50[peaks,10],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_10.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,11],height=2.2)
peaks = np.delete(peaks,(1,2,3,4))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_550mm = {DeltaF[10]} kHz")

ax.plot(amp_50[:,0],amp_50[:,11])
ax.plot(amp_50[peaks,0],amp_50[peaks,11],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_11.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,12],height=2)
peaks = np.delete(peaks,(0,2,3))
freq = unp.uarray(amp_50[peaks,0],np.ones(len(peaks)))
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_600mm = {DeltaF[11]} kHz")

ax.plot(amp_50[:,0],amp_50[:,12])
ax.plot(amp_50[peaks,0],amp_50[peaks,12],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
fig.savefig("50mm_12.pdf")

for i in range(len(DeltaF)):
    DeltaF[i] = DeltaF[i]**2

F_N = [DeltaF[i].n for i in range(len(DeltaF))]
F_S = [DeltaF[i].s for i in range(len(DeltaF))]

L = np.arange(50,601,50)
x = np.linspace(50,601)

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(1/(L**2),F_N,yerr=F_S,fmt="x",capsize=2,label=r"Messwerte: $\increment f$")
popt, pcov = curve_fit(lin,1/(L**2),F_N,sigma=F_S,absolute_sigma=True)
sigma = np.sqrt(np.diag(pcov))
a = ufloat(popt[0],sigma[0])
b = ufloat(popt[1],sigma[1])
print("--- Fit 50mm ---")
print(f"a = {a}\tb = {b}")
print(f"v = {2*unp.sqrt(a)}\n---------------------------------------------")
ax.plot(1/(x**2),lin(1/(x**2),*popt),label=r"Linearer Fit")
ax.set_xlabel(r"$L^{-2}[\mathrm{mm}^{-2}]$")
ax.set_ylabel(r"$\increment f$ [kHz]")
ax.grid()
ax.legend()
fig.savefig("50mm_dF.pdf")