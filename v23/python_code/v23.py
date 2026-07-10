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

plt.rcParams.update({'font.size': 14})

peaks, props = find_peaks(H180[:,1],height=5)
peaks = np.delete(peaks,(2,8,12))
peaks = np.delete(peaks,(4,6,8))
peaks[7] = 2150
fig, ax = plt.subplots(layout="constrained")
ax.plot(H180[:,0],H180[:,1],label=r"Messwerte")
ax.vlines(H180[peaks,0],ymin=0,ymax=26,ls="--",lw=1,color="red",label=r"Resonanzfrequenzen")
secax = ax.secondary_xaxis("top")
secax.set_xticks(H180[peaks,0])
secax.set_xticklabels(H180[peaks,0],rotation=45,ha="center")
secax.set_xlabel(r"$f$ [kHz]")
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
print("--- Verwendete Resonanzfrequenzen ---")
for i in range(len(f_res)):
    print(f"{f_res[i]:.3f}")
print("----------------------------------------------------------------")

resonance = [61,340,597,1089]
alpha = np.deg2rad(np.arange(0,181,5))
# alpha= np.arange(0,181,5)
theta = theta(alpha)

theta_theory = np.deg2rad(np.linspace(0,361,1000))
phi = np.zeros_like(theta_theory)

print(len(amplitude[resonance[0],:]))

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(1,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[0],1:])
ax.plot(theta,amplitude[resonance[0],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^0_1$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_2.31kHz.pdf")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(2,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[1],1:])
ax.plot(theta,amplitude[resonance[1],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^0_2$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_3.705kHz.pdf")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(3,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[2],1:])
ax.plot(theta,amplitude[resonance[2],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^0_3$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_4.99kHz.pdf")

fig, ax = plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
Y = np.abs(sph_harm_y(5,0,theta_theory,phi))
Y_N = Y/np.max(Y)*np.max(amplitude[resonance[3],1:])
ax.plot(theta,amplitude[resonance[3],1:],"x",label=r"Messwerte")
# ax.plot(theta_theory,1/2*(3*np.cos(theta_theory)**2-1))
ax.plot(theta_theory,Y_N,label=r"$Y^0_5$")
ax.grid(True)
ax.legend()

fig.savefig("Polar_7.45kHz.pdf")



amp_3 = np.genfromtxt("3mm_Aufspaltung.dat",delimiter=" ")
amp_3[:,0] /= 1000.0
amp_9 = np.genfromtxt("9mm_Aufspaltung.dat",delimiter=" ")
amp_9[:,0] /= 1000.0
amp_12 = np.genfromtxt("12mm_Aufspaltung.dat",delimiter=" ")
amp_12[:,0] /= 1000.0





fig, ax = plt.subplots(figsize=(5,5),layout="constrained")
peaks3, props = find_peaks(amp_3[:,1],height=4)
ax.plot(amp_3[:,0],amp_3[:,1],label=r"Aufspaltung 3mm Ring")
secax = ax.secondary_xaxis("top")
ax.vlines(amp_3[peaks3,0],ymin=0,ymax=5.9,ls="--",lw=1,color="red",label=r"Resonanzfrequenzen")
secax.set_xlabel(r"$f$ [kHz]")
secax.set_xticks(amp_3[peaks3,0])
secax.set_xticklabels(amp_3[peaks3,0],rotation=45,ha="center")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(1.8,2.6)
ax.set_ylim(0,5.9)
ax.grid()
ax.legend()
fig.savefig("3mm_Aufspaltung.pdf")



fig, ax = plt.subplots(figsize=(5,5),layout="constrained")

peaks9, props = find_peaks(amp_9[:,1],height=5)
ax.plot(amp_9[:,0],amp_9[:,1],label=r"Aufspaltung 9mm Ring")
secax = ax.secondary_xaxis("top")
ax.vlines(amp_9[peaks9,0],ymin=0,ymax=6.5,ls="--",lw=1,color="red",label=r"Resonanzfrequenzen")
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



fig, ax = plt.subplots(figsize=(5,5),layout="constrained")

peaks12, props = find_peaks(amp_12[:,1],height=3)
ax.plot(amp_12[:,0],amp_12[:,1],label=r"Aufspaltung 12mm Ring")
ax.vlines(amp_12[peaks12,0],ymax=6.9,ymin=0,ls="--",lw=1,color="red",label=r"Resonanzfrequenzen")
secax = ax.secondary_xaxis("top")
secax.set_xlabel(r"$f$ [kHz]")
secax.set_xticks(amp_12[peaks12,0])
secax.set_xticklabels(amp_12[peaks12,0],rotation=45,ha="center")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(1.8,2.6)
ax.set_ylim(0,6.9)
ax.grid()
ax.legend()
fig.savefig("12mm_Aufspaltung.pdf")

dF = [ (ufloat(amp_3[peaks3[1],0],0.01) - ufloat(amp_3[peaks3[0],0],0.01)) , (ufloat(amp_9[peaks9[1],0],0.01) - ufloat(amp_9[peaks9[0],0],0.01)) , (ufloat(amp_12[peaks12[1],0],0.01) - ufloat(amp_12[peaks12[0],0],0.01)) ]
print(dF)
F_N = [dF[i].n*1000.0 for i in range(len(dF))]
F_S = [dF[i].s*1000.0 for i in range(len(dF))]
d = np.linspace(0,50)
fig, ax = plt.subplots(layout="constrained")
ax.errorbar([3,9,12],F_N,yerr=F_S,fmt="x",capsize=2,label=r"Messwerte $\increment f$")
popt, pcov = curve_fit(lin,[3,9,12,],F_N,sigma=F_S,absolute_sigma=True)
ax.plot(d,lin(d,*popt),label=r"Linearer Fit")
ax.set_xlabel(r"$d$ [mm]")
ax.set_ylabel(r"$\increment f$ [Hz]")
ax.set_xlim(0,20)
ax.set_ylim(0,300)
ax.grid()
ax.legend()
fig.savefig("Ring_fit.pdf")


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
ax.plot(theta_theory,Y_N,label=r"$Y^0_1$")
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
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
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
ax.set_xlim(0,12)
ax.set_ylim(0,24.9)
ax.grid()
fig.savefig("50mm_1.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,2],height=5)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_100mm = {DeltaF[1]} kHz")

ax.plot(amp_50[:,0],amp_50[:,2])
ax.plot(amp_50[peaks,0],amp_50[peaks,2],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,29)
ax.grid()
fig.savefig("50mm_2.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,3],height=2)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_150mm = {DeltaF[2]} kHz")

ax.plot(amp_50[:,0],amp_50[:,3])
ax.plot(amp_50[peaks,0],amp_50[peaks,3],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,15.9)
ax.grid()
fig.savefig("50mm_3.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,4],height=2.5)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_200mm = {DeltaF[3]} kHz")

ax.plot(amp_50[:,0],amp_50[:,4])
ax.plot(amp_50[peaks,0],amp_50[peaks,4],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,19.9)
ax.grid()
fig.savefig("50mm_4.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,5],height=2.7)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_250mm = {DeltaF[4]} kHz")

ax.plot(amp_50[:,0],amp_50[:,5])
ax.plot(amp_50[peaks,0],amp_50[peaks,5],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,13.9)
ax.grid()
fig.savefig("50mm_5.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,6],height=2)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_300mm = {DeltaF[5]} kHz")

ax.plot(amp_50[:,0],amp_50[:,6])
ax.plot(amp_50[peaks,0],amp_50[peaks,6],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,22.5)
ax.grid()
fig.savefig("50mm_6.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,7],height=2)
peaks = np.delete(peaks,(0))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_350mm = {DeltaF[6]} kHz")

ax.plot(amp_50[:,0],amp_50[:,7])
ax.plot(amp_50[peaks,0],amp_50[peaks,7],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,27.5)
ax.grid()
fig.savefig("50mm_7.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,8],height=2)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_400mm = {DeltaF[7]} kHz")

ax.plot(amp_50[:,0],amp_50[:,8])
ax.plot(amp_50[peaks,0],amp_50[peaks,8],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,22.5)
ax.grid()
fig.savefig("50mm_8.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,9],height=2)
peaks = np.delete(peaks,(1))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_450mm = {DeltaF[8]} kHz")

ax.plot(amp_50[:,0],amp_50[:,9])
ax.plot(amp_50[peaks,0],amp_50[peaks,9],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,15.9)
ax.grid()
fig.savefig("50mm_9.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,10],height=1.9)
peaks = np.delete(peaks,(1,2,3))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_500mm = {DeltaF[9]} kHz")

ax.plot(amp_50[:,0],amp_50[:,10])
ax.plot(amp_50[peaks,0],amp_50[peaks,10],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,11.9)
ax.grid()
fig.savefig("50mm_10.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,11],height=2.2)
peaks = np.delete(peaks,(1,2,3,4))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_550mm = {DeltaF[10]} kHz")

ax.plot(amp_50[:,0],amp_50[:,11])
ax.plot(amp_50[peaks,0],amp_50[peaks,11],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,13)
ax.grid()
fig.savefig("50mm_11.pdf")

fig, ax = plt.subplots(layout="constrained")

peaks, props = find_peaks(amp_50[:,12],height=2)
peaks = np.delete(peaks,(0,2,3))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_50[peaks,0],ferr)

dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_600mm = {DeltaF[11]} kHz")

ax.plot(amp_50[:,0],amp_50[:,12])
ax.plot(amp_50[peaks,0],amp_50[peaks,12],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,15.9)
ax.grid()
fig.savefig("50mm_12.pdf")

# for i in range(len(DeltaF)):
#     DeltaF[i] = DeltaF[i]**2

F_N = [DeltaF[i].n for i in range(len(DeltaF))]
F_S = [DeltaF[i].s for i in range(len(DeltaF))]

L = np.arange(50,601,50)
x = np.linspace(50,601)

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(1/L,F_N,yerr=F_S,fmt="x",capsize=2,label=r"Messwerte: $\increment f$")
popt, pcov = curve_fit(lin,1/L,F_N,sigma=F_S,absolute_sigma=True)
sigma = np.sqrt(np.diag(pcov))
a = ufloat(popt[0],sigma[0])
b = ufloat(popt[1],sigma[1])
print("--- Fit 50mm ---")
print(f"a = {a}\tb = {b}")
print(f"v = {2*a}\n---------------------------------------------")
ax.plot(1/x,lin(1/x,*popt),label=r"Linearer Fit")
ax.set_xlabel(r"$L^{-1}\,[\mathrm{mm}^{-1}]$")
ax.set_ylabel(r"$\increment f\,[\mathrm{kHz}]$")
ax.ticklabel_format(axis="x",style="sci",scilimits=(-3,-3))
ax.grid()
ax.legend()
fig.savefig("50mm_dF.pdf")

### 75mm Zylinder ###

amp_75 = np.genfromtxt("75mm_Zylinder/1.dat",delimiter=" ")
for i in range(2,9):
    file = f"75mm_Zylinder/{i}.dat"
    data = np.genfromtxt(file,delimiter=" ")
    amp_75 = np.c_[amp_75,data[:,1]]

amp_75[:,0] /= 1000.0

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,1],height=4)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF = []
DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_75mm = {DeltaF[0]} kHz")

ax.plot(amp_75[:,0],amp_75[:,1])
ax.plot(amp_75[peaks,0],amp_75[peaks,1],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,19.9)
ax.grid()
fig.savefig("75mm_1.pdf")

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,2],height=2)
peaks = np.delete(peaks,(0))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_150mm = {DeltaF[1]} kHz")

ax.plot(amp_75[:,0],amp_75[:,2])
ax.plot(amp_75[peaks,0],amp_75[peaks,2],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,15.9)
ax.grid()
fig.savefig("75mm_2.pdf")

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,3],height=2.5)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_225mm = {DeltaF[2]} kHz")

ax.plot(amp_75[:,0],amp_75[:,3])
ax.plot(amp_75[peaks,0],amp_75[peaks,3],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,17.4)
ax.grid()
fig.savefig("75mm_3.pdf")

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,4],height=2.5)
peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_300mm = {DeltaF[3]} kHz")

ax.plot(amp_75[:,0],amp_75[:,4])
ax.plot(amp_75[peaks,0],amp_75[peaks,4],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,22.5)
ax.grid()
fig.savefig("75mm_4.pdf")

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,5],height=2.5)
peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_375mm = {DeltaF[4]} kHz")

ax.plot(amp_75[:,0],amp_75[:,5])
ax.plot(amp_75[peaks,0],amp_75[peaks,5],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,24.9)
ax.grid()
fig.savefig("75mm_5.pdf")

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,6],height=2)
peaks = np.delete(peaks,(0,2,3))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_450mm = {DeltaF[5]} kHz")

ax.plot(amp_75[:,0],amp_75[:,6])
ax.plot(amp_75[peaks,0],amp_75[peaks,6],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,15)
ax.grid()
fig.savefig("75mm_6.pdf")

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,7],height=2)
peaks = np.delete(peaks,(1,2,3,4,5))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_525mm = {DeltaF[6]} kHz")

ax.plot(amp_75[:,0],amp_75[:,7])
ax.plot(amp_75[peaks,0],amp_75[peaks,7],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,11.9)
ax.grid()
fig.savefig("75mm_7.pdf")

fig, ax = plt.subplots(layout="constrained")
peaks, props = find_peaks(amp_75[:,8],height=2)
# peaks = np.delete(peaks,(0,2))
ferr = np.empty(len(peaks))
ferr.fill(0.01)
freq = unp.uarray(amp_75[peaks,0],ferr)
dF = []
for i in range(len(peaks)-1):
    dF.append(freq[i+1]-freq[i])

DeltaF.append(sum(ufloat(dF[i].n,dF[i].s) for i in range(len(peaks)-1))/(len(peaks)-1))
print(f"DeltaF_600mm = {DeltaF[7]} kHz")

ax.plot(amp_75[:,0],amp_75[:,8])
ax.plot(amp_75[peaks,0],amp_75[peaks,8],"rx")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.set_xlim(0,12)
ax.set_ylim(0,15.9)
ax.grid()
fig.savefig("75mm_8.pdf")


# for i in range(len(DeltaF)):
#     DeltaF[i] = DeltaF[i]**2

F_N = [DeltaF[i].n for i in range(len(DeltaF))]
F_S = [DeltaF[i].s for i in range(len(DeltaF))]

L = np.arange(75,601,75)
x = np.linspace(75,601)

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(1/L,F_N,yerr=F_S,fmt="x",capsize=2,label=r"Messwerte: $\increment f^2$")
popt, pcov = curve_fit(lin,1/L,F_N,sigma=F_S,absolute_sigma=True)
sigma = np.sqrt(np.diag(pcov))
a = ufloat(popt[0],sigma[0])
b = ufloat(popt[1],sigma[1])
print("--- Fit 75mm ---")
print(f"a = {a}\tb = {b}")
print(f"v = {2*a}\n---------------------------------------------")
ax.plot(1/x,lin(1/x,*popt),label=r"Linearer Fit")
ax.set_xlabel(r"$L^{-1}\,[\mathrm{mm}^{-1}]$")
ax.set_ylabel(r"$\increment f\,[\mathrm{kHz}]$")
ax.ticklabel_format(axis="x",style="sci",scilimits=(-3,-3))
ax.grid()
ax.legend()
fig.savefig("75mm_dF.pdf")