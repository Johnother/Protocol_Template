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


### frequenzabhängigkeit ###


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
amp_9 = np.genfromtxt("9mm_Aufspaltung.dat",delimiter=" ")
amp_12 = np.genfromtxt("12mm_Aufspaltung.dat",delimiter=" ")

fig, ax = plt.subplots(layout="constrained")
ax.plot(amp_3[:,0]/1000.0,amp_3[:,1],label=r"Aufspalltung 3mm Ring")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
ax.legend()
fig.savefig("3mm_Aufspaltung.pdf")



fig, ax = plt.subplots(layout="constrained")
ax.plot(amp_9[:,0]/1000.0,amp_9[:,1],label=r"Aufspalltung 9mm Ring")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
ax.grid()
ax.legend()
fig.savefig("9mm_Aufspaltung.pdf")

peaks, props = find_peaks(amp_9[:,1],height=5)




fig, ax = plt.subplots(layout="constrained")
ax.plot(amp_12[:,0]/1000.0,amp_12[:,1],label=r"Aufspalltung 12mm Ring")
ax.set_xlabel(r"$f$ [kHz]")
ax.set_ylabel(r"Amplitude [a.u.]")
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
