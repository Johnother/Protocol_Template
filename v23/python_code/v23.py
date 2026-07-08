import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit, fsolve
from scipy.signal import find_peaks
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

resonance = [peaks[0],peaks[1],peaks[2],peaks[5]]
resonance = H180[resonance,0]
print(resonance)

alpha = np.arange(5,181,5)
amplitude = np.genfromtxt("Winkelmessung/0.dat",delimiter=" ")

for a in alpha:
    file = f"Winkelmessung/{a}.dat"
    data = np.genfromtxt(file,delimiter=" ")
    amplitude = np.c_[amplitude,data[:,1]]

fig, ax = plt.subplots(layout="constrained")
ax.plot(amplitude[:,0],amplitude[:,1:],label=r"Messwerte")
ax.grid()
# ax.legend()

fig.savefig("Test_Winkel.pdf")