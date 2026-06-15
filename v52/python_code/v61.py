import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit, fsolve
from uncertainties import ufloat
import uncertainties.unumpy as unp

### functions ###



def gauss(x,A,x0,w):
    return A*np.exp(-2*(x-x0)**2/(w**2))

def gaushermite(x,A,x0,w):
    return A*2*((x-x0)/(w))**2*np.exp(-2*(x-x0)**2/w**2)

def lin(x,a,b):
    return a*x+b

def quad(x,a,b,c):
    return a*x**2+b*x+c

def cos(phi,A,phi0):
    return  A*np.cos(phi-phi0)**2

dunkelAmp=0.009
dunkelWat=-0.037

####theoretische Stabilitaetsbed. ####

r=1.4

fig, ax = plt.subplots(layout="constrained")
#ax.errorbar(x,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

x_fit=np.linspace(-1,3,1000)

ax.plot(x_fit,quad(x_fit, (1/(r**2)), -2/r, 1), color="crimson",label="Stabilitätbedinung konkav-konkav")
ax.plot(x_fit,lin(x_fit, -1/r, 1),color="yellowgreen",label="Stabilitätbedinung plan-konkav")

ax.set_xlabel(r"$L$ [m]")
ax.set_ylabel(r"$g_1 g_2$ [1]")

ax.set_xlim(0,2.5)
ax.set_ylim(0,1.5)

ax.grid()
ax.legend()

fig.savefig("../content/images/theostab.pdf")




#####fit für tem00 mode####

tem00 = np.genfromtxt("TEM00.txt")
tem01 = np.genfromtxt("TEM01.txt")



x = tem00[:,0]
I = (tem00[:,1]-dunkelAmp)
I_err = tem00[:,2]



p0 = [np.max(I), x[np.argmax(I)], 2.0]


popt, pcov = curve_fit(gauss,x,I,p0=p0,sigma=I_err,absolute_sigma=True)

A, x0, w= popt
A_err, x0_err, w_err = np.sqrt(np.diag(pcov))


print("TEM00-Fit:")
print(f"A     = {A:.3f} ± {A_err:.3f}")
print(f"x0    = {x0:.3f} ± {x0_err:.3f} mm")
print(f"w = {w:.3f} ± {w_err:.3f} mm")
#print(f"C     = {C:.3f} ± {C_err:.3f}")

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(x,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

x_fit=np.linspace(np.min(x),np.max(x),1000)

ax.plot(x_fit,gauss(x_fit, *popt),label="Gaußförmiger-Fit")

ax.set_xlabel(r"$x$ [mm]")
ax.set_ylabel(r"$I$ [\mu A]")

ax.set_xlim(-16,16)
ax.set_ylim(0,8)

ax.grid()
ax.legend()

fig.savefig("../content/images/tem00.pdf")


######Fit für TEM01 mode######

x = tem01[:,0]
I = (tem01[:,1]-dunkelAmp)
I_err = tem01[:,2]


p0 = [np.max(I),0.5,5]

popt, pcov = curve_fit(gaushermite,x,I,p0=p0,sigma=I_err,absolute_sigma=True)

I0, x0, w= popt
I0_err, x0_err, w_err = np.sqrt(np.diag(pcov))

print("TEM01-Fit:")
print(f"I0 = {I0:.3f} \u00B1  {I0_err:.3f}")
print(f"x0 = {x0:.3f} \u00B1  {x0_err:.3f} mm")
print(f"w  = {w:.3f}  \u00B1 {w_err:.3f} mm")
#print(f"C  = {C:.3f}  \u00B1 {C_err:.3f}")


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(x,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

x_fit=np.linspace(np.min(x),np.max(x),1000)

ax.plot(x_fit,gaushermite(x_fit,*popt),label=r"Fit für $\mathrm{TEM}_{01}$-Mode")

ax.set_xlabel(r"$x$ [mm]")
ax.set_ylabel(r"$I$ [\mu A]")

ax.set_xlim(-16,16)
ax.set_ylim(0,2)

ax.grid()
ax.legend()

fig.savefig("../content/images/tem01.pdf")


####stabilitätsbedingung konkav konkav####

konkon = np.genfromtxt("konkon.txt")


L = konkon[:,0]
L_err = konkon[:,1]
I = (konkon[:,2]-dunkelWat)
I_err = konkon[:,3]


popt,pcov=curve_fit(quad,L,I,sigma=I_err,absolute_sigma=True)

a, b, c = popt
a_err, b_err, c_err = np.sqrt(np.diag(pcov))

print("konkon-Fit:")
print(f"a = {a:.4f} \u00B1 {a_err:.4f}")
print(f"b = {b:.4f} \u00B1 {b_err:.4f}")
print(f"c = {c:.4f} \u00B1 {c_err:.4f}")

fig, ax = plt.subplots(layout="constrained")
ax.errorbar(L,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

x_fit=np.linspace(np.min(x),np.max(x),1000)

ax.plot(x_fit,quad(x_fit, *popt),label="Quadratischer Fit")

ax.set_xlabel(r"$L$ [m]")
ax.set_ylabel(r"$I$ [mW]")

ax.set_xlim(np.min(L)-0.1,np.max(L)+0.1)
ax.set_ylim(0,7)

ax.grid()
ax.legend()

fig.savefig("../content/images/konkon.pdf")

####stabilitätsbedingung plan konkav####

plankon = np.genfromtxt("plankon.txt")


L = plankon[:,0]
L_err = plankon[:,1]
I = (plankon[:,2]-dunkelWat)
I_err = plankon[:,3]

popt,pcov=curve_fit(lin,L,I,sigma=I_err,absolute_sigma=True)

a, b = popt
a_err, b_err = np.sqrt(np.diag(pcov))

print("plankon-Fit:")
print(f"a = {a:.4f} \u00B1 {a_err:.4f}")
print(f"b = {b:.4f} \u00B1 {b_err:.4f}")


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(L,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

x_fit=np.linspace(np.min(x),np.max(x),1000)

ax.plot(x_fit,lin(x_fit, *popt),label="Linearer Fit")

ax.axvline(x=np.max(L),color="red",linestyle="--",linewidth=1.5,label=rf"$L_{{\max}} = {np.max(L):.3f}\,\mathrm{{m}}$")

ax.set_xlabel(r"$L$ [m]")
ax.set_ylabel(r"$I$ [mW]")

ax.set_xlim(np.min(L)-0.05,np.max(L)+0.05)
ax.set_ylim(0,5)

ax.grid()
ax.legend()

fig.savefig("../content/images/plankon.pdf")


#### Polarisations Fit ####

pol = np.genfromtxt("pol.txt")

phi = pol[:,0]
I = (pol[:,1]-dunkelWat)
I_err = pol[:,2]

n = len(phi)
print(n)

block = (n + 2) // 3  

with open("tab_pol.txt", "w") as f:
    for i in range(block):

        if i < n:
            left = f"{phi[i]} & {I[i]} & {I_err[i]}"
        else:
            left = " & "

        if i + block < n:
            mid = f"{phi[i+block]} & {I[i+block]} & {I_err[i+block]}"
        else:
            mid = " & "

        if i + 2*block < n:
            right = f"{phi[i+2*block]} & {I[i+2*block]} & {I_err[i+2*block]}"
        else:
            right = " & "

        f.write(f"{left} & {mid} & {right} \\\\\n")


phi = np.deg2rad(pol[:,0])


popt,pcov=curve_fit(cos,phi,I,sigma=I_err,absolute_sigma=True)

A, phi0= popt
a_err, phi0_err = np.sqrt(np.diag(pcov))

print("Polfit")
print(f"A    = {A:.3f}    \u00B1 {A_err:.3f} ")
print(f"phi0 = {phi0:.5f} \u00B1 {phi0_err:.5f} ")
print(f"phi0 = {np.rad2deg(phi0):.5f} \u00B1 {np.rad2deg(phi0_err):.5f} ")
#print(f"C    = {C:.3f}   \u00B1  {C_err:.3f} ")


fig, ax = plt.subplots(layout="constrained")
ax.errorbar(phi,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

phi_fit = np.linspace(phi.min(), phi.max(), 500)

ax.plot(phi_fit,cos(phi_fit, *popt),label=r"$A\cos^2$-Fit")

ax.set_xlabel(r"$\varphi$ [rad]")
ax.set_ylabel(r"$I$ [mW]")

ax.set_xlim(-16,16)
ax.set_ylim(0,8)

ax.grid()
ax.legend()

ax.set_xlim(np.min(phi)-0.02,np.max(phi)+0.02)
ax.set_ylim(0,7)

fig.savefig("../content/images/pol.pdf")



fig,ax=plt.subplots(subplot_kw={"projection": "polar"},layout="constrained")
ax.errorbar(phi,I,yerr=I_err,fmt="x",capsize=3,label="Messwerte")

phi_fit = np.linspace(phi.min(), phi.max(), 500)

ax.plot(phi_fit,cos(phi_fit, *popt),label=r"$A\cos^2$-Fit")

ax.set_xlabel(r"$\varphi$ [rad]")
ax.set_ylabel(r"$I$ [mW]")

ax.set_title("Polarisation")
ax.legend()



fig.savefig("../content/images/pol_polar.pdf")



####Wellenlaenge Gitter 80####

gitter = np.genfromtxt("abstand80.txt")


index = gitter[:,0]
ab = gitter[:,1]

abeff = ab/index

#print(abeff)
lamb=[]
g=1/80
d=72

for i in range(len(abeff)):
        La=(abeff[i]*g)/(np.sqrt(d**2 +ab[i]**2))
        lamb.append(La*1e6)



Lmean=np.mean(lamb)
std = np.std(lamb,ddof=1)

print("Gitter80")
print(lamb)
print(f"lambda    = {Lmean:.8f}    \u00B1 {std:.8f} ")
####Wellenlaenge Gitter 100####



gitter = np.genfromtxt("abstand100.txt")


index = gitter[:,0]
ab = gitter[:,1]

abeff = ab/index

#print(abeff)
lamb=[]
g=1/100
d=72

for i in range(len(abeff)):
        La=(abeff[i]*g)/(np.sqrt(d**2 +ab[i]**2))
        lamb.append(La*1e6)



Lmean=np.mean(lamb)
std = np.std(lamb,ddof=1)

print("Gitter100")
print(lamb)
print(f"lambda    = {Lmean:.8f}    \u00B1 {std:.8f} ")
