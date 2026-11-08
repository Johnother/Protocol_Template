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

bmean = bmean(T2)

print(bmean)



