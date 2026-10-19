import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('calibration.csv',delimiter=',',unpack=True)

calibration = {"theta": data[0], "N": data[1]}
print(calibration)
