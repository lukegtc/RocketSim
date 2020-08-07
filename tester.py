import matplotlib.pyplot as plt
from math import *
import numpy as np
g = 0.5
beta = np.linspace(0,np.pi,100)
rho1 = np.sqrt((g*np.sin(beta))**2 + (g**2*(np.cos(beta)-1)+1)**2)
g = 0.9
rho2 = np.sqrt((g*np.sin(beta))**2 + (g**2*(np.cos(beta)-1)+1)**2)
g = 1.1
rho3 = np.sqrt((g*np.sin(beta))**2 + (g**2*(np.cos(beta)-1)+1)**2)

plt.plot(beta,rho1)
plt.plot(beta,rho2)
plt.plot(beta,rho3)
plt.show()