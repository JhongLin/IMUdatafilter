import numpy as np
import math
import matplotlib.pylab as plt



x = np.linspace(0, 2*np.pi, 1000)
'''
sinwave = np.sin(x)
plt.plot(x, sinwave)
plt.xlabel('Angle [rad]')
plt.ylabel('sin(x)')
plt.axis('tight')
plt.show()
sin_integral = np.zeros(len(sinwave))
for i in range(len(sinwave)):
    sin_integral[i] = np.sum(sinwave[:i+1])

plt.plot(x, sin_integral)
plt.xlabel('Angle [rad]')
plt.ylabel('sin_integral')
plt.show()
'''


coswave = np.cos(x)
plt.plot(x, coswave)
plt.xlabel('Angle [rad]')
plt.ylabel('cos(x)')
#plt.axis('tight')
plt.show()
cos_integral = np.zeros(len(coswave))
for i in range(len(coswave)):
    cos_integral[i] = np.sum(coswave[:i+1])

plt.plot(x, cos_integral)
plt.xlabel('Angle [rad]')
plt.ylabel('cos_integral')
plt.show()
