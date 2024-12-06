import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-np.pi, np.pi, 0.001)

def f(x):
    return np.cosh(x)

def dfdx(f, x):
    DF = [(f[i+1] - f[i])/(x[i+1]-x[i]) for i in range(len(x)-1)]
    df = np.append(DF,DF[-1])
    return df


plt.plot(x,dfdx(f(x), x))
plt.plot(x,f(x))
plt.show()


