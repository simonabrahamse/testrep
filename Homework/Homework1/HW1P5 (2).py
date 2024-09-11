import numpy as np
import matplotlib.pyplot as plt
import pylab

X1 = np.pi
X2 = 10**6

exp = np.linspace(-16,0,17)
d = 10**exp

F1a = np.cos(X1+d) - np.cos(X1)
F1b =-2*np.sin((2*X1+d)/2)*np.sin(d/2)
F2a = np.cos(X2+d) - np.cos(X2)
F2b =  -2*np.sin((2*X2+d)/2)*np.sin(d/2)

F3a = d * -1* np.sin(X1) - ((d**2)/2)*np.cos(X1+d)
F3b = d * -1* np.sin(X2) - ((d**2)/2)*np.cos(X2+d)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

line, = ax.plot(d,F1a,label="Original")
line, = ax.plot(d,F1b,label="Sine Expansion")
line, = ax.plot(d,F3a,label="Taylor expansion",linestyle='dashed')
plt.legend(loc="upper left")
plt.title("Comparison for x=pi")
plt.xlabel("Delta")
plt.ylabel("Difference")

ax.set_xscale('log')
ax.set_yscale('log')

pylab.show()


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

line, = ax.plot(d,F2a,label="Original")
line, = ax.plot(d,F2b,label="Sine Expansion")
line, = ax.plot(d,F3b,label="Taylor series appproximation",linestyle='dashed')
plt.legend(loc="upper left")
plt.title("Comparison for x=10^6")
plt.xlabel("Delta")
plt.ylabel("Difference")

ax.set_xscale('log')
ax.set_yscale('log')

pylab.show()


