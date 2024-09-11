import numpy as np
import matplotlib.pyplot as plt

X = np.arange(1.920,2.080,0.001)
Y = X**9 -18*X**8 + 144*X**7-672*X**6+2016*X**5-4032*X**4+5376*X**3-4608*X**2+2304*X-512
plt.plot(X,Y)
plt.xlabel('x')
plt.ylabel('p(x)')
plt.show()

Y1 = (X-2)**9

plt.plot(X,Y1)
plt.xlabel('x')
plt.ylabel('p(x)')
plt.show()

Ydiff = Y1 - Y
plt.plot(X,Ydiff)
plt.xlabel('x')
plt.ylabel('difference')
plt.show()

X = 2.05
Y2 = X**9 -18*X**8 + 144*X**7-672*X**6+2016*X**5-4032*X**4+5376*X**3-4608*X**2+2304*X-512
Y3 = (X-2)**9
print(Y2)
print(Y3)