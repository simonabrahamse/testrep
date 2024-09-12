import numpy as np
import matplotlib.pyplot as plt

t=np.linspace(0,np.pi,31)
# print(np.pi/30)
# print(29*np.pi/30)
y= np.cos(t)
S = 0
for i in range(0,30):
    S = S + t[i]*y[i]
    
print('The sum is',S)

Theta = np.linspace(0,2*np.pi,1000)
R = 1.2
dr = 0.1
p= 0
f = 15

xtheta = R*(1+dr*np.sin(f*Theta+p))*np.cos(Theta)
ytheta = R*(1+dr*np.sin(f*Theta+p))*np.sin(Theta)

fig1 = plt.figure(figsize =(10, 10))
plt.plot(xtheta,ytheta)
plt.xlabel('x(theta)')
plt.ylabel('y(theta)')
plt.show()

fig2 = plt.figure(figsize=(10,8))
for i in range(1,11):
    R = i
    dr = 0.05
    f = 2+ i
    p = np.random.uniform(0.0,2.0)
    xtheta = R*(1+dr*np.sin(f*Theta+p))*np.cos(Theta)
    ytheta = R*(1+dr*np.sin(f*Theta+p))*np.sin(Theta)
    
    sub = plt.subplot(2,5,i)
    sub.plot(xtheta,ytheta)
    # plt.title('i = %i') % i
    plt.title('i = {}'.format(i))
plt.show()