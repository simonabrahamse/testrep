import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
from scipy.linalg import lu_solve
import time

def driver():

     ''' create  matrix for testing different ways of solving a square 
     linear system'''

     '''' N = size of system'''
     N = 500
 
     ''' Right hand side'''
     b = np.random.rand(N,1)
     A = np.random.rand(N,N)
     time1 = time.time()
     for i in range(100):
          x = scila.solve(A,b)
     
     time2 = time.time()
     time3 = (time2-time1)/100
     print(time3)
     test = np.matmul(A,x)
     
     r = la.norm(test-b)
     
     time1 = time.time()
     for i in range(100):
          lu, piv = scila.lu_factor(A)
     time2 = time.time()
     time3 = (time2-time1)/100
     print('Time taken for LU factorization:',time3)

     time1 = time.time()
     for i in range(100):
          x = lu_solve((lu,piv),b)
     time2 = time.time()
     time4 = (time2-time1)/100
     print('Time taken to do the LU solve:',time4)
     print('Total time taken for LU:',time3+time4)

     # print(r)

     ''' Create an ill-conditioned rectangular matrix '''
     N = 10
     M = 5
     A = create_rect(N,M)     
     b = np.random.rand(N,1)




     
def create_rect(N,M):
     ''' this subroutine creates an ill-conditioned rectangular matrix'''
     a = np.linspace(1,10,M)
     d = 10**(-a)
     
     D2 = np.zeros((N,M))
     for j in range(0,M):
        D2[j,j] = d[j]
     
     '''' create matrices needed to manufacture the low rank matrix'''
     A = np.random.rand(N,N)
     Q1, R = la.qr(A)
     test = np.matmul(Q1,R)
     A =    np.random.rand(M,M)
     Q2,R = la.qr(A)
     test = np.matmul(Q2,R)
     
     B = np.matmul(Q1,D2)
     B = np.matmul(B,Q2)
     return B     
          
  
if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()       
