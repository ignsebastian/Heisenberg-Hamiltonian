import numpy as np
import cmath
import matplotlib.pyplot as plt
import itertools
from scipy.sparse.linalg import eigsh

#Generating variables
#I use Jx = Jy = Jz = 1 for the coupling constant.
N = 5 #Number of vertices
n = 2**N #Hilbert space dimension
I = np.identity(2)
sx = [[0,0.5],[0.5,0]]
sy = [[0,-0.5j],[0.5j,0]]
sz = [[0.5,0],[0,-0.5]]
s = [sx,sy,sz]
s_big = np.zeros((3,N,n,n),dtype = 'complex_') # np.zeros((Sx (0) Sy (1) Sz (2), 10 vertices, hilbert space dimension, hilbert space dimension))
s_ladder = np.zeros((2,N,n,n),dtype = 'complex_') #np.zeros((plus (0) and minus (1), 10 vertices, hilbert space dimension,hilbert space dimension))

#Calculating the spin operator
for i in range(0,3): #Iterating over x,y,z
  for j in range(0,N): #Iterating over 10 spins
    A = np.copy(I)
    for k in range(0,N-1): #To calculate the tensor product
      if k == 0 and j == 0:
        A = np.kron(s[i],I) # To cover tensor product of spin in the beginning (sx * I * I ...)
      elif k == j-1:
        A = np.kron(A,s[i]) #To cover the tensor product of spin in the middle (I * ... * sx * I ... * I)
      else:
        A = np.kron(A,I) #To cover the identity tensor product
    s_big[i,j] = A #Combining S1 - S10 into one big array

#Calculating the spin ladder operator
for i in range(0,N): #iterating over 10 spins
  s_ladder[0,i] = s_big[0,i] + 1j * s_big[1,i]
  s_ladder[1,i] = s_big[0,i] - 1j * s_big[1,i]

#Preparing Operator O =  (Si)_plus(Si+1)_min - (Si)_min(Si+1)_plus
O = np.zeros((n,n))
for i in range(0,N):
  if i == N-1:
    O = O + (np.dot(s_ladder[0,i],s_ladder[1,0]) + np.dot(s_ladder[1,i],s_ladder[0,0]))
  else:
    O = O + (np.dot(s_ladder[0,i],s_ladder[1,i+1]) + np.dot(s_ladder[1,i],s_ladder[0,i+1]))

#Calculating the hamiltonian
#Note, we can make it more concise by using for i,j in itertools.product(range(0,3), range(0,N)):

H = np.zeros((n,n))
#For periodic boundary, use the following code
#for i in range(0,3):
#  for j in range(0,N):
#    if j == N-1: #For the spin on the edge, where we use periodic boundary condition
#      H = H + np.dot(s_big[i,j],s_big[i,0])
#    else:
#      H = H + np.dot(s_big[i,j],s_big[i,j+1])

#For non periodic boundary
for i in range(0,3):
  for j in range(0,N-1):
    H = H + np.dot(s_big[i,j],s_big[i,j+1])

eval, evec = np.linalg.eigh(H) #Calculating the eigenvalue-eval, and eigenvector, evac\
eval


#Calculating
O_new = np.zeros((n,n), dtype= 'complex_')
for i in range(0,n):
  for j in range(0,n):
    O_new[i,j] = np.dot(np.conj(evec[:,i]),np.dot(O,evec[:,j])) #Calculating the matrix elements of O =  <yn|(Si)^+(Si+1)^- - (Si)^-(Si+1)^+ |ym>

#Plot the histogram
plt.hist(v,bins = 20)
