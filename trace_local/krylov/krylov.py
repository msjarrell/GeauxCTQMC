from numpy import *
from numpy.random import rand
from numpy.linalg import norm
from scipy.linalg import expm

N=1024
M=3

def init_matrix(n):
    A=matrix(zeros([n,n]))
    
    for i in range(0,n-1):
        A[i,i+1]=-1
        A[i+1,i]=-1

    A[0,n-1]=-1
    A[n-1,0]=-1

    for i in range(n):
        A[i,i]=2.0*rand()-1.0

    return A

def Arnoldi(vec,A,m):
    n=len(vec)
    V=zeros([n,m])
    H=zeros([m,m])
    vec=vec/norm(vec)
    V[:,0]=vec
    for j in range(m-1):
        w=array(A.dot(vec))[0]
        for i in range(j+1):
            h=dot(w,V[:,i])
            H[i,j]=h
            w=w-h*V[:,i]

        vec=w/norm(w)
        V[:,j+1]=vec
        H[j+1,j]=norm(w)

    j=m-1
    w=array(A.dot(vec))[0]
    for i in range(j+1):
        h=dot(w,V[:,i])
        H[i,j]=h
        w=w-h*V[:,i]
        
    return [matrix(V),matrix(H)]

def krylov(vec,A,t):
    V,H=Arnoldi(vec,A,M)
    return norm(vec)*V.dot(expm(H*t)[0])

A=init_matrix(N)
vec=rand(N)
tau=0.1
error=sum((array(krylov(vec,A,tau))[0] - dot(expm(A*tau),vec))**2)
print error

