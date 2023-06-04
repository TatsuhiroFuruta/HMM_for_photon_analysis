#coding: utf-8
#(c) Kazuma Nakamura,Tatsuhiro Furuta
import sys,os
import numpy as np

def distribution_density(a,A):
    #a=0.01
    d=0
    #l=np.arange(np.min(A),np.max(A)+a,a)
    l=np.arange(0.0,np.max(A)+2.0*a,a)
    l=l.reshape(l.shape[0],1)
    #print'l.shape=',l.shape
    #print 

    D_mat=np.zeros((2,l.shape[0]))
    D_mat=D_mat.reshape(2,l.shape[0])
    #print'D_mat.shape=',D_mat.shape
    #print'min(A)',np.min(A)
    #for i in range(l.shape[0]):
       #print 'l[i][0]=',l[i][0]
    D_mat[0,:]=l[:,0]

    for j in range(l.shape[0]-1):
        for n in range(A.shape[1]):
            #if D_mat[0,j]-0.5*a<=A[0,n]<D_mat[0,j]+0.5*a:
            if D_mat[0,j]<A[0,n]<=D_mat[0,j+1]:
                D_mat[1,j+1]=D_mat[1,j+1]+1
    for i in range(l.shape[0]):
        d=d+D_mat[1,i]*a
           #print'D_mat[1,j]=',D_mat[1,j]
    #D_mat[1,:]=D_mat[1,:]/A.shape[1]
    D_mat[1,:]=D_mat[1,:]/d

    #print'd/A=',d/A.shape[1]

    return D_mat

def distribution_density_for_few_photon_data(a,A):
    #a=0.01
    d=0
    #l=np.arange(np.min(A),np.max(A)+a,a)
    l=np.arange(0.0,np.max(A)+2.0*a,a)
    N_tau=l.shape[0]
    l=l.reshape(N_tau,1)
    #print'l.shape=',l.shape
    #print 

    D_mat=np.zeros((2,N_tau))
    #print'D_mat.shape=',D_mat.shape
    #print'min(A)',np.min(A)
    #for i in range(l.shape[0]):
       #print 'l[i][0]=',l[i][0]
    D_mat[0,:]=l[:,0]
    N_T=A.shape[1]
    for j in range(N_tau-1):
        for n in range(N_T):
            #if D_mat[0,j]-0.5*a<=A[0,n]<D_mat[0,j]+0.5*a:
            if D_mat[0,j]<A[0,n]<=D_mat[0,j+1]:
                D_mat[1,j+1]=D_mat[1,j+1]+1

    k=1
    while int(D_mat[1][k-1])==0:
        k=k+1
    m=1
    while int(D_mat[1][N_tau-m])==0:
        m=m+1
    P_mat=np.zeros((2,N_tau-k-m))
    P_mat[0,:]=l[k:N_tau-m,0]
    
    for i in range(k,N_tau-m):
        j=1
        while int(D_mat[1][i-j])==0:
            if i-j==k-1:
                break
            j=j+1
        f=i-j
     
        j=1
        while int(D_mat[1,i+j])==0:
            if i+j==N_tau-m:
                break
            j=j+1
        b=i+j
        tau_f=D_mat[0,i]-D_mat[0,f]
        tau_b=D_mat[0,b]-D_mat[0,i]
        tau=(tau_f+tau_b)/2.0
        P_mat[1,i-k]=float(D_mat[1,i])/tau
 
    P_hat_mat=P_mat[1,:]*a
    d=P_hat_mat.sum()
    P_mat[1,:]=P_mat[1,:]/d

    return P_mat

def distribution_density_weighted(a,A):
    T_s_mat=np.sort(A)
    T_s_mat=T_s_mat.flatten()
    #print(T_s_mat.shape)
    #sys.exit()
    T_g=np.arange(0.0,np.max(A)+2.0*a,a)
    #print(T_g.shape)
    #sys.exit()

    N=T_s_mat.shape[0]
    M=T_g.shape[0]
    print('# N=',N)
    print('# M=',M)
    
    T_tilde_mat=np.zeros(N)
    f_tilde_mat=np.zeros(N)
    N_irr=0

    for i in range(N):
        for m in range(M-1):
            if T_g[m]<T_s_mat[i]<=T_g[m+1]:
                if i==0:
                    N_irr+=1
                    T_tilde_mat[i]=T_g[m+1]
                else:
                    if T_g[m]<T_s_mat[i-1]<=T_g[m+1]:
                        f_tilde_mat[i]+=1
                    else:
                        N_irr+=1
                        T_tilde_mat[i]=T_g[m+1]
              
    T_irr_mat=np.zeros(N_irr)
    f_irr_mat=np.ones(N_irr)
    
    l=0
    for i in range(N):
        k=1
        if f_tilde_mat[i]==0:
            T_irr_mat[i-l]=T_tilde_mat[i]
            if i!=N-1:
                while f_tilde_mat[i+k]==1:
                    f_irr_mat[i-l]+=1
                    k+=1
                    if i+k==N:
                        break
            l+=k-1

    D_mat=np.zeros(N_irr-2)
    
    for i in range(1,N_irr-1):
        T=T_irr_mat[i+1]-T_irr_mat[i-1]
        w=2.0/T
        D_mat[i-1]=w*float(f_irr_mat[i])

    D_hat_mat=D_mat*a
    d=D_hat_mat.sum()
    print('# d=',d)
    D_mat=D_mat/d
    
    P_mat=np.zeros((2,N_irr-2))
    P_mat[0,:]=T_irr_mat[1:N_irr-1]
    P_mat[1,:]=D_mat

    return P_mat
    
    
    
