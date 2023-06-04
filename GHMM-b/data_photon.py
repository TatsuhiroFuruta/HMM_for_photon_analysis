# coding:utf-8
#(c) Kazuma Nakamura,Tatsuhiro Furuta
import sys, os
import numpy as np
import scipy 

def mktrain_x_mat(input_file,s):
    f=open(input_file,'r') 
    #f=open('photon.data01.txt','r') 
    #f=open('photon.data02.txt','r') 
    #f=open('wet-009.txt','r') 
    #f=open('vac-001.txt','r') 
    #f=open('T.csv','r') 
    #f=open('T-50.csv','r') 
    line=f.readlines()#[1:] 
    f.close()
    #tmp=np.genfromtxt(line,delimiter=",",dtype="str")
    tmp=np.loadtxt(line)
    #print tmp.shape 
    #sys.exit() 
    #
    Nd=tmp.shape[0] 
    N=tmp.shape[1] 
    #print Nd,N 
    #sys.exit() 
    #
    x=np.zeros((N,Nd))
    for n in range(Nd): 
        for d in range(N): 
            x[d][n]=float(tmp[n][d]) 
    #print x
    #print x.shape 
    #
    D=N-1
    x_time_mat=np.zeros((D,Nd))
    x_time_mat=x_time_mat.reshape(D,Nd)
    for n in range(Nd):
        x_time_mat[D-1][n]=x[0][n] 
    #print D
    x_mat=np.zeros((D,Nd))
    x_mat=x_mat.reshape(D,Nd)
    for n in range(Nd):
        x_mat[D-1][n]=x[1][n] 
    #a=120
    #x_time_mat=np.zeros((D,Nd-a))
    #x_time_mat=x_time_mat.reshape(D,Nd-a)
    #for n in range(Nd-a):
        #x_time_mat[D-1][n]=x[0][n+a] 
    #print D
    #x_mat=np.zeros((D,Nd-a))
    #x_mat=x_mat.reshape(D,Nd-a)
    #for n in range(Nd-a):
        #x_mat[D-1][n]=x[1][n+a] 
    #for n in range(Nd): 
        #print x[0][n],x[1][n] 
    #Nd=Nd-a

    x_mean=0
    for n in range(Nd):
        x_mean=x_mean+x_mat[D-1][n]
    x_mean=x_mean/Nd
     
    x=0
    sigma=0
    for n in range(Nd):
        x=x_mat[D-1][n]-x_mean
        sigma=sigma+x**2
    sigma=sigma/Nd
    sigma=np.sqrt(sigma)

    #s=0
    N_renormalize=0
    for n in range(Nd):
        if np.abs(x_mat[D-1][n]-x_mean)>s*sigma:
            continue
        else:
            N_renormalize=N_renormalize+1
     
    x_renormalize_time_mat=np.zeros((D,N_renormalize))
    x_renormalize_time_mat=x_renormalize_time_mat.reshape(D,N_renormalize)     
    x_renormalize_mat=np.zeros((D,N_renormalize))     
    x_renormalize_mat=x_renormalize_mat.reshape(D,N_renormalize)
    n_renormalize=0
    for n in range(N_renormalize):
        x_renormalize_time_mat[D-1][n]=x_time_mat[D-1][n]

    for n in range(Nd):
        if np.abs(x_mat[D-1][n]-x_mean)>s*sigma:
            n_renormalize=n_renormalize+1
        else:
            #x_renormalize_time_mat[D-1][n-n_renormalize]=x_time_mat[D-1][n]
            x_renormalize_mat[D-1][n-n_renormalize]=x_mat[D-1][n]
         
    return Nd,N_renormalize,D,x_time_mat,x_mat,x_renormalize_time_mat,x_renormalize_mat,x_mean,sigma

def x_integer(N,x_time_mat,x_mat,x_int_term):
    x_min=np.min(x_mat)
    x_int_mat=np.zeros((1,N))
    for n in range(N):
        x_int_mat[0][n]=(x_mat[0][n]-x_min)*x_int_term+1.0
        x_int_mat[0][n]=int(x_int_mat[0][n])
    x_mat=x_int_mat
    f=open('x_int_mat.csv','w')
    for i in range(N): 
        f.write(str(x_time_mat[0,i]))
        f.write('\t') 
        f.write(str(x_mat[0,i]))      
        f.write('\n')
    f.close()
    return x_mat

def x_scaling(N,normalized,x_scale,x_shift,x_time_mat,x_mat):
    x_max=np.max(x_mat)
    x_min=np.min(x_mat)
    if normalized==True:
        x_mat=(x_mat-x_min*np.ones((1,N)))/(x_max-x_min)
    if x_scale!=1.0:
        x_mat=x_scale*x_mat
    if x_shift!=0.0:
        x_mat=x_mat+x_shift
    f=open('x_scaling_mat.csv','w')
    for i in range(N): 
        f.write(str(x_time_mat[0,i]))
        f.write('\t') 
        f.write(str(x_mat[0,i]))      
        f.write('\n')
    f.close()
    return x_mat

def x_norm(N,x_time_mat,x_mat):
    x_mat=x_mat.flatten()
    x_bar=np.mean(x_mat)
    x_mean=x_bar*np.ones(N)
    x_sd=np.std(x_mat)
    x_mat=x_mat-x_mean
    x_mat=x_mat/x_sd
    x_mat=x_mat.reshape(1,N)
    f=open('x_norm_mat.csv','w')
    for i in range(N): 
        f.write(str(x_time_mat[0,i]))
        f.write('\t') 
        f.write(str(x_mat[0,i]))      
        f.write('\n')
    f.close()
    return x_mat,x_mean,x_sd
    
   
    
