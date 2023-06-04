#coding: utf-8
#(c) Tatsuhiro Furuta
import sys,os
import numpy as np
import argparse

def input_file(file_name):
    f=open(file_name,'r')
    line=f.readlines()
    f.close()

    tmp=np.loadtxt(line)
    #print 'tmp.shape=',tmp.shape

    N=tmp.shape[0]
    Nd=tmp.shape[1]
    #print('N,Nd=',N,Nd)

    x=np.zeros((Nd,N))
    for i in range(Nd):
        for j in range(N):
            x[i][j]=float(tmp[j][i])

    z_mat=np.zeros((1,N)) 
    t_mat=np.zeros((1,N))
    
    t_mat=x[0,:]
    t_mat=t_mat.reshape(1,N)
    z_mat=x[1,:] 
    z_mat=z_mat.reshape(1,N)
    #print(z_mat.shape)
    #print(z_mat[0][0])

    return N,z_mat
