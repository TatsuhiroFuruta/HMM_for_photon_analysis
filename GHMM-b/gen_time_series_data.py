# coding: utf-8
import numpy as np
from numpy.random import *

def gen_time_series(Qon,Qoff,K,level_ON,level_OFF): 
    L_timeON=np.zeros(K,dtype=int)
    L_timeOFF=np.zeros(K,dtype=int)
    N=len(Qon) 
    for event in range(K): 
        r=np.random.rand() 
        for n in range(N): 
            if(r<Qon[n]): 
               break
        r=np.random.rand() 
        for m in range(N): 
            if(r<Qoff[m]): 
               break
        L_timeON[event]=n+1
        L_timeOFF[event]=m+1
    for event in range(K):
        print('event,L_timeON',event,L_timeON[event])
        print('event,L_timeOFF',event,L_timeOFF[event])
    #Ntot=np.sum(L_timeON)+np.sum(L_timeOFF)-K 
    Ntot=np.sum(L_timeON)+np.sum(L_timeOFF) 
    print(Ntot)
    time=np.zeros(Ntot)
    i=0
    for event in range(K): 
        L=L_timeON[event] 
        M=L_timeOFF[event] 
        #time[i:i+L]=2.0 
        time[i:i+L]=level_ON 
        #time[i+L:i+L+M]=1.0 
        time[i+L:i+L+M]=level_OFF 
        i=i+L+M 
        print(L,M,i) 
    return time,L_timeON,L_timeOFF  

def add_noise(noise,time): 
    L=len(time) 
    for i in range(L): 
        #r=noise*(np.random.rand()-0.5) 
        #r=np.random.normal(time[i],noise,1) 
        r=np.random.normal(0,noise,1) 
        time[i]=time[i]+r 
    return time  
