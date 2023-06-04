# coding: utf-8
import numpy as np
from numpy.random import *

def P_ON_exponential(A,alp,dlt,N):
    x=np.zeros(N)
    for n in range(N):
        x[n]=A*np.exp(-n/alp*dlt)
    return x

def P_ON_expo_power(A,R_ON,mul_ON,dlt,N,cut_off):
    x=np.zeros(N)
    eps=0.0000001
    for n in range(N):
        if(n<cut_off):#100):
        #if(n<10):#100):
        #if(n<50):#100):
            x[n]=0.0
        else: 
            #x[n]=A*(((float(n)+eps)*dlt)**(-mul_ON))*np.exp(-float(n)/R_ON*dlt)
	    #x[n]=A*((float(n+1)*dlt)**(-mul_ON))*np.exp(-float(n+1)/R_ON*dlt)
            x[n]=A*((float(n)*dlt)**(-mul_ON))*np.exp(-float(n)/R_ON*dlt)
	    #x[n]=A*np.exp(-float(n)/R_ON*dlt)
    return x

def P_OFF_exponential(B,bet,dlt,N):
    y=np.zeros(N)
    for n in range(N): 
        y[n]=B*np.exp(-n/bet*dlt) 
    return y 

def P_OFF_expo_power(B,R_OFF,mul_OFF,dlt,N,cut_off):
    y=np.zeros(N)
    for n in range(N):
        if(n<cut_off):#100):
        #if(n<10):#100):
        #if(n<50):#100):
            y[n]=0.0
        else:  
            #y[n]=B*((float(n+1)*dlt)**(-mul_OFF))
            y[n]=B*((float(n)*dlt)**(-mul_OFF))
            #y[n]=B*((float(n+1)*dlt)**(-mul_OFF))*np.exp(-float(n+1)/R_OFF*dlt)
    return y

def normalize(dlt,x): 
    N=len(x)  
    s=0.0
    for n in range(N):
        s=s+x[n]*dlt 
    #print(s) 
    #x=x*dlt/s 
    #x2=x/dlt
    #
    x=x/s
    x2=x 
    #
    q=np.zeros(N)
    qd=np.zeros(N)
    s=0.0
    for n in range(N):
        s=s+x2[n]*dlt 
        q[n]=s
    return x,q

def check(dlt,QX):
    N=len(QX)
    Nevent=1*N #number of events
    hist=np.zeros(Nevent)
    for event in range(Nevent):
         r=np.random.random()
         for n in range(N):
             if r <= QX[n]:
                  hist[event]=float(n)*dlt
                  break
    #
    Xdata=np.zeros(N)
    for event in range(Nevent):
        event_n=int(round(hist[event]/dlt))
        Xdata[event_n]=Xdata[event_n]+1.0
    #
    Xdata,QXdata=normalize(dlt,Xdata)
    return Xdata
