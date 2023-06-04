# coding: utf-8
import numpy as np
from numpy.random import *

def wrt_Pon_Poff(dlt,Pon,Poff): 
    f=open('P.csv','w')
    N=len(Pon)  
    for n in range(N): 
        f.write(str(dlt*float(n)))
        f.write('\t')
        f.write(str(Pon[n]))
        f.write('\t')
        f.write(str(Poff[n]))
        f.write('\n')
    f.close()

def wrt_Qon_Qoff(dlt,Qon,Qoff): 
    f=open('Q.csv','w')
    N=len(Qon)  
    for n in range(N): 
        f.write(str(dlt*float(n)))
        f.write('\t')
        f.write(str(Qon[n]))
        f.write('\t')
        f.write(str(Qoff[n]))
        f.write('\n')
    f.close()

def wrt_QDon_QDoff(dlt,QDon,QDoff): 
    f=open('QD.csv','w')
    N=len(QDon)  
    for n in range(N): 
        f.write(str(dlt*float(n)))
        f.write('\t')
        f.write(str(QDon[n]))
        f.write('\t')
        f.write(str(QDoff[n]))
        f.write('\n')
    f.close()

def wrt_check_data(dlt,Xdata): 
    f=open('check.csv','w')
    N=len(Xdata)  
    for event_n in range(N): 
        f.write(str(dlt*float(event_n)))
        f.write('\t')
        f.write(str(Xdata[event_n]))
        f.write('\n')
    f.close()
    return 

def wrt_time_series_no_noise(dlt,time): 
    f=open('T_no_noise.csv','w')
    Ntot=len(time) 
    for n in range(Ntot): 
        if n % 1 ==0: 
            f.write(str(dlt*float(n)))
            f.write('\t')
            f.write(str(time[n]))
            f.write('\n')
    f.close()

def wrt_time_series(dlt,time): 
    f=open('T.csv','w')
    Ntot=len(time) 
    for n in range(Ntot): 
        if n % 1 ==0: 
            f.write(str(dlt*float(n)))
            f.write('\t')
            f.write(str(time[n]))
            f.write('\n')
    f.close()

def wrt_T_ON(K,L_timeON,L_timeOFF): 
    f=open('T_ON.txt','w')
    M=0
    for event in range(K):
        f.write(str(M))
        f.write('\t')
        f.write(str(M+L_timeON[event]))
        f.write('\n')
        M=M+L_timeON[event]+L_timeOFF[event]
    f.close()

def wrt_T_OFF(K,L_timeON,L_timeOFF): 
    f=open('T_OFF.txt','w')
    M=0
    for event in range(K):
        f.write(str(M+L_timeON[event]))
        f.write('\t')
        M=M+L_timeON[event]+L_timeOFF[event]
        f.write(str(M))
        f.write('\n')
    f.close()
