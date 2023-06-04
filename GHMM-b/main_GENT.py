#coding: utf-8
#(c) Kazuma Nakamura
from __future__ import print_function
import sys,os 
import numpy as np
import argparse
from functions import *
from wrt_GENT import *
from gen_time_series_data import *

#read input from command line 
parser=argparse.ArgumentParser()
parser.add_argument('--N',type=int,default=10000,help='Number of duration-time grids: default=10000')
parser.add_argument('--dlt',type=float,default=0.1,help='Time separation on time-series data [sec]: default=0.1') 
parser.add_argument('--K',type=int,default=20,help='Number of events of ON/OFF in time-series data: default=20') 
parser.add_argument('--noise',type=float,default=0.5,help='Noise amplitude in time-series data: default=0.5') 
parser.add_argument('--cut_off',type=float,default=10,help='Cut off in time-series data: default=10') 
parser.add_argument('--level_ON',type=float,default=101.5,help='ON level of time-series data: default=2.0') 
parser.add_argument('--level_OFF',type=float,default=101.0,help='OFF level of time-series data: default=1.0') 
args = parser.parse_args()
N=args.N  
K=args.K  
dlt=args.dlt 
noise=args.noise 
cut_off=args.cut_off 
level_ON=args.level_ON 
level_OFF=args.level_OFF 

#make Pon and Qon 
if False: 
#if True:
    A_ON=1.0 #Amplitude of P_ON(td) 
    R_ON=20.0 #Relaxation time of P_ON(td) [sec]
    Pon=P_ON_exponential(A_ON,R_ON,dlt,N) 
    A_OFF=1.0 #Amplitude of P_OFF(td) 
    R_OFF=40.0 #Relaxation time of P_OFF(td) [sec] 
    Poff=P_OFF_exponential(A_OFF,R_OFF,dlt,N) 
    print("EXP")
else: 
    A_ON=1.0 #Amplitude of P_ON(td) 
    R_ON=1800.0#1800.0 #7.0 #Relaxation time of P_ON(td) [sec] 
    mul_ON=1.30#0.50 #Multiplier of P_ON(td)
    Pon=P_ON_expo_power(A_ON,R_ON,mul_ON,dlt,N,cut_off)
    #
    A_OFF=1.0 #Amplitude of P_OFF(td) 
    R_OFF=1800.0#1400.0 #7.0 #Relaxation time of P_ON(td) [sec] 
    mul_OFF=1.7#0.50#1.4#1.5#1.6#1.70 #1.40 #Multiplier of P_ON(td)
    Poff=P_OFF_expo_power(A_OFF,R_OFF,mul_OFF,dlt,N,cut_off)
    print("POWER")
Pon,Qon=normalize(dlt,Pon) 
Poff,Qoff=normalize(dlt,Poff) 

#write 
wrt_Pon_Poff(dlt,Pon,Poff)
wrt_Qon_Qoff(dlt,Qon,Qoff)

#check 
XQ=Qon 
#XQ=Qoff 
Xdata=check(dlt,XQ) 
wrt_check_data(dlt,Xdata)
#sys.exit() 

#make time-series data 
(time,L_timeON,L_timeOFF)=gen_time_series(Qon,Qoff,K,level_ON,level_OFF) 
#write 
wrt_time_series_no_noise(dlt,time) 

#superpose noise on time-series data 
time=add_noise(noise,time) 

#write 
wrt_time_series(dlt,time) 
wrt_T_ON(K,L_timeON,L_timeOFF) 
wrt_T_OFF(K,L_timeON,L_timeOFF) 
