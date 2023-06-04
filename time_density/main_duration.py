#coding: utf-8
#(c) Tatsuhiro Furuta
import sys,os
import numpy as np
import argparse
import scipy.interpolate as interp
from distutils.util import strtobool
from duration_time import *
from answer_check import *
from wrt_duration_time import *
import time 

t1=time.time() 
parser=argparse.ArgumentParser()
parser.add_argument('--normal_mode',type=strtobool,default=False,help='The parameter for duration time function.(default: False)') 
parser.add_argument('--M',type=int,default=0,help='file numbers of s_mat-ON.csv(default: 0)') 
parser.add_argument('--x',type=float,default=0.5,help='The boundary between ON and OFF(default: 0.5)') 
parser.add_argument('--T_x',type=float,default=0.5,help='The boundary between ON and OFF for answer check(default: 0.5)') 
parser.add_argument('--ON_line',type=float,default=0.9,help='The boundary of ON in duration_time (default: 0.9)') 
parser.add_argument('--OFF_line',type=float,default=0.1,help='The boundary of OFF in duration_time (default: 0.1)') 
parser.add_argument('--event_error',type=str,default='OFF',help='If event_error == ON,the event number error program turns on.(default: OFF)') 
parser.add_argument('--event_number',type=int,default=2,help='The event number of time_seris_data (default: 2)') 
parser.add_argument('--out_ON_OFF',type=str,default='OFF',help='The parameter for outputs of duration time files.(default: OFF)') 
parser.add_argument('--out_duration_mat',type=str,default='OFF',help='The parameter for output of duration time.csv.(default: OFF)') 
parser.add_argument('--answer_check',type=str,default='OFF',help='If answer_check == ON,the event check program turns on.(default: OFF)') 
parser.add_argument('--out_ans',type=str,default='OFF',help='The parameter for outputs of T_ON_mat and T_OFF_mat.(default: OFF)') 
parser.add_argument('--T_ON_line',type=float,default=0.9,help='The boundary of T_ON in duration_time (default: 0.9)') 
parser.add_argument('--T_OFF_line',type=float,default=0.1,help='The boundary of T_OFF in duration_time (default: 0.1)') 
parser.add_argument('--t_grid',type=float,default=1.0,help='The parameter which determines the minimum value of the range for duration time on the answer check. (default: 1.0)') 
parser.add_argument('--N_eve',type=int,default=3,help='The parameter where t_grid*N_eve determines the maximum value of range for duration time on the answer check. (default: 3)') 
parser.add_argument('--T_sig',type=float,default=0.3,help='The boundary of ON_mat-T_ON_mat or OFF_mat-T_OFF_mat in the answer check program, where the parameter represents the allowable error per 1 sec. (default: 0.3)') 
parser.add_argument('--sim_eva',type=strtobool,default=True,help='The parameter for num fuction. (default: True)') 
parser.add_argument('--tau_max',type=float,default=1000.0,help='The parameter for num function. (default: 1000.0)') 
parser.add_argument('--tau_min',type=float,default=1.0,help='The parameter for num function. (default: 1.0)') 

args = parser.parse_args()
normal_mode=args.normal_mode
M=args.M
x=args.x
T_x=args.T_x
ON_line=args.ON_line
OFF_line=args.OFF_line
event_error=args.event_error
event_number=args.event_number
out_ON_OFF=args.out_ON_OFF
out_duration_mat=args.out_duration_mat
answer_check=args.answer_check
out_ans=args.out_ans
T_ON_line=args.T_ON_line
T_OFF_line=args.T_OFF_line
t_grid=args.t_grid
N_eve=args.N_eve
T_sig=args.T_sig
sim_eva=args.sim_eva
tau_max=args.tau_max
tau_min=args.tau_min
#x=0.5

#ON_line=0.6
#OFF_line=0.4

print('# normal_mode=',bool(normal_mode))
print('# M=',M)
print('# x=',x)
print('# T_x=',T_x)
print('# ON_line=',ON_line)
print('# OFF_line=',OFF_line)
print('# event_error=',event_error)
print('# event_number=',event_number)
print('# out_ON_OFF=',out_ON_OFF)
print('# out_duration_mat=',out_duration_mat)
print('# answer_check=',answer_check)
print('# out_ans=',out_ans)
print('# T_ON_line=',T_ON_line)
print('# T_OFF_line=',T_OFF_line)
print('# t_grid=',t_grid)
print('# N_eve=',N_eve)
print('# T_sig=',T_sig)
print('# sim_eva=',bool(sim_eva))
print('# tau_max=',tau_max)
print('# tau_min=',tau_min)

#T_sig_tilde=T_sig*t_grid
#print 'T_sig_tilde=',T_sig_tilde

#num_mat=num(t_grid,N_eve)
num_mat=num(sim_eva,t_grid,N_eve,tau_max,tau_min)
#T_sig_mat=Generate_T_sig_mat(T_sig_tilde,N_eve)

W_pre_sum_tot_mat=np.zeros((1,N_eve+2))
W_sum_tot_mat=np.zeros((1,N_eve+2))
W_sum_ans_mat=np.zeros((1,N_eve+2))
W_sum_suc_mat=np.zeros((1,N_eve+2))
V_sum_tot_mat=np.zeros((1,N_eve+2))
V_sum_ans_mat=np.zeros((1,N_eve+2))
V_sum_suc_mat=np.zeros((1,N_eve+2))

print('# num_mat=',num_mat)
#print "T_sig_mat=",T_sig_mat

if out_ON_OFF=='OFF' and answer_check=='OFF':
    print('The values of parameter are out_ON_OFF==OFF and answer_check==OFF. If you carry out the program, you have to turn on either or both of out_ON_OFF and answer_check.')
else:
    for m in range(M):
        f=open('s_mat-ON-'+str(m)+'.csv','r')
        line=f.readlines()
        f.close()

        tmp=np.loadtxt(line)
        #print 'tmp.shape=',tmp.shape

        N=tmp.shape[0]
        Nd=tmp.shape[1]

        z=np.zeros((Nd,N))
        z=z.reshape(Nd,N)
        for i in range(Nd):
            for j in range(N):
                z[i][j]=float(tmp[j][i])

        z_mat=np.zeros(N) 
        z_mat=z_mat.reshape(1,N) 
        t_mat=np.zeros(N)
        t_mat=t_mat.reshape(1,N)

        for n in range(N):
            z_mat[0][n]=z[1][n]
            t_mat[0][n]=z[0][n]

        '''
        Ns=int(N/20)
        tt=np.linspace(min(t_mat),max(t_mat),Ns)
        zz=interp.spline(t_mat.z_mat,tt)
        '''
        if normal_mode==True:
            (ON_mat,OFF_mat,duration_mat,ON_ini_mat,ON_end_mat,OFF_ini_mat,OFF_end_mat)=duration_time(N,x,z_mat,t_mat)
        else:
            (ON_mat,OFF_mat,duration_mat,ON_ini_mat,ON_end_mat,OFF_ini_mat,OFF_end_mat)=duration_time_kai(N,ON_line,OFF_line,z_mat,t_mat)
        #for i in range(ON_mat.shape[1]):
        #   print ON_ini_mat[i] 
        #   print ON_end_mat[i] 
        #   print OFF_ini_mat[i] 
        #   print OFF_end_mat[i] 
        #for i in range(OFF_mat.shape[1]): 
        #(ON_mat,OFF_mat,duration_mat)=duration_time(N,x,zz,tt)
        #print "OFF_mat.shape[1]=",OFF_mat.shape[1]
        #print "OFF_mat.shape=",OFF_mat.shape
        if answer_check=='ON':
            f=open('T_no_noise-'+str(m)+'.csv','r')
            line=f.readlines()
            f.close()

            tmp=np.loadtxt(line)
            #print 'tmp.shape=',tmp.shape

            N_T=tmp.shape[0]
            Nd_T=tmp.shape[1]

            T_z=np.zeros((Nd_T,N_T))
            T_z=T_z.reshape(Nd_T,N_T)
            for i in range(Nd_T):
                for j in range(N_T):
                    T_z[i][j]=float(tmp[j][i])

            T_mat=np.zeros(N_T) 
            T_mat=T_mat.reshape(1,N_T) 
            T_t_mat=np.zeros(N_T)
            T_t_mat=T_t_mat.reshape(1,N_T)

            for n in range(N_T):
                T_mat[0][n]=T_z[1][n]
                T_t_mat[0][n]=T_z[0][n]

            #(T_ON_mat,T_OFF_mat,T_duration_mat,T_ON_ini_mat,T_ON_end_mat,T_OFF_ini_mat,T_OFF_end_mat)=duration_time_kai(N_T,T_ON_line,T_OFF_line,T_mat,T_t_mat)
            (T_ON_mat,T_OFF_mat,T_duration_mat,T_ON_ini_mat,T_ON_end_mat,T_OFF_ini_mat,T_OFF_end_mat)=duration_time(N_T,T_x,T_mat,T_t_mat)
            #for i in range(T_ON_mat.shape[1]):
            #    print(T_ON_ini_mat[i])
            #    print(T_ON_end_mat[i])
            #for i in range(T_OFF_mat.shape[1]): 
            #    print(T_OFF_ini_mat[i])
            #    print(T_OFF_end_mat[i]) 

            W_tot_mat=count_events(N_eve,num_mat,ON_mat)
            W_ans_mat=count_events(N_eve,num_mat,T_ON_mat)
            V_tot_mat=count_events(N_eve,num_mat,OFF_mat)
            V_ans_mat=count_events(N_eve,num_mat,T_OFF_mat)

            W_pre_sum_tot_mat=W_pre_sum_tot_mat+W_tot_mat
            (W_tot_mat,W_suc_mat)=success_events(N_eve,T_sig,num_mat,W_tot_mat,ON_ini_mat,ON_end_mat,T_ON_ini_mat,T_ON_end_mat,ON_mat,T_ON_mat)
            (V_tot_mat,V_suc_mat)=success_events(N_eve,T_sig,num_mat,V_tot_mat,OFF_ini_mat,OFF_end_mat,T_OFF_ini_mat,T_OFF_end_mat,OFF_mat,T_OFF_mat)

            W_sum_tot_mat=W_sum_tot_mat+W_tot_mat
            W_sum_ans_mat=W_sum_ans_mat+W_ans_mat
            W_sum_suc_mat=W_sum_suc_mat+W_suc_mat
            V_sum_tot_mat=V_sum_tot_mat+V_tot_mat
            V_sum_ans_mat=V_sum_ans_mat+V_ans_mat
            V_sum_suc_mat=V_sum_suc_mat+V_suc_mat

            if out_ans=='ON':
                wrt_T_ON_m(m,T_ON_mat)
                wrt_T_OFF_m(m,T_OFF_mat)
            #if m==999:
            #   print 'm=',m 
            #   print 'Total_ON_event:',w_tot
            #   print 'ON_event_success:',w_suc
            #   print 'Total_OFF_event:',v_tot
            #   print 'OFF_event_success:',v_suc
            #elif m==4999:
            #   print 'm=',m 
            #   print 'Total_ON_event:',w_tot
            #   print 'ON_event_success:',w_suc
            #   print 'Total_OFF_event:',v_tot
            #   print 'OFF_event_success:',v_suc
            #elif m==9999:
            #   print 'm=',m 
            #   print 'Total_ON_event:',w_tot
            #   print 'ON_event_success:',w_suc
            #   print 'Total_OFF_event:',v_tot
            #   print 'OFF_event_success:',v_suc

        #wrt_x_line(m,N,t_mat,x)
        if out_ON_OFF=='ON':
            wrt_ON_m(m,ON_mat)
            wrt_OFF_m(m,OFF_mat)
            if out_duration_mat=='ON':
                wrt_duration(m,t_mat,duration_mat)
        #wrt_ON_m(m,T_ON_mat)
        #wrt_OFF_m(m,T_OFF_mat)

    if out_ON_OFF=='ON':
        wrt_ON(M,event_error,event_number)
        wrt_OFF(M,event_error,event_number)
    if answer_check=='ON':
        wrt_ON_hist(N_eve,num_mat,W_sum_ans_mat,W_sum_tot_mat,W_sum_suc_mat)
        wrt_OFF_hist(N_eve,num_mat,V_sum_ans_mat,V_sum_tot_mat,V_sum_suc_mat)
        wrt_success_rate(N_eve,num_mat,W_sum_suc_mat,W_sum_ans_mat,V_sum_suc_mat,V_sum_ans_mat)
        wrt_precision(N_eve,num_mat,W_sum_suc_mat,W_sum_tot_mat,V_sum_suc_mat,V_sum_tot_mat)
    if out_ans=='ON':
        wrt_T_ON(M,event_error,event_number)
        wrt_T_OFF(M,event_error,event_number)

    #print 'W_pre_sum_tot_mat=',W_pre_sum_tot_mat
    #print 'W_sum_tot_mat=',W_sum_tot_mat

t2=time.time() 
elapsed_time=t2-t1
print('# elapsed_time=',elapsed_time)
#for i in range(ON_mat.shape[1]):
   #print'ON_mat[0]['+str(i)+']',ON_mat[0][i]
#for j in range(OFF_mat.shape[1]):
   #print'OFF_mat[0]['+str(j)+']',OFF_mat[0][j]
   
