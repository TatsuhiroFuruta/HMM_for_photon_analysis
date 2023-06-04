#coding: utf-8
#(c) Tatsuhiro Furuta
import sys,os
import numpy as np

def wrt_duration(m,t_mat,duration_mat):
    f=open(str(m)+'-duration_mat.csv','w')
    for n in range(duration_mat.shape[1]):
        f.write(str(t_mat[0][n]))
        f.write('\t')
        f.write(str(duration_mat[0][n]))
        f.write('\n')
    f.close()

def wrt_x_line(m,N,t_mat,x):
    f=open(str(m)+'-x_line.csv','w')
    for n in range(N):
        f.write(str(t_mat[0][n]))
        f.write('\t')
        f.write(str(x))
        f.write('\t')
        f.write(str(-x))
        f.write('\n')
    f.close()

def wrt_ON_m(m,ON_mat):
    f=open(str(m)+'-ON_mat.csv','w')
    for i in range(ON_mat.shape[1]):
        f.write(str(ON_mat[0][i]))
        f.write('\n')
    f.close()

def wrt_OFF_m(m,OFF_mat):
    f=open(str(m)+'-OFF_mat.csv','w')
    for j in range(OFF_mat.shape[1]):
        f.write(str(OFF_mat[0][j]))
        f.write('\n')
    f.close()

def wrt_ON(M,event_error,event_number):
    #event_number=2
    f=open('ON_mat.csv','w')
    for m in range(M):
        g=open(str(m)+'-ON_mat.csv','r')
        line=g.readlines()
        g.close()

        #print 'm=',m
        if line==[]:
            print('# ON:the number of []=',m)
        else:
            tmp=np.loadtxt(line)
            if tmp.shape==():
                if event_error=='ON':
                    if event_number!=1:
                        print('# ON:error number=',m)
                f.write(str(tmp))
                f.write('\n')
            else:
                if event_error=='ON':
                    if tmp.shape[0]!=event_number:
                        print('# ON:error number=',m)
                for n in range(tmp.shape[0]):
                    f.write(str(tmp[n]))
                    f.write('\n')
        #print tmp.shape
        #sys.exit()

    f.close() 

def wrt_OFF(M,event_error,event_number):
    #event_number=2
    f=open('OFF_mat.csv','w')
    for m in range(M):
        #print 'm=',m
        g=open(str(m)+'-OFF_mat.csv','r')
        line=g.readlines()
        g.close()

        if line==[]:
            print('# OFF:the number of []=',m)
        else:
            tmp=np.loadtxt(line)
            #print tmp.shape
            if tmp.shape==():
                if event_error=='ON':
                    if event_number!=1:
                        print('# OFF:error number=',m)
                f.write(str(tmp))
                f.write('\n')
            else:
                if event_error=='ON':
                    if tmp.shape[0]!=event_number:
                        print('# OFF:error number=',m)
                for n in range(tmp.shape[0]):
                    f.write(str(tmp[n]))
                    f.write('\n')
        #sys.exit()

    f.close() 

def wrt_ON_hist(N_eve,num_mat,W_ans_mat,W_tot_mat,W_suc_mat):
    w_ans=W_ans_mat.sum(axis=1)
    w_tot=W_tot_mat.sum(axis=1)
    w_suc=W_suc_mat.sum(axis=1)
    f=open('ON_hist.csv','w')
    f.write('Duration_time')
    f.write('\t')
    f.write('Answer_events')
    f.write('\t')
    f.write('Simulation_events')
    f.write('\t')
    f.write('Success_events')
    f.write('\n')
    f.write('0-'+str(num_mat[0]))
    f.write('\t')
    f.write(str(int(W_ans_mat[0][0])))
    f.write('\t')
    f.write(str(int(W_tot_mat[0][0])))
    f.write('\t')
    f.write(str(int(W_suc_mat[0][0])))
    f.write('\n')
    for i in range(N_eve):
        f.write(str(num_mat[i])+'-'+str(num_mat[i+1]))
        f.write('\t')
        f.write(str(int(W_ans_mat[0][i+1])))
        f.write('\t')
        f.write(str(int(W_tot_mat[0][i+1])))
        f.write('\t')
        f.write(str(int(W_suc_mat[0][i+1])))
        f.write('\n')
    f.write(str(num_mat[N_eve])+'-')
    f.write('\t')
    f.write(str(int(W_ans_mat[0][N_eve+1])))
    f.write('\t')
    f.write(str(int(W_tot_mat[0][N_eve+1])))
    f.write('\t')
    f.write(str(int(W_suc_mat[0][N_eve+1])))
    f.write('\n')
    f.write('\n')
    f.write('Total_of_answer_events:')
    f.write('\t')
    f.write(str(int(w_ans)))
    f.write('\n')
    f.write('Total_of_simulation_events:')
    f.write('\t')
    f.write(str(int(w_tot)))
    f.write('\n')
    f.write('Total_of_success_events:')
    f.write('\t')
    f.write(str(int(w_suc)))
    f.write('\n')
    f.close()

def wrt_OFF_hist(N_eve,num_mat,V_ans_mat,V_tot_mat,V_suc_mat):
    v_ans=V_ans_mat.sum(axis=1)
    v_tot=V_tot_mat.sum(axis=1)
    v_suc=V_suc_mat.sum(axis=1)
    f=open('OFF_hist.csv','w')
    f.write('Duration_time')
    f.write('\t')
    f.write('Answer_events')
    f.write('\t')
    f.write('Simulation_events')
    f.write('\t')
    f.write('Success_events')
    f.write('\n')
    f.write('0-'+str(num_mat[0]))
    f.write('\t')
    f.write(str(int(V_ans_mat[0][0])))
    f.write('\t')
    f.write(str(int(V_tot_mat[0][0])))
    f.write('\t')
    f.write(str(int(V_suc_mat[0][0])))
    f.write('\n')
    for i in range(N_eve):
        f.write(str(num_mat[i])+'-'+str(num_mat[i+1]))
        f.write('\t')
        f.write(str(int(V_ans_mat[0][i+1])))
        f.write('\t')
        f.write(str(int(V_tot_mat[0][i+1])))
        f.write('\t')
        f.write(str(int(V_suc_mat[0][i+1])))
        f.write('\n')
    f.write(str(num_mat[N_eve])+'-')
    f.write('\t')
    f.write(str(int(V_ans_mat[0][N_eve+1])))
    f.write('\t')
    f.write(str(int(V_tot_mat[0][N_eve+1])))
    f.write('\t')
    f.write(str(int(V_suc_mat[0][N_eve+1])))
    f.write('\n')
    f.write('\n')
    f.write('Total_of_answer_events:')
    f.write('\t')
    f.write(str(int(v_ans)))
    f.write('\n')
    f.write('Total_of_simulation_events:')
    f.write('\t')
    f.write(str(int(v_tot)))
    f.write('\n')
    f.write('Total_of_success_events:')
    f.write('\t')
    f.write(str(int(v_suc)))
    f.write('\n')
    f.close()

def wrt_success_rate(N_eve,num_mat,W_suc_mat,W_ans_mat,V_suc_mat,V_ans_mat):
    P_ON_mat=np.zeros(N_eve)
    P_OFF_mat=np.zeros(N_eve)
    for i in range(N_eve):
        if W_ans_mat[0][i+1]==0:
            P_ON_mat[i]=0
        else:
            P_ON_mat[i]=float(W_suc_mat[0][i+1])/float(W_ans_mat[0][i+1])
        if V_ans_mat[0][i+1]==0:
            P_OFF_mat[i]=0
        else:
            P_OFF_mat[i]=float(V_suc_mat[0][i+1])/float(V_ans_mat[0][i+1])
    f=open('success_rate.csv','w')
    f.write('Duration_time')
    f.write('\t')
    f.write('ON')
    f.write('\t')
    f.write('OFF')
    f.write('\n')
    for i in range(N_eve):
        f.write(str(num_mat[i])+'-'+str(num_mat[i+1]))
        f.write('\t')
        f.write(str(P_ON_mat[i]))
        f.write('\t')
        f.write(str(P_OFF_mat[i]))
        f.write('\n')
    f.close() 

def wrt_precision(N_eve,num_mat,W_suc_mat,W_tot_mat,V_suc_mat,V_tot_mat):
    P_ON_mat=np.zeros(N_eve+2)
    P_OFF_mat=np.zeros(N_eve+2)
    for i in range(N_eve+2):
        if W_tot_mat[0][i]==0:
            P_ON_mat[i]=0
        else:
            P_ON_mat[i]=float(W_suc_mat[0][i])/float(W_tot_mat[0][i])
        if V_tot_mat[0][i]==0:
            P_OFF_mat[i]=0
        else:
            P_OFF_mat[i]=float(V_suc_mat[0][i])/float(V_tot_mat[0][i])
    f=open('precision.csv','w')
    f.write('Duration_time')
    f.write('\t')
    f.write('ON')
    f.write('\t')
    f.write('OFF')
    f.write('\n')
    f.write('0-'+str(num_mat[0]))
    f.write('\t')
    f.write(str(P_ON_mat[0]))
    f.write('\t')
    f.write(str(P_OFF_mat[0]))
    f.write('\n')
    for i in range(N_eve):
        f.write(str(num_mat[i])+'-'+str(num_mat[i+1]))
        f.write('\t')
        f.write(str(P_ON_mat[i+1]))
        f.write('\t')
        f.write(str(P_OFF_mat[i+1]))
        f.write('\n')
    f.write(str(num_mat[N_eve])+'-')
    f.write('\t')
    f.write(str(P_ON_mat[N_eve+1]))
    f.write('\t')
    f.write(str(P_OFF_mat[N_eve+1]))
    f.write('\n')
    f.close() 
   
def wrt_T_ON_m(m,T_ON_mat):
    f=open(str(m)+'-T_ON_mat.csv','w')
    for i in range(T_ON_mat.shape[1]):
        f.write(str(T_ON_mat[0][i]))
        f.write('\n')
    f.close()

def wrt_T_OFF_m(m,T_OFF_mat):
    f=open(str(m)+'-T_OFF_mat.csv','w')
    for j in range(T_OFF_mat.shape[1]):
        f.write(str(T_OFF_mat[0][j]))
        f.write('\n')
    f.close()

def wrt_T_ON(M,event_error,event_number):
    #event_number=2
    f=open('T_ON_mat.csv','w')
    for m in range(M):
        g=open(str(m)+'-T_ON_mat.csv','r')
        line=g.readlines()
        g.close()

        #print 'm=',m
        if line==[]:
            print('# ON:the number of []=',m)
        else:
            tmp=np.loadtxt(line)
            if tmp.shape==():
                if event_error=='ON':
                    if event_number!=1:
                        print('# ON:error number=',m)
                f.write(str(tmp))
                f.write('\n')
            else:
                if event_error=='ON':
                    if tmp.shape[0]!=event_number:
                        print('# ON:error number=',m)
                for n in range(tmp.shape[0]):
                    f.write(str(tmp[n]))
                    f.write('\n')
        #print tmp.shape
        #sys.exit()

    f.close() 

def wrt_T_OFF(M,event_error,event_number):
    #event_number=2
    f=open('T_OFF_mat.csv','w')
    for m in range(M):
        #print 'm=',m
        g=open(str(m)+'-T_OFF_mat.csv','r')
        line=g.readlines()
        g.close()

        if line==[]:
            print('# OFF:the number of []=',m)
        else:
            tmp=np.loadtxt(line)
            #print tmp.shape
            if tmp.shape==():
                if event_error=='ON':
                    if event_number!=1:
                        print('# OFF:error number=',m)
                f.write(str(tmp))
                f.write('\n')
            else:
                if event_error=='ON':
                    if tmp.shape[0]!=event_number:
                        print('# OFF:error number=',m)
                for n in range(tmp.shape[0]):
                    f.write(str(tmp[n]))
                    f.write('\n')
        #sys.exit()

    f.close() 
