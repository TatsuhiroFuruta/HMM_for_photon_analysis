#coding: utf-8
#(c) Tatsuhiro Furuta
import sys,os
import numpy as np

def num(sim_eva,t_grid,N_eve,tau_max,tau_min):
    num_mat=np.zeros(N_eve+1)
    if sim_eva==True:
        for i in range(N_eve+1):
            num_mat[i]=t_grid*10**i
    else:
        num_mat[0]=tau_min
        for i in range(1,N_eve):
            num_mat[i]=t_grid*10**(i-1)
        num_mat[N_eve]=tau_max
    return num_mat

def Generate_T_sig_mat(T_sig,N_eve):
    T_sig_mat=np.zeros(N_eve)
    for i in range(N_eve):
        T_sig_mat[i]=T_sig*10**i
    return T_sig_mat

def count_events(N_eve,num_mat,T_mat):
    X_mat=np.zeros((1,N_eve+2))
    for i in range(T_mat.shape[1]):
        if T_mat[0][i] < num_mat[0]:
            X_mat[0][0]=X_mat[0][0]+1
        for j in range(N_eve):
            if num_mat[j] <= T_mat[0][i] < num_mat[j+1]:
                X_mat[0][j+1]=X_mat[0][j+1]+1
        if T_mat[0][i] >= num_mat[N_eve]:
            X_mat[0][N_eve+1]=X_mat[0][N_eve+1]+1
    return X_mat

def success_events(N_eve,T_sig,num_mat,C_mat,T_ini_mat,T_end_mat,T_ans_ini_mat,T_ans_end_mat,T_mat,T_ans_mat):
    Y_mat=np.zeros((1,N_eve+2))
    T_sig_hat=0.5*T_sig
    #print'T_sig_hat=',T_sig_hat
    for i in range(T_ans_mat.shape[1]):
        T_ans=T_sig_hat*T_ans_mat[0][i]
        for j in range(T_mat.shape[1]):
            if 0.0 <= T_ans_mat[0][i] < num_mat[0]:
                if np.abs(T_ans_ini_mat[i]-T_ini_mat[j]) <= T_ans and np.abs(T_ans_end_mat[i]-T_end_mat[j]) <= T_ans:
                    Y_mat[0][0]=Y_mat[0][0]+1
                    if num_mat[0] < T_mat[0][j] < num_mat[1]:
                        C_mat[0][1]=C_mat[0][1]-1
                        C_mat[0][0]=C_mat[0][0]+1
            elif T_ans_mat[0][i] >= num_mat[N_eve]:
                if np.abs(T_ans_ini_mat[i]-T_ini_mat[j]) <= T_ans and np.abs(T_ans_end_mat[i]-T_end_mat[j]) <= T_ans:
                    Y_mat[0][N_eve+1]=Y_mat[0][N_eve+1]+1
                    if num_mat[N_eve-1] < T_mat[0][j] < num_mat[N_eve]:
                        C_mat[0][N_eve]=C_mat[0][N_eve]-1
                        C_mat[0][N_eve+1]=C_mat[0][N_eve+1]+1
            
            else:
                for k in range(N_eve):
                    if num_mat[k] <= T_ans_mat[0][i] < num_mat[k+1]:
                        if np.abs(T_ans_ini_mat[i]-T_ini_mat[j]) <= T_ans and np.abs(T_ans_end_mat[i]-T_end_mat[j]) <= T_ans:
                            Y_mat[0][k+1]=Y_mat[0][k+1]+1
                            #if k==1:
                            #    print'T_ans_mat[0][i]=',T_ans_mat[0][i]
                            #    print'T_mat[0][j]=',T_mat[0][j]
                            if k==0:
                                if 0.0 <= T_mat[0][j] < num_mat[k]:
                                    C_mat[0][k]=C_mat[0][k]-1
                                    C_mat[0][k+1]=C_mat[0][k+1]+1
                                elif num_mat[k+1] <= T_mat[0][j] < num_mat[k+2]:
                                    C_mat[0][k+2]=C_mat[0][k+2]-1
                                    C_mat[0][k+1]=C_mat[0][k+1]+1
                            elif k==N_eve-1:
                                if num_mat[k-1] <= T_mat[0][j] < num_mat[k]:
                                    C_mat[0][k]=C_mat[0][k]-1
                                    C_mat[0][k+1]=C_mat[0][k+1]+1
                                elif num_mat[k+1] <= T_mat[0][j]:
                                    C_mat[0][k+2]=C_mat[0][k+2]-1
                                    C_mat[0][k+1]=C_mat[0][k+1]+1
                            else:
                                if num_mat[k-1] <= T_mat[0][j] < num_mat[k]:
                                    C_mat[0][k]=C_mat[0][k]-1
                                    C_mat[0][k+1]=C_mat[0][k+1]+1
                                    #if k==4:
                                    #    print '(i,T_ans_mat[0][i])=',(i,T_ans_mat[0][i])
                                    #    print '(j,T_mat[0][j])=',(j,T_mat[0][j])
                                #if k==4:
                                #    if T_mat[0][j] < num_mat[k]:
                                #        print'T_ans_mat[0][i]=',T_ans_mat[0][i]
                                #        print'T_mat[0][j]=',T_mat[0][j]
                                elif num_mat[k+1] <= T_mat[0][j] < num_mat[k+2]:
                                    C_mat[0][k+2]=C_mat[0][k+2]-1
                                    C_mat[0][k+1]=C_mat[0][k+1]+1
    return C_mat,Y_mat
