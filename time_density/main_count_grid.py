#coding: utf-8
#(c) Tatsuhiro Furuta
import sys,os
import numpy as np
import argparse
import pandas as pd
from input_file import *

parser=argparse.ArgumentParser()
parser.add_argument('--M',type=int,default=0,help='file numbers of s_mat-ON.csv(default: 0)') 
parser.add_argument('--s_b',type=float,default=0.5,help='The boundary between ON and OFF in hidden variable.(default: 0.5)') 
parser.add_argument('--T_ON_line',type=float,default=0.5,help='The boundary between ON and OFF in answer data.(default: 0.5)') 

args = parser.parse_args()
M=args.M
s_b=args.s_b
T_ON_line=args.T_ON_line

print('M=',M)
print('s_b=',s_b)
print('T_ON_line=',T_ON_line)

N_mat=np.zeros((M,1))
C_mat=np.zeros((M,1))
C_R_mat=np.zeros((M,1))

for m in range(M):
    file_name_s='s_mat-ON-'+str(m)+'.csv'
    (Ns,s_mat)=input_file(file_name_s)
    file_name_T='T_no_noise-'+str(m)+'.csv'
    (N_T,T_N_mat)=input_file(file_name_T)

    if Ns!=N_T:
        print('m=',m)
        print('Ns=',Ns)
        print('N_T=',N_T)
        print('The program doesn\'t work due to Ns!=N_T.')
        sys.exit()
    else:
        N=N_T
        s_b_mat=np.where(s_mat>=s_b,1,0)
        #print(s_b_mat.shape)
        T_Nb_mat=np.where(T_N_mat>=T_ON_line,1,0)
        #print(T_Nb_mat.shape)
        A_C_mat=np.where(T_Nb_mat==s_b_mat,1,0)
        #print(A_C_mat.shape)
        
        C=A_C_mat.sum()
        C_R=float(C)/float(N)
        
        N_mat[m][0]=N
        C_mat[m][0]=C
        C_R_mat[m][0]=C_R

T_a_mat=np.concatenate([N_mat,C_mat,C_R_mat],1)
print('T_a_mat.shape=',T_a_mat.shape)

df=pd.DataFrame(T_a_mat,columns=['Total of grid','Total of correct grid','correct ratio'])
df['Total of grid']=df['Total of grid'].map(lambda x: int(x))
df['Total of correct grid']=df['Total of correct grid'].map(lambda x: int(x))
df.to_csv('grid_number_in_'+str(M)+'_data_set.csv',sep='\t')
df.to_excel('grid_number_in_'+str(M)+'_data_set.xlsx')

C_R_mean=np.mean(C_R_mat)
print('The mean of correct ratio=',C_R_mean)
