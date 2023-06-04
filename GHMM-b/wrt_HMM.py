#coding: utf-8
#(c) Kazuma Nakamura,Tatsuhiro Furuta	
import sys,os 
import numpy as np
import matplotlib.pyplot as plt #use /usr/bin/python2.6 main.py
import matplotlib.cm as cm #use /usr/bin/python2.6 main.py

def s_ON_OFF(N,x_mat,s_mat,iter_mat):
     X=0;Y=0;Z=0;W=0;a=0.5
     for n in range(N):
         if s_mat[0][n]>=a:
             X=X+x_mat[0][n]
             Z=Z+1
         else:
             Y=Y+x_mat[0][n]
             W=W+1

     if Z==0:
         Z=1
     elif W==0:
         W=1

     x_mean_mat=np.zeros((2,1))
     x_mean_mat=x_mean_mat.reshape(2,1)

     x_mean_mat[0][0]=X/Z
     x_mean_mat[1][0]=Y/W

     #print'x_mean_mat=',x_mean_mat
     s_ON_mat=np.zeros(N)
     s_ON_mat=s_ON_mat.reshape(1,N)
     s_OFF_mat=np.zeros(N)
     s_OFF_mat=s_OFF_mat.reshape(1,N)
     #print s_OFF_mat.shape

     if x_mean_mat[0][0]<=x_mean_mat[1][0]:
         s_OFF_mat=iter_mat[0,:]
         s_ON_mat=iter_mat[1,:] 
     else:
         s_OFF_mat=iter_mat[1,:]
         s_ON_mat=iter_mat[0,:] 
     
     return s_OFF_mat,s_ON_mat
     
def wrt_x_mat(N,D,x_time_mat,x_mat):  
     f=open('x_mat.csv','w')
     for i in range(N): 
         f.write(str(x_time_mat[0,i]))
         f.write('\t') 
         f.write(str(x_mat[0,i]))      
         f.write('\n')
     f.close()

def wrt_x_restored_mat(N,x_time_mat,x_mat,x_mean,x_sd):  
    x_mat=x_mat.flatten()
    x_mat=x_sd*x_mat
    x_mat=x_mat+x_mean
    x_mat=x_mat.reshape(1,N)
    f=open('x_restored_mat.csv','w')
    for i in range(N): 
        f.write(str(x_time_mat[0,i]))
        f.write('\t') 
        f.write(str(x_mat[0,i]))      
        f.write('\n')
    f.close()

def wrt_s_ave(N,K,calc,x_time_mat,x_mat,s_mat,iter_mat): 
     if K==2:
         (s_OFF_mat,s_ON_mat)=s_ON_OFF(N,x_mat,s_mat,iter_mat)
         s_ON_mat=s_ON_mat.reshape(1,N)
         s_OFF_mat=s_OFF_mat.reshape(1,N)
         #print s_ON_mat.shape
         #print s_OFF_mat.shape
         #for n in range(N):
            #print s_ON_mat[0,n] 
         f=open('s_mat-ON-'+str(calc)+'.csv','w')
         for n in range(N): 
             #if n%10 == 0:
             f.write(str(x_time_mat[0,n]))
             f.write('\t') 
             f.write(str(s_ON_mat[0,n]))      
             f.write('\n')
             #else: 
                 #continue 
         f.close()
         f=open('s_mat-OFF-'+str(calc)+'.csv','w')
         for n in range(N): 
             #if n%10 == 0:
             f.write(str(x_time_mat[0,n]))
             f.write('\t') 
             f.write(str(s_OFF_mat[0,n]))      
             f.write('\n')
             #else: 
                 #continue 
         f.close()
     for k in range(K):
         f=open('s_mat-'+str(k)+'-'+str(calc)+'.csv','w')
         for n in range(N): 
             #if n%10 == 0:
             f.write(str(x_time_mat[0,n]))
             f.write('\t') 
             f.write(str(iter_mat[k,n]))      
             f.write('\n')
             #else: 
                 #continue 
         f.close()

def wrt_x_renormalize_mat(N,D,x_renormalize_time_mat,x_renormalize_mat):  
     f=open('x_renormalize_mat.csv','w')
     for i in range(N): 
         f.write(str(x_renormalize_time_mat[0,i]))
         f.write('\t') 
         f.write(str(x_renormalize_mat[0,i]))      
         f.write('\n')
     f.close()

def wrt_s_ave_itr(itr,N,K,calc,x_time_mat,x_mat,s_mat,iter_mat):  
     if K==2:
         (s_OFF_mat,s_ON_mat)=s_ON_OFF(N,x_mat,s_mat,iter_mat)
         s_ON_mat=s_ON_mat.reshape(1,N)
         s_OFF_mat=s_OFF_mat.reshape(1,N)
         #for n in range(N):
            #print s_ON_mat[0,n] 
         f=open('s_mat-ON-itr'+str(itr)+'-'+str(calc)+'.csv','w')
         for n in range(N): 
             #if n%10 == 0:
             f.write(str(x_time_mat[0,n]))
             f.write('\t') 
             f.write(str(s_ON_mat[0,n]))      
             f.write('\n')
             #else:
                 #continue 
         f.close()
         f=open('s_mat-OFF-itr'+str(itr)+'-'+str(calc)+'.csv','w')
         for n in range(N): 
             #if n%10 == 0:
             f.write(str(x_time_mat[0,n]))
             f.write('\t') 
             f.write(str(s_OFF_mat[0,n]))      
             f.write('\n')
             #else:
                 #continue 
         f.close()
     for k in range(K):
         f=open('s_mat-'+str(k)+'-itr'+str(itr)+'-'+str(calc)+'.csv','w')
         for n in range(N): 
             #if n%10 == 0:
             f.write(str(x_time_mat[0,n]))
             f.write('\t') 
             f.write(str(iter_mat[k,n]))      
             f.write('\n')
             #else:
                 #continue 
         f.close()

def wrt_mu(K,MAXITER,calc,mu_mat):
     f=open('mu_mat-'+str(calc)+'.csv','w')
     for itr in range(MAXITER+1): 
         f.write(str(itr))
         for k in range(K):
             f.write('\t') 
             f.write(str(mu_mat[k][itr]))      
         f.write('\n')
     f.close()

def wrt_lamda(K,MAXITER,calc,lamda_mat):
     f=open('lamda_mat-'+str(calc)+'.csv','w')
     for itr in range(MAXITER+1): 
         f.write(str(itr))
         for k in range(K):
             f.write('\t') 
             f.write(str(lamda_mat[k][itr]))      
         f.write('\n')
     f.close()

def wrt_pi(K,MAXITER,calc,pi_mat):
     f=open('pi_mat-'+str(calc)+'.csv','w')
     for itr in range(MAXITER+1): 
         f.write(str(itr))
         for k in range(K):
             f.write('\t') 
             f.write(str(pi_mat[k][itr]))      
         f.write('\n')
     f.close()
      
def wrt_A(K,MAXITER,calc,A_3d_mat):
     f=open('A_3d_mat-'+str(calc)+'.csv','w')
     for itr in range(MAXITER+1): 
         f.write(str(itr))
         for i in range(K):
             for j in range(K):
                 f.write('\t') 
                 f.write(str(A_3d_mat[i][j][itr]))      
         f.write('\n')
     f.close()

def wrt_iter(itr,N,K,calc,iter_mat,a_mat,d_mat):     
     f=open('iter_mat-'+str(itr)+'-'+str(calc)+'.csv','w')
     f.write('data number')
     f.write('\t') 
     f.write('iter_mat[0,:]')
     f.write('\t') 
     f.write('iter_mat[1,:]')
     f.write('\t') 
     f.write('a_mat[0,:]')
     f.write('\t') 
     f.write('a_mat[1,:]')
     f.write('\t') 
     f.write('d_mat[0,:]')
     f.write('\t') 
     f.write('d_mat[1,:]')
     f.write('\n') 
     for n in range(N):
         f.write(str(n))
         f.write('\t') 
         for k in range(K):
             f.write(str(iter_mat[k][n]))
             f.write('\t') 
         for k in range(K):
             f.write(str(a_mat[k][n]))
             f.write('\t') 
         for k in range(K):
             f.write(str(d_mat[k][n]))
             f.write('\t') 
         f.write('\n')
     f.close()
         
def wrt_s_continue(N,K,calc,s_mat_ave):
     f=open('s_mat-'+str(calc)+'.csv','w')
     for n in range(N):
         for k in range(K):
             f.write(str(s_mat_ave[k][n]))              
             f.write('\t') 
         f.write('\n')
     f.close()

def wrt_mu_continue(K,calc,mu_vec):
     f=open('mu_vec-'+str(calc)+'.csv','w')
     for k in range(K):
         f.write(str(mu_vec[k]))              
         f.write('\t') 
     f.close()

def wrt_lamda_continue(K,calc,lamda_vec):
     f=open('lamda_vec-'+str(calc)+'.csv','w')
     for k in range(K):
         f.write(str(lamda_vec[k]))              
         f.write('\t') 
     f.close()

def wrt_pi_continue(K,calc,pi_vec):
     f=open('pi_vec-'+str(calc)+'.csv','w')
     for k in range(K):
         f.write(str(pi_vec[k]))              
         f.write('\t') 
     f.close()
     
def wrt_A_continue(K,calc,A_mat):
     f=open('A_mat-'+str(calc)+'.csv','w')
     for i in range(K):
         for j in range(K):
             f.write(str(A_mat[i][j]))              
             f.write('\t') 
         f.write('\n')
     f.close()
