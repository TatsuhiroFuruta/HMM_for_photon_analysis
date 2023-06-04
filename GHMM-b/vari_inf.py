#coding: utf-8
# (c) Tatsuhiro Furuta
#from __future__ import print_function
import sys,os
import numpy as np
import scipy.special as ss

def variational_inference_s(K,N,x_mat,s_mat,lamda_vec,lamda_log_vec,pi_log_vec,A_log_mat,mu_lamda_vec,mu_pow_lamda_vec):
     #print '# x_mat.shape=',x_mat.shape
     a_vec=np.zeros((K,1))
     s_hat_mat=np.zeros((K,N))
     iter_mat=np.zeros((K,N))
     pi_log_vec=pi_log_vec.reshape(K,1)
     #A_log_mat=A_log_mat.reshape(K,K)
     a_mat=np.zeros((K,N))
     d_mat=np.zeros((K,N))
     for n in range(N):
         #print '# n=',n 
         #print '# lamda_vec=',lamda_vec
         #print '# lamda_log_vec=',lamda_log_vec
         #print '# x_mat[n]=',x_mat[n]
         #print '# x_mat[n].shape=',x_mat[n].shape
         x=x_mat[0][n]
         #print 'x=',x
         x_pow=np.power(x,2)
         #print 'x_pow=',x_pow
         #print '# x=',x
         a=0
         for k in range(K):
             lamda=lamda_vec[k]
             lamda_log=lamda_log_vec[k]
             mu_lamda=mu_lamda_vec[k]
             mu_pow_lamda=mu_pow_lamda_vec[k]
             #print 'mu_lamda=',mu_lamda
             #print 'mu_pow_lamda=',mu_pow_lamda
             a=-0.5*lamda*x_pow
             #print 'a=',a
             a=a+mu_lamda*x
             #print 'a=',a
             a=a-0.5*mu_pow_lamda
             #print 'a=',a
             a=a+0.5*lamda_log
             #print 'a=',a
             a_vec[k][0]=a
         a_vec=a_vec.flatten()
         a_mat[:,n]=a_vec
         a_vec=a_vec.reshape(K,1)
         #print '# a_vec=',a_vec
         #sys.exit()
         if n==0:
             #print '# pi_log_vec=',pi_log_vec
             a_vec=a_vec+pi_log_vec
             s=s_mat[:,n+1]
             s=s.reshape(K,1)
             #print '# s_vec=',s
             #print '# A_log_mat=',A_log_mat
             b_vec=np.dot(s.T,A_log_mat)
             b_vec=b_vec.reshape(K,1)
             #print '# b_vec=',b_vec
             a_vec=a_vec+b_vec
             #print '# a_vec=',a_vec
             #sys.exit()
             b_vec=b_vec.flatten()
             pi_log_vec=pi_log_vec.flatten()
             d_mat[:,n]=b_vec+pi_log_vec
             pi_log_vec=pi_log_vec.reshape(K,1)
         elif n==N-1:
             s=s_mat[:,n-1]
             s=s.reshape(K,1)
             #print '# s_vec=',s
             #print '# A_log_mat=',A_log_mat
             c_vec=np.dot(A_log_mat,s)
             #print '# c_vec=',c_vec
             a_vec=a_vec+c_vec
             #print '# a_vec=',a_vec
             c_vec=c_vec.flatten()
             d_mat[:,n]=c_vec
         else:    
             s=s_mat[:,n+1]
             s=s.reshape(K,1)
             #print '# s_vec=',s
             #print '# A_log_mat=',A_log_mat
             b_vec=np.dot(s.T,A_log_mat)
             b_vec=b_vec.reshape(K,1)
             #print '# b_vec=',b_vec
             a_vec=a_vec+b_vec
             #print '# a_vec=',a_vec
             s=s_mat[:,n-1]
             s=s.reshape(K,1)
             #print '# s_vec=',s
             #print '# A_log_mat=',A_log_mat
             c_vec=np.dot(A_log_mat,s)
             #print '# c_vec=',c_vec
             a_vec=a_vec+c_vec
             #print '# a_vec=',a_vec
             b_vec=b_vec.flatten()
             c_vec=c_vec.flatten()
             d_vec=b_vec+c_vec
             d_mat[:,n]=d_vec
         c=0
         iter_vec=np.zeros((K,1))
         iter_vec=np.exp(a_vec)
         for k in range(K):
             #iter_vec[k][0]=np.exp(a_vec[k][0])
             c=c+iter_vec[k][0]
         iter_vec=iter_vec/c
         iter_vec=iter_vec.flatten()
         #print '# iter_vec=',iter_vec
         #print '# iter_vec.shape=',iter_vec.shape
         #print '# s_mean_mat[:,n].shape=',s_mean_mat[:,n].shape
         iter_mat[:,n]=iter_vec
         s_hat_mat[:,n]=iter_vec
         #sys.exit()
     s_mat=s_hat_mat
     return a_mat,d_mat,iter_mat,s_mat

def variational_inference_mu_lamda(K,x_mat,s_mat,new,m,a,b):
     mu_vec=np.zeros(K) 
     lamda_vec=np.zeros(K) 
     lamda_log_vec=np.zeros(K) 
     mu_lamda_vec=np.zeros(K) 
     mu_pow_lamda_vec=np.zeros(K) 
     #x_mat=x_mat.reshape(1,N)
     new_hat_vec=np.zeros((K,1))
     m_hat_vec=np.zeros((K,1))
     a_hat_vec=np.zeros((K,1))
     b_hat_vec=np.zeros((K,1))
     new_hat_vec=np.sum(s_mat,axis=1)
     new_hat_vec=new_hat_vec.reshape(K,1)
     m_hat_vec=np.dot(s_mat,x_mat.T)
     a_hat_vec=np.sum(s_mat,axis=1)
     a_hat_vec=0.5*a_hat_vec
     a_hat_vec=a_hat_vec.reshape(K,1)
     x_power_mat=np.power(x_mat,2)
     b_hat_vec=np.dot(s_mat,x_power_mat.T)
     b_hat_vec=0.5*b_hat_vec
     for k in range(K):
         #print new_hat_vec[k].shape
         new_hat=new_hat_vec[k][0]
         m_hat=m_hat_vec[k][0]
         a_hat=a_hat_vec[k][0]
         b_hat=b_hat_vec[k][0]
         new_hat=new_hat+new
         m_hat=m_hat+new*m
         m_hat=m_hat/new_hat
         a_hat=a_hat+a
         b_hat=b_hat+0.5*new*m**2
         b_hat=b_hat-0.5*new_hat*m_hat**2
         b_hat=b_hat+b
         new_hat_vec[k][0]=new_hat
         m_hat_vec[k][0]=m_hat
         a_hat_vec[k][0]=a_hat
         b_hat_vec[k][0]=b_hat
     new_hat_vec=new_hat_vec.flatten()
     m_hat_vec=m_hat_vec.flatten()
     a_hat_vec=a_hat_vec.flatten()
     b_hat_vec=b_hat_vec.flatten()
     for k in range(K):
         #print '# k=',k
         new_hat=new_hat_vec[k]
         m_hat=m_hat_vec[k]
         a_hat=a_hat_vec[k]
         b_hat=b_hat_vec[k]
         mu_mean=m_hat
         lamda_mean=a_hat/b_hat
         lamda_log_mean=ss.digamma(a_hat)-np.log(b_hat)
         mu_lamda_mean=(a_hat*m_hat)/b_hat
         mu_pow_lamda_mean=(a_hat*m_hat**2)/b_hat+new_hat**(-1)
         mu_vec[k]=mu_mean
         lamda_vec[k]=lamda_mean
         lamda_log_vec[k]=lamda_log_mean
         mu_lamda_vec[k]=mu_lamda_mean 
         mu_pow_lamda_vec[k]=mu_pow_lamda_mean 
         #print '# lamda_vec['+str(k)+']=',lamda_vec[k]
         #print '# ss.digamma(a_hat_vec['+str(k)+'][0])=',ss.digamma(a_hat_vec[k][0])
         #print '# np.log(b_hat_vec['+str(k)+'][0])=',np.log(b_hat_vec[k][0])
         #print '# lamda_log_vec['+str(k)+']=',lamda_log_vec[k]
     #print'# a_hat_vec=',a_hat_vec
     #print'# b_hat_vec=',b_hat_vec
     #sys.exit()
     return mu_vec,lamda_vec,lamda_log_vec,mu_lamda_vec,mu_pow_lamda_vec

def variational_inference_pi(K,s_mat,alpha_vec):
     pi_vec=np.zeros(K) 
     pi_log_vec=np.zeros(K) 
     #alpha_vec=alpha_vec.reshape(K,1)
     alpha_hat_vec=np.zeros(K)
     a=0
     s=s_mat[:,0]
     #s=s.reshape(K,1)
     #print '# s=',s
     alpha_hat_vec=s+alpha_vec
     for k in range(K):
         alpha_hat=alpha_hat_vec[k]
         a=a+alpha_hat
     #print '# alpha_hat_vec=',alpha_hat_vec
     #a=alpha_hat_vec.sum(axis=0)
     #print '# a=',a
     pi_vec=alpha_hat_vec/a
     pi_vec=pi_vec.flatten()
     for k in range(K):
         alpha_hat=alpha_hat_vec[k]
         pi_log_vec[k]=ss.digamma(alpha_hat)-ss.digamma(a)
     return pi_vec,pi_log_vec
     
def variational_inference_A(K,N,s_mat,beta_mat):
     A_mat=np.zeros((K,K)) 
     A_log_mat=np.zeros((K,K)) 
     s_1=s_mat[:,1:N]
     s_1=s_1.reshape(K,N-1) 
     s_0=s_mat[:,0:N-1]
     s_0=s_0.reshape(K,N-1)
     beta_hat_mat=np.zeros((K,K))
     beta_hat_mat=np.dot(s_1,s_0.T)+beta_mat 
     for k in range(K):
         beta_vec=beta_hat_mat[:,k]
         beta_vec=beta_vec.reshape(K,1)
         b=beta_vec.sum(axis=0)
         A_mean=beta_vec/b
         A_mean=A_mean.flatten()
         #print '# A_mean.shape=',A_mean.shape
         A_mat[:,k]=A_mean
         for i in range(K):
             c=beta_vec[i][0]
             A_log_mat[i][k]=ss.digamma(c)-ss.digamma(b)
     #print'# (beta_hat_mat[0][0],beta_hat_mat[0][1],beta_hat_mat[1][0],beta_hat_mat[1][1])=',(beta_hat_mat[0][0],beta_hat_mat[0][1],beta_hat_mat[1][0],beta_hat_mat[1][1])
     return A_mat,A_log_mat
     
