#coding: utf-8
import sys,os
import numpy as np
import scipy.stats as ss

def first_phase_s(K,N):
     s_mat=np.zeros((K,N))
     s_mat=s_mat.reshape(K,N)
     for n in range(0,N):
         #alpha=np.random.random_sample((K,))
         alpha=np.ones(K)
         #print'alpha.shape',alpha.shape  
         pi=ss.dirichlet.rvs(alpha)
         #print'pi.shape',pi.shape  
         pi=pi.flatten() 
         #print'pi.shape',pi.shape  
         s=ss.multinomial.rvs(1,pi)
         #print's.shape',s.shape  
         #sys.exit() 
         s_mat[:,n]=s
     return s_mat

def first_phase_s_mean(K,N):
     s_mat=0.50*np.ones((K,N))
     return s_mat

def first_phase_lamda(K,a,b):
     lamda_vec=np.zeros(K)
     lamda_vec=lamda_vec.reshape(K,1)
     b=b**(-1)
     for k in range(0,K):
         lamda_vec[k,:]=np.random.gamma(a,b)
     return lamda_vec

def generate_init_mu_lamda(K,new,m,a,b):
     mu_vec=np.zeros(K)
     lamda_vec=np.zeros(K)
     b=b**(-1)
     for k in range(0,K):
         lamda=np.random.gamma(a,b)
         lamda_vec[k]=lamda
         lamda_inv=new*lamda
         lamda_inv=lamda_inv**(-1)
         lamda_inv=np.sqrt(lamda_inv)
         mu=np.random.normal(m,lamda_inv)
         mu_vec[k]=mu
     return mu_vec,lamda_vec

def first_phase_pi(K):
     #alpha=np.random.random_sample((K,))
     alpha=np.ones(K)
     pi_vec=ss.dirichlet.rvs(alpha)
     pi_vec=pi_vec.flatten() 
     return pi_vec

def generate_init_pi(K,alpha_vec):
     pi_vec=ss.dirichlet.rvs(alpha_vec)
     pi_vec=pi_vec.flatten() 
     return pi_vec

def first_phase_A(K,beta_mat):
     A_mat=np.zeros((K,K))
     for k in range(0,K):
         beta_vec=beta_mat[:,k]
         A_vec=np.random.dirichlet(beta_vec)
         A_mat[:,k]=A_vec
     A_mat=A_mat.reshape(K,K)
     return A_mat

def first_phase_beta(K,N):
     D=1.0
     gam=0.01
     #print'# gam=',gam
     print('# gam=',gam)
     C=(1-2*gam)*N+4*(1-gam)*D
     C=C/(4*gam)
     #D=1.0
     beta_mat=np.zeros((K,K))
     for i in range(0,K):
         for j in range(0,K):
            if i==j:
                beta_mat[i,j]=C
            else:
                beta_mat[i,j]=D
     return beta_mat
                
def moving_average(K,N,term,x_time_mat,x_mat): 
     s_mat=np.zeros((K,N))
     x_mean_mat=np.zeros(N)
     I_vec=np.ones(N)
     for i in range(N):
         if i > N-term:
             x_mean_mat[i]=x_mean_mat[N-term]
         else:
             for j in range(term):
                 x_mean_mat[i]=x_mean_mat[i]+x_mat[0][i+j]
             x_mean_mat[i]=x_mean_mat[i]/term
     a=np.max(x_mean_mat)
     b=np.min(x_mean_mat)
     for i in range(N):
         x_mean_mat[i]=(x_mean_mat[i]-b)/(a-b)
     s_mat[0,:]=I_vec-x_mean_mat
     s_mat[1,:]=x_mean_mat
     f=open('s_mat_init.csv','w')
     for i in range(N):
         f.write(str(x_time_mat[0][i]))
         f.write('\t')
         f.write(str(x_mean_mat[i]))
         f.write('\n')
     f.close()
     #for i in range(N):
         #print s_mat[0][i]+s_mat[1][i]
     #sys.exit()
     return s_mat

#def first_phase(CONTINUE,calc,N,K,x_mat):
def first_phase(CONTINUE,calc,N,K,x_time_mat,x_mat,term,first_moving):
     if CONTINUE=='ON':
         f=open('s_mat.csv','r')
         line=f.readlines()
         tmp=np.loadtxt(line) 
         #print'# tmp.shape=',tmp.shape
         #sys.exit()
         s_mat=np.zeros((K,N))
         for n in range(N):
             for k in range(K):
                 s_mat[k][n]=float(tmp[n][k])
         #print s_mat[1][2]
         #sys.exit()
         
         f=open('mu_vec.csv','r') 
         line=f.readlines()
         tmp=np.loadtxt(line) 
         mu_vec=np.zeros(K)
         for k in range(K):
             mu_vec[k]=float(tmp[k])

         f=open('lamda_vec.csv','r') 
         line=f.readlines()
         tmp=np.loadtxt(line) 
         lamda_vec=np.zeros(K)
         for k in range(K):
             lamda_vec[k]=float(tmp[k])

         f=open('pi_vec.csv','r') 
         line=f.readlines()
         tmp=np.loadtxt(line) 
         pi_vec=np.zeros(K)
         for k in range(K):
             pi_vec[k]=float(tmp[k])
       
         f=open('A_mat.csv','r')
         line=f.readlines()
         tmp=np.loadtxt(line) 
         A_mat=np.zeros((K,K))
         for i in range(K):
             for j in range(K):
                 A_mat[i][j]=float(tmp[i][j])
         #print'# (A_mat[0][0],A_mat[0][1],A_mat[1][0],A_mat[1][1])=',(A_mat[0][0],A_mat[0][1],A_mat[1][0],A_mat[1][1])
       
     else:
         if K==2:
             if first_moving=='ON':
                 s_mat=moving_average(K,N,term,x_time_mat,x_mat) 
             else:
                 if calc=='Comp':
                     s_mat=first_phase_s_mean(K,N)
                 else:
                     s_mat=first_phase_s(K,N)
         else:
             if calc=='Comp':
                 s_mat=first_phase_s_mean(K,N)
             else:
                 s_mat=first_phase_s(K,N)
         mu_vec=np.linspace(np.min(x_mat),np.max(x_mat),int(K))
         #lamda_vec=np.linspace(np.min(x_mat),np.max(x_mat),int(K))
         #lamda_vec=np.zeros(K)
         x_sig=np.var(x_mat)
         lamda=x_sig**(-1)
         #lamda_vec=x_sig*np.ones(K)
         lamda_vec=lamda*np.ones(K)
         #lamda_vec[0]=0.5
         #lamda_vec[1]=1.5
         pi_vec=first_phase_pi(K)
         A_mat=np.eye(K)
         for i in range(0,K):
             for j in range(0,K):
                if i==j:
                    A_mat[i,j]=A_mat[i,j]-0.01
                else:
                    A_mat[i,j]=A_mat[i,j]+0.01
         
     return s_mat,mu_vec,lamda_vec,pi_vec,A_mat    

def first_phase_mu_lamda(K,s_mat,mu_vec,lamda_vec,new):  
     new_hat_vec=np.zeros(K)        
     mu_lamda_vec=np.zeros(K)        
     mu_pow_lamda_vec=np.zeros(K)   
     new_hat_vec=np.sum(s_mat,axis=1)        
     for k in range(K):
         mu=mu_vec[k]     
         lamda=lamda_vec[k]     
         new_hat=new_hat_vec[k]     
         new_hat=new_hat+new     
         mu_lamda=mu*lamda     
         mu_pow_lamda=lamda*mu**2+new_hat**(-1)   
     return mu_lamda_vec,mu_pow_lamda_vec  
