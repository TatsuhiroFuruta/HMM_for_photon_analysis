#coding: utf-8
import sys,os
import numpy as np
import scipy.stats as ss
import math

def gibbs_sampling_s(K,N,x_mat,s_mat,mu_vec,lamda_vec,pi_vec,A_mat):
     a_vec=np.zeros((K,1))
     s_hat_mat=np.zeros((K,N))
     s_hat_mat=s_hat_mat.reshape(K,N)
     iter_mat=np.zeros((K,N))
     iter_mat=iter_mat.reshape(K,N)
     lamda_log_vec=np.log(lamda_vec)
     lamda_log_vec=lamda_log_vec.reshape(K,1)
     A_log_mat=np.log(A_mat)
     pi_log_vec=np.log(pi_vec)
     pi_log_vec=pi_log_vec.reshape(K,1)
     a_mat=np.zeros((K,N))
     d_mat=np.zeros((K,N))
     for n in range(0,N):
         #print'n=',n
         x=x_mat[:,n]
         a=0
         for k in range(0,K):
             #print'k=',k
             lamda=lamda_vec[k]
             mu=mu_vec[k]
             a=(x-mu)**2
             a=-0.5*lamda*a
             a_vec[k][0]=a
         a_vec=a_vec+0.5*lamda_log_vec
         a_vec=a_vec.flatten()
         a_mat[:,n]=a_vec
         a_vec=a_vec.reshape(K,1)
         if n==0:
             a_vec=a_vec+pi_log_vec
             #print'A_mat=',A_mat
             #print'log(A_mat)=',np.log(A_mat)
             #print'A_vec=',A_vec
             #sys.exit()
             s_vec=s_mat[:,n+1]
             s_vec=s_vec.reshape(K,1)
             #print's_mat=',s_mat
             #print's_vec=',s_vec
             b_vec=np.dot(s_vec.T,A_log_mat) 
             b_vec=b_vec.reshape(K,1)
             a_vec=a_vec+b_vec
             d_vec=pi_log_vec+b_vec
             d_vec=d_vec.flatten()
             d_mat[:,n]=d_vec
         elif n==N-1:
             s_vec=s_mat[:,n-1]
             s_vec=s_vec.reshape(K,1)
             c_vec=np.dot(A_log_mat,s_vec) 
             a_vec=a_vec+c_vec
             d_vec=c_vec
             d_vec=d_vec.flatten()
             d_mat[:,n]=d_vec
         else:
             s_vec=s_mat[:,n+1]
             s_vec=s_vec.reshape(K,1)
             b_vec=np.dot(s_vec.T,A_log_mat) 
             b_vec=b_vec.reshape(K,1)
             a_vec=a_vec+b_vec
             s_vec=s_mat[:,n-1]
             s_vec=s_vec.reshape(K,1)
             c_vec=np.dot(A_log_mat,s_vec) 
             a_vec=a_vec+c_vec
             d_vec=b_vec+c_vec
             d_vec=d_vec.flatten()
             d_mat[:,n]=d_vec
             #sys.exit()
             #print'b=',b
             #sys.exit()
         c=0
         iter_vec=np.exp(a_vec) 
         for k in range(K):
             c=c+iter_vec[k][0]
         iter_vec=iter_vec/c 
         iter_vec=iter_vec.flatten()
         iter_mat[:,n]=iter_vec
         s_hat_vec=np.random.multinomial(1,iter_vec)
         s_hat_mat[:,n]=s_hat_vec
     s_mat=s_hat_mat
     #print's_mat=',s_mat
     return a_mat,d_mat,iter_mat,s_mat

def gibbs_sampling_s_simple_mixing(N,D,K,x_mat,s_mat,lamda_vec,pi_vec,A_mat,alpha,iter_mat_tilde):
     s_hat_mat=np.zeros((K,N))
     s_hat_mat=s_hat_mat.reshape(K,N)
     iter_mat=np.zeros((K,N))
     iter_mat=iter_mat.reshape(K,N)
     a_mat=np.zeros((K,N))
     d_mat=np.zeros((K,N))
     for n in range(0,N):
        #print'n=',n
        c=0
        iter_vec=np.zeros(K)
        iter_vec_tilde=np.zeros(K)
        if n==0:
           for k in range(0,K):
              #print'k=',k
              A_vec=np.log(A_mat[:,k])
              A_vec=A_vec.reshape(K,1)
              #print'A_mat=',A_mat
              #print'log(A_mat)=',np.log(A_mat)
              #print'A_vec=',A_vec
              #sys.exit()
              s_vec=s_mat[:,1]
              s_vec=s_vec.reshape(K,1)
              #print's_mat=',s_mat
              #print's_vec=',s_vec
              a=np.dot(s_vec.T,A_vec)
              #print'a=',a
              #sys.exit()
              x_vec=x_mat[:,0]
              x_vec=x_vec.reshape(D,1)
              lamda=lamda_vec[k]
              p=pi_vec[k]
              #print'x_vec=',x_vec
              #print'lamda=',lamda
              #print'p=',p
              #print'x_vec*np.log(lamda)=',x_vec*np.log(lamda)
              #print'log(p)=',np.log(p)
              b=x_vec*np.log(lamda)-lamda+np.log(p)+a
              a_mat[k][n]=x_vec*np.log(lamda)-lamda
              d_mat[k][n]=np.log(p)+a
              #print'b=',b
              #sys.exit()
              iter_vec[k]=np.exp(b)
              c=c+iter_vec[k]
           iter_vec=iter_vec/c
           iter_mat[:,0]=iter_vec
           #print'iter_vec=',iter_vec
           #sys.exit()
           iter_vec_tilde=iter_mat_tilde[:,n]
           iter_vec=alpha*iter_vec+(1.0-alpha)*iter_vec_tilde
           #print'iter_vec=',iter_vec
           s1_vec=np.random.multinomial(1,iter_vec)
           #sys.exit()
           s_hat_mat[:,0]=s1_vec
           #s_mat[:,0]=s1_vec
        elif n>=1 and N-2>=n:
           for k in range(0,K):
              #print'k=',k
              A_vec=np.log(A_mat[:,k])
              A_vec=A_vec.reshape(K,1)
              #print'A_vec=',A_vec
              s_vec=s_mat[:,n+1]
              s_vec=s_vec.reshape(K,1)
              #print's_vec=',s_vec
              a=np.dot(s_vec.T,A_vec)
              #print'a=',a
              A_vec=np.log(A_mat[k,:])
              A_vec=A_vec.reshape(K,1)
              #print'A_vec=',A_vec
              s_vec=s_mat[:,n-1]
              s_vec=s_vec.reshape(K,1)
              #print's_vec=',s_vec
              a_dash=np.dot(s_vec.T,A_vec)
              #print'a_dash=',a_dash
              x_vec=x_mat[:,n]
              x_vec=x_vec.reshape(D,1)
              #print'x_vec=',x_vec
              lamda=lamda_vec[k]
              b=x_vec*np.log(lamda)-lamda+a+a_dash
              a_mat[k][n]=x_vec*np.log(lamda)-lamda
              d_mat[k][n]=a+a_dash
              #print'lamda=',lamda
              #print'log(lamda)=',np.log(lamda)
              #print'b=',b
              #print'exp(b)=',np.exp(b)
              #sys.exit()
              iter_vec[k]=np.exp(b)
              c=c+iter_vec[k]
           iter_vec=iter_vec/c
           iter_mat[:,n]=iter_vec
           #print'iter_vec=',iter_vec
           #sys.exit()
           iter_vec_tilde=iter_mat_tilde[:,n]
           iter_vec=alpha*iter_vec+(1.0-alpha)*iter_vec_tilde
           sn_vec=np.random.multinomial(1,iter_vec)
           s_hat_mat[:,n]=sn_vec
           #s_mat[:,n]=sn_vec
        elif n==N-1:
           for k in range(0,K):
              A_vec=np.log(A_mat[k,:])
              A_vec=A_vec.reshape(K,1)
              s_vec=s_mat[:,N-2]
              s_vec=s_vec.reshape(K,1)
              a=np.dot(s_vec.T,A_vec)
              x_vec=x_mat[:,N-1]
              x_vec=x_vec.reshape(D,1)
              lamda=lamda_vec[k]
              b=x_vec*np.log(lamda)-lamda+a
              a_mat[k][n]=x_vec*np.log(lamda)-lamda
              d_mat[k][n]=a
              iter_vec[k]=np.exp(b)
              c=c+iter_vec[k]
           iter_vec=iter_vec/c
           iter_mat[:,N-1]=iter_vec
           #print'iter_vec=',iter_vec
           iter_vec_tilde=iter_mat_tilde[:,n]
           iter_vec=alpha*iter_vec+(1.0-alpha)*iter_vec_tilde
           sN_vec=np.random.multinomial(1,iter_vec)
           s_hat_mat[:,N-1]=sN_vec
           #s_mat[:,N-1]=sN_vec
     s_mat=s_hat_mat
     #print's_mat=',s_mat
     #sys.exit()
     return a_mat,d_mat,iter_mat,s_mat

def gibbs_sampling_mu_lamda(K,x_mat,s_mat,new,m,a,b):
     mu_vec=np.zeros(K) 
     lamda_vec=np.zeros(K) 
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
     for k in range(0,K):
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
         #print'# (new_hat,m_hat)=',(new_hat,m_hat)
         #print'# (a_hat,b_hat)=',(a_hat,b_hat)
         new_hat_vec[k][0]=new_hat
         m_hat_vec[k][0]=m_hat
         a_hat_vec[k][0]=a_hat
         b_hat_vec[k][0]=b_hat
     new_hat_vec=new_hat_vec.flatten()
     m_hat_vec=m_hat_vec.flatten()
     a_hat_vec=a_hat_vec.flatten()
     b_hat_vec=b_hat_vec.flatten()
     for k in range(0,K):
         #print '# k=',k
         new_hat=new_hat_vec[k]
         m_hat=m_hat_vec[k]
         a_hat=a_hat_vec[k]
         b_hat=b_hat_vec[k]
         b_hat_inv=b_hat**(-1)
         lamda=np.random.gamma(a_hat,b_hat_inv)
         lamda_vec[k]=lamda
         lamda_inv=new_hat*lamda
         lamda_inv=lamda_inv**(-1)
         lamda_inv=np.sqrt(lamda_inv)
         mu=np.random.normal(m_hat,lamda_inv)
         mu_vec[k]=mu
         #print'b_hat=',b_hat
         #sys.exit()
     #sys.exit()
     return mu_vec,lamda_vec

def gibbs_sampling_pi(s_mat,alpha_vec):
     s_vec=s_mat[:,0]
     alpha_hat_vec=s_vec+alpha_vec
     pi_vec=np.random.dirichlet(alpha_hat_vec)
     return pi_vec

def gibbs_sampling_A(K,N,s_mat,beta_mat):
     A_mat=np.zeros((K,K))
     beta_hat_mat=np.zeros((K,K))
     beta_hat_mat=beta_hat_mat.reshape(K,K)
     s_1=s_mat[:,1:N]
     s_1=s_1.reshape(K,N-1) 
     s_0=s_mat[:,0:N-1]
     s_0=s_0.reshape(K,N-1)
     beta_hat_mat=np.dot(s_1,s_0.T)+beta_mat
     for j in range(K):
         beta_hat_vec=beta_hat_mat[:,j]
         A_vec=np.random.dirichlet(beta_hat_vec)
         A_mat[:,j]=A_vec
     #print'# (beta_hat_mat[0][0],beta_hat_mat[0][1],beta_hat_mat[1][0],beta_hat_mat[1][1])=',(beta_hat_mat[0][0],beta_hat_mat[0][1],beta_hat_mat[1][0],beta_hat_mat[1][1])
     #sys.exit()
     return A_mat
 
def blocking_gibbs_sampling_s(N,D,K,x_mat,mu_vec,lamda_vec,pi_vec,A_mat):
    s_mat=np.zeros((K,N))
    f_mat=np.zeros((K,N))
    b_mat=np.zeros((K,N))
    p_mat=np.zeros((K,N))
    b_last_vec=np.ones(K)

    f_vec=np.zeros(K)
    c=0.0
    x=x_mat[:,0]
    for k in range(K):
        mu=mu_vec[k]
        lamda=lamda_vec[k]
        lamda_inv=1.0/lamda
        sigma=np.sqrt(lamda_inv)
        g=ss.norm.pdf(x=x,loc=mu,scale=sigma)
        p=pi_vec[k]
        f=g*p
        f_vec[k]=f
        c=c+f
    f_hat_vec=f_vec/c
    f_mat[:,0]=f_hat_vec
        #print '# f=',f
        #print '# f_mat[k][0]=',f_mat[k][0]

    for n in range(1,N):
        x=x_mat[:,n]
        c=0.0
        f_pre_vec=f_mat[:,n-1]
        #print '# f_pre_vec.shape=',f_pre_vec.shape
        for k in range(K):
            mu=mu_vec[k]
            lamda=lamda_vec[k]
            lamda_inv=1.0/lamda
            sigma=np.sqrt(lamda_inv)
            g=ss.norm.pdf(x=x,loc=mu,scale=sigma)
            A_vec=A_mat[k,:]
            #print '# A_vec.shape=',A_vec.shape
            C=np.dot(A_vec.T,f_pre_vec)
            #print '# C.shape=',C.shape
            f=g*C
            #print '# f.shape=',f.shape
            f_vec[k]=f
            c=c+f
        f_hat_vec=f_vec/c
        f_mat[:,n]=f_hat_vec
    
    #print '# f_mat[:,0]=',f_mat[:,0]
    #print '# f_mat[:,1]=',f_mat[:,1]
    #print '# f_mat[:,2]=',f_mat[:,2]
    #print '# f_mat[:,3]=',f_mat[:,3]
    #print '# f_mat[:,4]=',f_mat[:,4]
    #print '# f_mat[:,N-1]=',f_mat[:,N-1]

    d=0.0
    d=b_last_vec.sum()
    b_hat_vec=b_last_vec/d
    b_mat[:,N-1]=b_hat_vec

    b_vec=np.zeros(K)
    for n in range(N-2,-1,-1):
        d=0.0
        x=x_mat[:,n+1]
        for k in range(K):
            b=0.0
            for i in range(K):
                mu=mu_vec[i]
                lamda=lamda_vec[i]
                lamda_inv=1.0/lamda
                sigma=np.sqrt(lamda_inv)
                g=ss.norm.pdf(x=x,loc=mu,scale=sigma)
                A=A_mat[i][k]
                b_back=b_mat[i][n+1]
                D=g*A*b_back
                b=b+D
            b_vec[k]=b
            d=d+b
        b_hat_vec=b_vec/d
        b_mat[:,n]=b_hat_vec
    
    p_vec=np.zeros(K)
    for n in range(N):
        p_sum=0.0
        for k in range(K):
            f=f_mat[k][n]
            b=b_mat[k][n]
            p=b*f
            p_vec[k]=p
            p_sum=p_sum+p
        p_vec=p_vec/p_sum
        p_mat[:,n]=p_vec
        s_vec=np.random.multinomial(1,p_vec)
        s_mat[:,n]=s_vec
    return p_mat,s_mat 

     
