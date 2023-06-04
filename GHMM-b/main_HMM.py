#coding: utf-8
#from __future__ import print_function
import sys,os
import numpy as np
from distutils.util import strtobool
from gibbs import *
from vari_inf import *
from data_photon import *
from wrt_HMM  import * 
from first_phase import *
import scipy.stats as ss
import argparse
import time 

t1=time.time() 
parser=argparse.ArgumentParser()
parser.add_argument('--input_file',type=str,default='T.csv',help='The parameter for file name you would like to input.(default: T.csv)')
parser.add_argument('--K',type=int,default=2,help='The number of cluster(default: 2)') 
parser.add_argument('--MAXITER',type=int,default=0,help='The number of Gibbs sampling(default: 0)') 
parser.add_argument('--sampling_term',type=int,default=1,help='The number of Gibbs sampling step for wrt_s_mat(default: 1)') 
parser.add_argument('--s',type=int,default=0,help='The permissible range parameter of sigma (default: 0)')
parser.add_argument('--s_pri',type=str,default='OFF',help='If s_pri==ON, print parameter of x_mean and sigma (default: OFF)')
parser.add_argument('--calc',type=str,default='Gibbs',help='calc is the string value of Gibbs or Block,Comp.Gibbs is gibbs sampling, Block is blocking gibbs sampling and Comp is completely factorized variational inference.(default: Gibbs)')
parser.add_argument('--CONTINUE',type=str,default='OFF',help='If CONTINUE==ON, you can read files, which is s_mat.csv, lamda_vec.csv, pi_vec.csv and A_mat.csv (default: OFF)')
parser.add_argument('--term',type=int,default=2,help='The parameter of first_phase s_mat(if K==2)(default: 2)') 
parser.add_argument('--alpha',type=float,default=0.1,help='The parameter of simple mixing. The value is range from 0.0 to 1.0. (default: 0.1)') 
parser.add_argument('--mixing',type=str,default='OFF',help='If mixing==ON, simple mixing program of gibbs_sampling_s turns on.(default: OFF)')
parser.add_argument('--first_moving',type=str,default='OFF',help='If first_moving==ON and K==2, moving_average program of first_phase turns on.(default: OFF)')
parser.add_argument('--x_int',type=str,default='OFF',help='If x_int==ON, x_integer program of data_photon turns on.(default: OFF)')
parser.add_argument('--x_int_term',type=int,default=100,help='a prameter for x_integer program of data_photon.(default: 100)')
parser.add_argument('--scaling',type=strtobool,default=False,help='If scaling==True, x_scaling function in data_photon.py works.(default: False)')
parser.add_argument('--normalized',type=strtobool,default=False,help='a prameter for x_scaling function in data_photon.py.(default: False)')
parser.add_argument('--x_scale',type=float,default=1.0,help='a prameter for x_scaling function in data_photon.py.(default: 1.0)')
parser.add_argument('--x_shift',type=float,default=0.0,help='a prameter for x_scaling function in data_photon.py.(default: 0.0)')
#parser.add_argument('--A',type=int,default=50,help='The parameter of noise data (defalut: 50)') 
#parser.add_argument('--B',type=int,default=50,help='The parameter of noise data (defalut: 50)') 
#parser.add_argument('--alp',type=int,default=100,help='The parameter of center (default: 100)') 
#parser.add_argument('--bet',type=int,default=105,help='The parameter of center (default: 105)') 
args = parser.parse_args()

input_file=args.input_file
K=args.K
MAXITER=args.MAXITER
sampling_term=args.sampling_term
s=args.s
s_pri=args.s_pri
calc=args.calc
CONTINUE=args.CONTINUE
term=args.term
alpha=args.alpha
mixing=args.mixing
first_moving=args.first_moving
x_int=args.x_int
x_int_term=args.x_int_term
scaling=args.scaling
normalized=args.normalized
x_scale=args.x_scale
x_shift=args.x_shift
#A=args.A
#B=args.B
#alp=args.alp
#bet=args.bet
 
print('# input_file=',input_file)
print('# calc=',calc)
print('# CONTINUE=',CONTINUE)
print('# alpha=',alpha)
print('# mixing=',mixing)
print('# first_moving=',first_moving)
print('# scaling=',bool(scaling))
print('# normalized=',bool(normalized))
print('# x_scale=',x_scale)
print('# x_shift=',x_shift)

#new=1.0
#m=1.0
#a=8.0
#b=1.0
alpha_vec=100*np.ones(K)
#beta_mat=np.random.random_sample((K,K))
#beta_mat=np.zeros((K,K))
#beta_mat=np.array([[100000,100],[100,100000]])
#beta_mat=100*np.ones((K,K))
#beta_mat=first_phase_beta(K)
#beta_mat=beta_mat.reshape(K,K)
#print'# beta_mat=',beta_mat
#print'# (beta_mat[0][0],beta_mat[0][1],beta_mat[1][0],beta_mat[1][1])=',(beta_mat[0][0],beta_mat[0][1],beta_mat[1][0],beta_mat[1][1])
(N,N_renormalize,D,x_time_mat,x_mat,x_renormalize_time_mat,x_renormalize_mat,x_mean,sigma)=mktrain_x_mat(input_file,s)

#print N,D
#sys.exit() 
#print x_mat
#print x_mat.shape
print('# x_int=',x_int)
print('# x_int_term=',x_int_term)
if s==0:
    wrt_x_mat(N,D,x_time_mat,x_mat)
    if scaling==True:
        x_mat=x_scaling(N,normalized,x_scale,x_shift,x_time_mat,x_mat)
    beta_mat=first_phase_beta(K,N)
    beta_mat=beta_mat.reshape(K,K)
    if x_int=='ON':
        x_mat=x_integer(N,x_time_mat,x_mat,x_int_term)
    (x_mat,x_mean_vec,x_sd)=x_norm(N,x_time_mat,x_mat)
    #(s_mat,lamda_vec,pi_vec,A_mat)=first_phase(CONTINUE,calc,N,K,x_mat)
    #(s_mat,mu_vec,lamda_vec,pi_vec,A_mat)=first_phase(CONTINUE,calc,N,K,x_time_mat,x_mat,term,first_moving)
    new=1.0
    #m=np.mean(x_mat)
    m=0.0
    a=1.0
    #b=np.var(x_mat)
    b=1.0
    print('# (new,m)=',(new,m))
    print('# (a,b)=',(a,b))
    (mu_vec,lamda_vec)=generate_init_mu_lamda(K,new,m,a,b)
    pi_vec=generate_init_pi(K,alpha_vec)
    A_mat=first_phase_A(K,beta_mat)
    print('# mu_vec=',mu_vec)
    print('# lamda_vec=',lamda_vec)
    print('# pi_vec=',pi_vec)
    print('# (A_mat[0][0],A_mat[0][1],A_mat[1][0],A_mat[1][1])=',(A_mat[0][0],A_mat[0][1],A_mat[1][0],A_mat[1][1]))
    iter_mat_tilde=0.50*np.ones((K,N))
    #if calc=='Comp':
      #s_mat=first_phase_s_mean(K,N)
    #else:
      #s_mat=first_phase_s(K,N)
    #s_mat=first_phase_s(K,N)
    #sum_s_mat=np.zeros((K,N)) 
    sum_s_mat=np.zeros((K,N)) 
else:
    wrt_x_mat(N,D,x_time_mat,x_mat)
    wrt_x_renormalize_mat(N_renormalize,D,x_renormalize_time_mat,x_renormalize_mat)
    if scaling==True:
        x_renormalize_mat=x_scaling(N_renormalize,normalized,x_scale,x_shift,x_time_mat,x_mat)
    beta_mat=first_phase_beta(K,N_renormalize)
    beta_mat=beta_mat.reshape(K,K)
    if x_int=='ON':
        x_renormalize_mat=x_integer(N_renormalize,x_time_mat,x_renormalize_mat,x_int_term)
    (x_renormalize_mat,x_mean_vec,x_sd)=x_norm(N,x_time_mat,x_renormalize_mat)
    #(s_mat,lamda_vec,pi_vec,A_mat)=first_phase(CONTINUE,calc,N_renormalize,K,x_renormalize_mat,)
    #(s_mat,mu_vec,lamda_vec,pi_vec,A_mat)=first_phase(CONTINUE,calc,N_renormalize,K,x_time_mat,x_renormalize_mat,term,first_moving)
    new=1.0
    m=np.mean(x_renormalize_mat)
    a=1.0
    b=np.var(x_renormalize_mat)
    print('# (new,m)=',(new,m))
    print('# (a,b)=',(a,b))
    (mu_vec,lamda_vec)=generate_init_mu_lamda(K,new,m,a,b)
    pi_vec=generate_init_pi(K,alpha_vec)
    A_mat=first_phase_A(K,beta_mat)
    print('# mu_vec=',mu_vec)
    print('# lamda_vec=',lamda_vec)
    print('# pi_vec=',pi_vec)
    print('# (A_mat[0][0],A_mat[0][1],A_mat[1][0],A_mat[1][1])=',(A_mat[0][0],A_mat[0][1],A_mat[1][0],A_mat[1][1]))
    iter_mat_tilde=0.50*np.ones((K,N_renormalize))
    #if calc=='Comp':
       #s_mat=first_phase_s_mean(K,N_renormalize)
    #else:
       #s_mat=first_phase_s(K,N_renormalize)
    #s_mat=first_phase_s(K,N_renormalize)
    sum_s_mat=np.zeros((K,N_renormalize))

#s_mat_tilde=s_mat
#lamda_vec_tilde=lamda_vec
#pi_vec_tilde=pi_vec
#A_mat_tilde=A_mat
 
print('# (beta_mat[0][0],beta_mat[0][1],beta_mat[1][0],beta_mat[1][1])=',(beta_mat[0][0],beta_mat[0][1],beta_mat[1][0],beta_mat[1][1]))

if s_pri=='ON':
    print('# x_mean=',x_mean)
    print('# sigma=',sigma)
    sys.exit()

#sys.exit() 

lamda_log_vec=np.log(lamda_vec)
#print '# lamda_log_vec=',lamda_log_vec
pi_log_vec=np.log(pi_vec)
A_log_mat=np.log(A_mat)
if calc=='Comp':
    (mu_lamda_vec,mu_pow_lamda_vec)=first_phase_mu_lamda(K,s_mat,mu_vec,lamda_vec,new)

#s_mat=first_phase_s(K,N)
#lamda_vec=first_phase_lamda(K,a,b)
#lamda_vec=np.ones(K)
#lamda_vec=np.zeros(K)
#lamda_vec[0]=np.min(x_mat)
#lamda_vec[1]=np.max(x_mat)
#lamda_log_vec=np.log(lamda_vec)
#print '#lamda_log_vec=',lamda_log_vec
#lamda_vec=lamda_vec.reshape(K,1)
#pi_vec=first_phase_pi(K)
#pi_log_vec=np.log(pi_vec)
#A_mat=first_phase_A(K,beta_mat)
#A_mat=0.50*np.ones((K,K))
#A_mat=np.array([[0.999,0.001],[0.001,0.999]]) 
#A_log_mat=np.log(A_mat)
#sys.exit()

mu_mat=np.zeros((K,MAXITER+1))
lamda_mat=np.zeros((K,MAXITER+1))
lamda_mat=lamda_mat.reshape(K,MAXITER+1)
pi_mat=np.zeros((K,MAXITER+1))
pi_mat=pi_mat.reshape(K,MAXITER+1)
A_3d_mat=np.zeros((K,K,MAXITER+1))
A_3d_mat=A_3d_mat.reshape(K,K,MAXITER+1)

mu_mat[:,0]=mu_vec
lamda_mat[:,0]=lamda_vec
pi_mat[:,0]=pi_vec
A_3d_mat[:,:,0]=A_mat

#print's_mat=',s_mat
#print'lamda_vec=',lamda_vec
#sys.exit()
#print'A_mat=',A_mat
#print'log(A_mat)=',np.log(A_mat)
#sys.exit()

print('# sampling_term=',sampling_term)

for itr in range(0,MAXITER):
    if s==0:
        #print('# itr=',itr)
        if calc=='Gibbs':
            (a_mat,d_mat,iter_mat,s_mat)=gibbs_sampling_s(K,N,x_mat,s_mat,mu_vec,lamda_vec,pi_vec,A_mat)
            #if mixing=='ON':
                #(a_mat,d_mat,iter_mat,s_mat)=gibbs_sampling_s_simple_mixing(N,D,K,x_mat,s_mat,lamda_vec,pi_vec,A_mat,alpha,iter_mat_tilde)
                #iter_mat_tilde=iter_mat
            #else:
                #(a_mat,d_mat,iter_mat,s_mat)=gibbs_sampling_s(N,D,K,x_mat,s_mat,lamda_vec,pi_vec,A_mat)
            #if itr%sampling_term == 0:
                #wrt_iter(itr,N,K,calc,iter_mat,a_mat,d_mat)     
            #s_mat=alpha*s_mat+(1.0-alpha)*s_mat_tilde 
            #s_mat_tilde=s_mat
            (mu_vec,lamda_vec)=gibbs_sampling_mu_lamda(K,x_mat,s_mat,new,m,a,b)
            #lamda_vec=alpha*lamda_vec+(1.0-alpha)*lamda_vec_tilde 
            #lamda_vec_tilde=lamda_vec
            mu_mat[:,itr+1]=mu_vec
            #wrt_mu(K,MAXITER,calc,mu_mat)
            lamda_mat[:,itr+1]=lamda_vec
            #wrt_lamda(K,MAXITER,calc,lamda_mat)
            pi_vec=gibbs_sampling_pi(s_mat,alpha_vec)
            #pi_vec=alpha*pi_vec+(1.0-alpha)*pi_vec_tilde 
            #pi_vec_tilde=pi_vec
            pi_mat[:,itr+1]=pi_vec
            #wrt_pi(K,MAXITER,calc,pi_mat)
            A_mat=gibbs_sampling_A(K,N,s_mat,beta_mat)
            #A_mat=alpha*A_mat+(1.0-alpha)*A_mat_tilde 
            #A_mat_tilde=A_mat
            A_3d_mat[:,:,itr+1]=A_mat
            #wrt_A(K,MAXITER,calc,A_3d_mat)
            #sum_s_mat+=iter_mat 
            sum_s_mat+=s_mat 
        elif calc=='Block':
            (p_mat,s_mat)=blocking_gibbs_sampling_s(N,D,K,x_mat,mu_vec,lamda_vec,pi_vec,A_mat)
            (mu_vec,lamda_vec)=gibbs_sampling_mu_lamda(K,x_mat,s_mat,new,m,a,b)
            mu_mat[:,itr+1]=mu_vec
            #wrt_mu(K,MAXITER,calc,mu_mat)
            lamda_mat[:,itr+1]=lamda_vec
            #wrt_lamda(K,MAXITER,calc,lamda_mat)
            pi_vec=gibbs_sampling_pi(s_mat,alpha_vec)
            pi_mat[:,itr+1]=pi_vec
            #wrt_pi(K,MAXITER,calc,pi_mat)
            A_mat=gibbs_sampling_A(K,N,s_mat,beta_mat)
            A_3d_mat[:,:,itr+1]=A_mat
            #wrt_A(K,MAXITER,calc,A_3d_mat)
            sum_s_mat+=s_mat 
        elif calc=='Comp':
            (a_mat,d_mat,iter_mat,s_mat)=variational_inference_s(K,N,x_mat,s_mat,lamda_vec,lamda_log_vec,pi_log_vec,A_log_mat,mu_lamda_vec,mu_pow_lamda_vec)
            if itr%sampling_term == 0:
                wrt_iter(itr,N,K,calc,iter_mat,a_mat,d_mat)     
            (mu_vec,lamda_vec,lamda_log_vec,mu_lamda_vec,mu_pow_lamda_vec)=variational_inference_mu_lamda(K,x_mat,s_mat,new,m,a,b)
            mu_mat[:,itr+1]=mu_vec
            #wrt_mu(K,MAXITER,calc,mu_mat)
            lamda_mat[:,itr+1]=lamda_vec
            #wrt_lamda(K,MAXITER,calc,lamda_mat)
            (pi_vec,pi_log_vec)=variational_inference_pi(K,s_mat,alpha_vec)
            pi_mat[:,itr+1]=pi_vec
            #wrt_pi(K,MAXITER,calc,pi_mat)
            (A_mat,A_log_mat)=variational_inference_A(K,N,s_mat,beta_mat)
            A_3d_mat[:,:,itr+1]=A_mat
            #wrt_A(K,MAXITER,calc,A_3d_mat)
            sum_s_mat+=iter_mat 
 
    else:
        #print'# itr=',itr
        #lamda_vec=gibbs_sampling_lamda(N,D,K,x_mat,s_mat,lamda_vec,a,b)
        if calc=='Gibbs':
            (a_mat,d_mat,iter_mat,s_mat)=gibbs_sampling_s(K,N_renormalize,x_renormalize_mat,s_mat,mu_vec,lamda_vec,pi_vec,A_mat)
            #if mixing=='ON':
                #(a_mat,d_mat,iter_mat,s_mat)=gibbs_sampling_s_simple_mixing(N_renormalize,D,K,x_renormalize_mat,s_mat,lamda_vec,pi_vec,A_mat,alpha,iter_mat_tilde)
                #iter_mat_tilde=iter_mat 
            #else:
                #(a_mat,d_mat,iter_mat,s_mat)=gibbs_sampling_s(N_renormalize,D,K,x_renormalize_mat,s_mat,lamda_vec,pi_vec,A_mat)
            if itr%sampling_term == 0:
                wrt_iter(itr,N_renormalize,K,calc,iter_mat,a_mat,d_mat)     
            #s_mat=alpha*s_mat+(1.0-alpha)*s_mat_tilde 
            #s_mat_tilde=s_mat
            (mu_vec,lamda_vec)=gibbs_sampling_mu_lamda(K,x_renormalize_mat,s_mat,new,m,a,b)
            #lamda_vec=alpha*lamda_vec+(1.0-alpha)*lamda_vec_tilde 
            #lamda_vec_tilde=lamda_vec
            mu_mat[:,itr+1]=mu_vec
            #wrt_mu(K,MAXITER,calc,mu_mat)
            lamda_mat[:,itr+1]=lamda_vec
            #wrt_lamda(K,MAXITER,calc,lamda_mat)
            pi_vec=gibbs_sampling_pi(s_mat,alpha_vec)
            #pi_vec=alpha*pi_vec+(1.0-alpha)*pi_vec_tilde 
            #pi_vec_tilde=pi_vec
            pi_mat[:,itr+1]=pi_vec
            #wrt_pi(K,MAXITER,calc,pi_mat)
            A_mat=gibbs_sampling_A(K,N_renormalize,s_mat,beta_mat)
            #A_mat=alpha*A_mat+(1.0-alpha)*A_mat_tilde 
            #A_mat_tilde=A_mat
            A_3d_mat[:,:,itr+1]=A_mat
            #wrt_A(K,MAXITER,calc,A_3d_mat)
            sum_s_mat+=s_mat 
        elif calc=='Block':
            (p_mat,s_mat)=blocking_gibbs_sampling_s(N_renormalize,D,K,x_renormalize_mat,mu_vec,lamda_vec,pi_vec,A_mat)
            (mu_vec,lamda_vec)=gibbs_sampling_mu_lamda(K,x_renormalize_mat,s_mat,new,m,a,b)
            mu_mat[:,itr+1]=mu_vec
            #wrt_mu(K,MAXITER,calc,mu_mat)
            lamda_mat[:,itr+1]=lamda_vec
            #wrt_lamda(K,MAXITER,calc,lamda_mat)
            pi_vec=gibbs_sampling_pi(s_mat,alpha_vec)
            pi_mat[:,itr+1]=pi_vec
            #wrt_pi(K,MAXITER,calc,pi_mat)
            A_mat=gibbs_sampling_A(K,N_renormalize,s_mat,beta_mat)
            A_3d_mat[:,:,itr+1]=A_mat
            #wrt_A(K,MAXITER,calc,A_3d_mat)
            sum_s_mat+=s_mat 
        elif calc=='Comp':
            (a_mat,d_mat,iter_mat,s_mat)=variational_inference_s(K,N_renormalize,x_renormalize_mat,s_mat,lamda_vec,lamda_log_vec,pi_log_vec,A_log_mat,mu_lamda_vec,mu_pow_lamda_vec)
            if itr%sampling_term == 0:
                wrt_iter(itr,N_renormalize,K,calc,iter_mat,a_mat,d_mat)     
            (mu_vec,lamda_vec,lamda_log_vec,mu_lamda_vec,mu_pow_lamda_vec)=variational_inference_mu_lamda(K,x_renormalize_mat,s_mat,new,m,a,b)
            mu_mat[:,itr+1]=mu_vec
            wrt_mu(K,MAXITER,calc,mu_mat)
            lamda_mat[:,itr+1]=lamda_vec
            wrt_lamda(K,MAXITER,calc,lamda_mat)
            (pi_vec,pi_log_vec)=variational_inference_pi(K,s_mat,alpha_vec)
            pi_mat[:,itr+1]=pi_vec
            wrt_pi(K,MAXITER,calc,pi_mat)
            (A_mat,A_log_mat)=variational_inference_A(K,N_renormalize,s_mat,beta_mat)
            A_3d_mat[:,:,itr+1]=A_mat
            wrt_A(K,MAXITER,calc,A_3d_mat)
            sum_s_mat+=iter_mat 

    #output 
    if itr%sampling_term == 0:
        print('# itr=',itr)
    #if itr%10 == 0:
    #if itr%1 == 0:
    #if itr%100 == 0:
        #if itr==0: 
            #continue 
        s_mat_ave=sum_s_mat/(itr+1)
        if s==0:
            #wrt_s_ave_itr(itr,N,K,calc,x_time_mat,x_mat,s_mat,iter_mat)
            wrt_s_ave_itr(itr,N,K,calc,x_time_mat,x_mat,s_mat,s_mat_ave)
        else:
            #wrt_s_ave_itr(itr,N_renormalize,K,calc,x_renormalize_time_mat,x_renormalize_mat,s_mat,iter_mat)
            wrt_s_ave_itr(itr,N_renormalize,K,calc,x_renormalize_time_mat,x_renormalize_mat,s_mat,s_mat_ave)
        t2=time.time() 
        elapsed_time=t2-t1
        print('# elapsed_time=',elapsed_time)
  

#print'A_vec=',A_vec
#print'A_mat=',A_mat
#print'lamda_vec=',lamda_vec
#print's_mat=',s_mat

s_mat_ave=sum_s_mat/MAXITER 
#print's_mat_ave=',s_mat_ave
if s==0:
    #wrt_s_ave(N,K,calc,x_time_mat,x_mat,s_mat,iter_mat)
    wrt_s_ave(N,K,calc,x_time_mat,x_mat,s_mat,s_mat_ave)
    wrt_s_continue(N,K,calc,s_mat_ave)
    wrt_x_restored_mat(N,x_time_mat,x_mat,x_mean_vec,x_sd)  
else:
    #wrt_s_ave(N_renormalize,K,calc,x_renormalize_time_mat,x_renormalize_mat,s_mat,iter_mat)
    wrt_s_ave(N_renormalize,K,calc,x_renormalize_time_mat,x_renormalize_mat,s_mat,s_mat_ave)
    wrt_s_continue(N_renormalize,K,calc,s_mat_ave)
    wrt_x_restored_mat(N_renormalize,x_time_mat,x_renormalize_mat,x_mean_vec,x_sd)  

wrt_mu_continue(K,calc,mu_vec)
wrt_lamda_continue(K,calc,lamda_vec)
wrt_pi_continue(K,calc,pi_vec)
wrt_A_continue(K,calc,A_mat)
#wrt_lamda(K,MAXITER,calc,lamda_mat)
#wrt_pi(K,MAXITER,calc,pi_mat)
#wrt_A(K,MAXITER,calc,A_3d_mat)
#wrt_s_ave(N,K,s_mat)
t2=time.time() 
elapsed_time=t2-t1
print('# elapsed_time=',elapsed_time)
