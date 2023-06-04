#coding: utf-8
import sys,os
import numpy as np
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--s_mu',type=float,default=0.01,help='Permission range of difference from mean(default: 0.01)') 
parser.add_argument('--sigma',type=float,default=0.01,help='Permission range of standard deviation(default: 0.01)') 
args = parser.parse_args()

s_mu=args.s_mu
sigma=args.sigma
print'# (s_mu,sigma)=',(s_mu,sigma)
#sys.exit()

f=open('s_mat-ON.csv','r') 
line=f.readlines()#[1:] 
f.close()
tmp=np.loadtxt(line)
#print tmp.shape 
#sys.exit() 
#
N=tmp.shape[0] 
K=tmp.shape[1] 
#sys.exit() 
#
s_mat=np.zeros((K,N))
for n in range(N):
    for k in range(K):
        s_mat[k][n]=float(tmp[n][k])
#print s_mat.shape

f=open('T_ON.txt','r') 
line=f.readlines()#[1:] 
f.close()
tmp=np.loadtxt(line)
#print tmp.shape 
#sys.exit() 
#
event1=tmp.shape[0] 
K1=tmp.shape[1] 
#sys.exit() 
#
ON=np.zeros((K1,event1))
for n in range(event1):
    for k in range(K1):
        ON[k][n]=float(tmp[n][k])

f=open('T_OFF.txt','r') 
line=f.readlines()#[1:] 
f.close()
tmp=np.loadtxt(line)
#print tmp.shape 
#sys.exit() 
#
event2=tmp.shape[0] 
K2=tmp.shape[1] 
#sys.exit() 
#
OFF=np.zeros((K2,event2))
for n in range(event2):
    for k in range(K2):
        OFF[k][n]=float(tmp[n][k])

event=event1

s_ON_mean=np.zeros((1,event))
s_ON_mean=s_ON_mean.reshape(1,event)

s_ON_sig=np.zeros((1,event))
s_ON_sig=s_ON_sig.reshape(1,event)

s_OFF_mean=np.zeros((1,event))
s_OFF_mean=s_OFF_mean.reshape(1,event)

s_OFF_sig=np.zeros((1,event))
s_OFF_sig=s_OFF_sig.reshape(1,event)

#for i in range(int(ON[0][1]),int(ON[1][1])):
#    print (s_mat[0][i],s_mat[1][i])
#for i in range(int(OFF[0][1]),int(OFF[1][1])):
#    print (s_mat[0][i],s_mat[1][i])
#print s_mat[0,int(ON[0][0]):int(ON[1][0])]
#print s_mat[0,int(ON[0][1]):int(ON[1][1])]
#print s_mat[0,int(OFF[0][0]):int(OFF[1][0])]
#print s_mat[0,int(OFF[0][1]):int(OFF[1][1])]

for j in range(event):
    s_ON_mean[0][j]=np.mean(s_mat[1,int(ON[0][j]):int(ON[1][j])])
    s_ON_sig[0][j]=np.std(s_mat[1,int(ON[0][j]):int(ON[1][j])])
    #print 's_ON_mean,s_ON_sig',s_ON_mean[0][j],s_ON_sig[0][j]
    s_OFF_mean[0][j]=np.mean(s_mat[1,int(OFF[0][j]):int(OFF[1][j])])
    s_OFF_sig[0][j]=np.std(s_mat[1,int(OFF[0][j]):int(OFF[1][j])])
    #print 's_OFF_mean,s_OFF_sig',s_OFF_mean[0][j],s_OFF_sig[0][j]

result=np.array([['Success'],['Fail']])
#print result.shape

result_ON=[]
result_OFF=[]

mu=np.array([[1.0],[0.0]])
#print mu.shape
#print mu
#s_mu=0.20
#sigma=0.20

for j in range(event):
    if np.abs(s_ON_mean[0][j]-mu[0][0])<s_mu:
        if s_ON_sig[0][j]<sigma:
            result_ON.append(result[0][0])
        else:
            result_ON.append(result[1][0])
    else:
        result_ON.append(result[1][0])
    if np.abs(s_OFF_mean[0][j]-mu[1][0])<s_mu:
        if s_OFF_sig[0][j]<sigma:
            result_OFF.append(result[0][0])
        else:
            result_OFF.append(result[1][0])
    else:
        result_OFF.append(result[1][0])
#print result_ON
#print result_OFF

a='ON'
b='OFF' 
#for j in range(event):
#    print '(event,event_number,result,mean,standard_deviation)=',(a,j,result_ON[j],s_ON_mean[0][j],s_ON_sig[0][j])               
#for j in range(event):
#    print '(event,event_number,result,mean,standard_deviation)=',(b,j,result_OFF[j],s_OFF_mean[0][j],s_OFF_sig[0][j])               

f=open('result_ON.txt','w')
for j in range(event):
    f.write('(event,event_number,result,mean,standard_deviation)=('+str(a)+','+str(j)+','+str(result_ON[j])+','+str(s_ON_mean[0][j])+','+str(s_ON_sig[0][j])+')')               
    f.write('\n')
f.close()
#f.write('event')
#f.write('\t')
#f.write('event_number')
#f.write('\t')
#f.write('result')
#f.write('\t')
#f.write('mean')
#f.write('\t')
#f.write('standard_deviation')
#f.write('\n')
#for j in range(event):
#    f.write(str(a))               
#    f.write('\t')
#    f.write(str(j))               
#    f.write('\t')
#    f.write(str(result_ON[j]))               
#    f.write('\t')
#    f.write(str(s_ON_mean[0][j]))               
#    f.write('\t')
#    f.write(str(s_ON_sig[0][j]))               
#    f.write('\n')
#f.close()

f=open('result_OFF.txt','w')
for j in range(event):
    f.write('(event,event_number,result,mean,standard_deviation)=('+str(b)+','+str(j)+','+str(result_OFF[j])+','+str(s_OFF_mean[0][j])+','+str(s_OFF_sig[0][j])+')')               
    f.write('\n')
f.close()
#f.write('event')
#f.write('\t')
#f.write('event_number')
#f.write('\t')
#f.write('result')
#f.write('\t')
#f.write('mean')
#f.write('\t')
#f.write('standard_deviation')
#f.write('\n')
#for j in range(event):
#    f.write(str(b))               
#    f.write('\t')
#    f.write(str(j))               
#    f.write('\t')
#    f.write(str(result_OFF[j]))               
#    f.write('\t')
#    f.write(str(s_OFF_mean[0][j]))               
#    f.write('\t')
#    f.write(str(s_OFF_sig[0][j]))               
#    f.write('\n')
#f.close()

c=0
d=0
for j in range(event):
    if result_ON[j]=='Fail':
        c=c+1

for j in range(event):
    if result_OFF[j]=='Fail':
        d=d+1
    
#print c,d

f=open('fail.txt','w')
f.write(str(c))
f.write('\t')
f.write(str(d))
f.close()
