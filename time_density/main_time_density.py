#coding: utf-8
#(c) Kazuma Nakamura,Tatsuhiro Furuta
import sys,os
import numpy as np
import argparse
from distutils.util import strtobool
from density import *

parser=argparse.ArgumentParser()
parser.add_argument('--few_photon_data',type=strtobool,default=False,help='If few_photon_data==True,distribution_density_for_few_photon_data function works instead of distribution_density function.(default: False)') 
parser.add_argument('--a',type=float,default=0.01,help='Grid number of duration time(default: 0.01)') 

args = parser.parse_args()
a=args.a
few_photon_data=args.few_photon_data

print('# a=',a)
print('# few_photon_data=',bool(few_photon_data))

f=open('ON_mat.csv','r')
line=f.readlines()
f.close()

tmp=np.loadtxt(line)
#print tmp.shape 
#sys.exit() 
     
Nd=tmp.shape[0] 
N=1
print('# Nd=',Nd)
print('# N=',N)
x=np.zeros((N,Nd))
for n in range(Nd): 
    x[0][n]=float(tmp[n]) 
ON_mat=x

if few_photon_data==True:
    #ON_dens_mat=distribution_density_for_few_photon_data(a,ON_mat)
    ON_dens_mat=distribution_density_weighted(a,ON_mat)
else:
    ON_dens_mat=distribution_density(a,ON_mat)

f=open('ON_distribution_density.csv','w')
for j in range(ON_dens_mat.shape[1]):
    f.write(str(ON_dens_mat[0,j]))
    f.write('\t')
    f.write(str(ON_dens_mat[1,j]))
    f.write('\n')
f.close()

f=open('OFF_mat.csv','r')
line=f.readlines()
f.close()

tmp=np.loadtxt(line)
#print tmp.shape 
#sys.exit() 
     
Nd=tmp.shape[0] 
N=1
print('# Nd=',Nd)
print('# N=',N)
y=np.zeros((N,Nd))
for n in range(Nd): 
    y[0][n]=float(tmp[n]) 
OFF_mat=y

if few_photon_data==True:
    #OFF_dens_mat=distribution_density_for_few_photon_data(a,OFF_mat)
    OFF_dens_mat=distribution_density_weighted(a,OFF_mat)
else:
    OFF_dens_mat=distribution_density(a,OFF_mat)

f=open('OFF_distribution_density.csv','w')
for j in range(OFF_dens_mat.shape[1]):
    f.write(str(OFF_dens_mat[0,j]))
    f.write('\t')
    f.write(str(OFF_dens_mat[1,j]))
    f.write('\n')
f.close()
