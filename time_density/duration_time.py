#coding: utf-8
#(c) Tatsuhiro Furuta
import sys,os
import numpy as np

def A_function(N,z_mat,t_mat):
     A_mat=np.zeros(N-1)
     A_mat=A_mat.reshape(1,N-1)
     a=0
     for n in range(N-1):
        a=(z_mat[0][n+1]-z_mat[0][n])/(t_mat[0][n+1]-t_mat[0][n])
        A_mat[0][n]=a
        #print a
     #print 'A_mat=',A_mat
     #sys.exit()
     return A_mat

def A_mean_function(A_mat):
     a=1
     b=0
     A=0
     for n in range(A_mat.shape[1]):
        if np.abs(A_mat[0][n]) >= a:
           A=A+np.abs(A_mat[0][n])
           b=b+1
     A_mean=A/b
     return A_mean

def duration_time(N,x,z_mat,t_mat):
     duration_mat=np.zeros(N)
     duration_mat=duration_mat.reshape(1,N)

     #x=20

     #event_count
     v=0; w=0; m=0
     for n in range(1,N):
         if z_mat[0][m]<x and z_mat[0][n]>=x:
             v=v+1
             m=n
         elif z_mat[0][m]>=x and z_mat[0][n]<x:
             w=w+1
             m=n

     if m==N-1:
         pass
     else:
         if z_mat[0][m]<x:
             v=v+1
         elif z_mat[0][m]>=x:
             w=w+1
 
     ON_mat=np.zeros(w)
     OFF_mat=np.zeros(v)
     ON_mat=ON_mat.reshape(1,w)
     OFF_mat=OFF_mat.reshape(1,v)
     ON_ini_mat=np.zeros(w)
     ON_end_mat=np.zeros(w)
     OFF_ini_mat=np.zeros(v)
     OFF_end_mat=np.zeros(v)
     #print("# w=",w)
     #print("# v=",v)

     #duration_time
     m=0; d=0; e=0;w_i=0;w_e=0;v_i=0;v_e=0

     if z_mat[0][m]<x:
         OFF_ini_mat[v_i]=t_mat[0][m]
         v_i=v_i+1
     else:
         ON_ini_mat[w_i]=t_mat[0][m]
         w_i=w_i+1

     for n in range(1,N):
         if z_mat[0][m]<x and z_mat[0][n]>=x:
             OFF_mat[0][d]=t_mat[0][n]-t_mat[0][m]
             for i in range(m,n):
                 duration_mat[0][i]=0.0
                 #print 'OFF_mat[0]['+str(d)+']=',OFF_mat[0][d]
             m=n
             d=d+1
             OFF_end_mat[v_e]=t_mat[0][n]
             v_e=v_e+1
             if n==N-1:
                 pass
             else:
                 ON_ini_mat[w_i]=t_mat[0][n]
                 w_i=w_i+1
             #print("# w_i=",w_i)
             #print("# ON_ini_mat[w_i]=",ON_ini_mat[w_i])
         elif z_mat[0][m]>=x and z_mat[0][n]<x:
             ON_mat[0][e]=t_mat[0][n]-t_mat[0][m]
             for i in range(m,n):
                 duration_mat[0][i]=1.0
                 #print 'ON_mat[0]['+str(e)+']=',ON_mat[0][e]
             #print m
             m=n
             e=e+1
             ON_end_mat[w_e]=t_mat[0][n]
             w_e=w_e+1
             if n==N-1:
                 pass
             else:
                 OFF_ini_mat[v_i]=t_mat[0][n]
                 v_i=v_i+1
             #print("# v_i=",v_i)
             #print("# OFF_ini_mat[v_i]=",OFF_ini_mat[v_i])
    
     if m==N-1:
         pass
     else:
         if z_mat[0][m]<x:
             OFF_mat[0][d]=t_mat[0][N-1]-t_mat[0][m]
             OFF_end_mat[v_e]=t_mat[0][N-1]
             for i in range(m,n):
                 duration_mat[0][i]=0.0
         elif z_mat[0][m]>=x:
             ON_mat[0][e]=t_mat[0][N-1]-t_mat[0][m]
             ON_end_mat[w_e]=t_mat[0][N-1]
             for i in range(m,n):
                 duration_mat[0][i]=1.0
         
     #return ON_mat,OFF_mat,duration_mat
     return ON_mat,OFF_mat,duration_mat,ON_ini_mat,ON_end_mat,OFF_ini_mat,OFF_end_mat
'''
     v=0; w=0; m=0
     for n in range(N):
        if n==0:
           if z_mat[0][0]>=x:
              if z_mat[0][1]<=x and z_mat[0][2]>x:
                 v=v+1
                 m=1
           elif z_mat[0][0]<=x:
              if z_mat[0][1]>=x and z_mat[0][2]<x:
                 w=w+1
                 m=1
        elif n==N-1:
           if z_mat[0][m]>=x:
              v=v+1
           elif z_mat[0][m]<=x:
              w=w+1 
        else:   
           if z_mat[0][n]<=x:
              if z_mat[0][n+1]>x: 
                 v=v+1
                 m=n
           elif z_mat[0][n]>=x:
              if z_mat[0][n+1]<x: 
                 w=w+1
                 m=n

     #print'v=',v
     #print'w=',w
     #print'm=',m
     #print'b=',b
     #print'c=',c
     #print'f=',f
     #print'g=',g
     #print'm=',m
     #print'A_mat[0][m]=',A_mat[0][m]


     #ON=ON+v
     #ON=ON+w
     #OFF=OFF+w
     #OFF=OFF+v
     #print 'ON=',ON
     #print 'OFF=',OFF

     ON_mat=np.zeros(w)
     OFF_mat=np.zeros(v)
     ON_mat=ON_mat.reshape(1,w)
     OFF_mat=OFF_mat.reshape(1,v)

     m=0; d=0; e=0
     for n in range(N):
        if n==0:
           if z_mat[0][0]>=x:
              if z_mat[0][1]<=x and z_mat[0][2]>x:
                 OFF_mat[0][0]=t_mat[0][1]-t_mat[0][0]
                 duration_mat[0][1]=0.0
                 #print 'ON_mat[0]['+str(e)+']=',ON_mat[0][e]
                 m=1
                 d=d+1
           elif z_mat[0][0]<=x:
              if z_mat[0][1]>=x and z_mat[0][2]<x:
                 ON_mat[0][0]=t_mat[0][1]-t_mat[0][0]
                 duration_mat[0][1]=1.0
                 #print 'ON_mat[0]['+str(e)+']=',ON_mat[0][e]
                 m=1
                 e=e+1
        elif n==N-1:
           if z_mat[0][m]>=x:
              OFF_mat[0][d]=t_mat[0][N-1]-t_mat[0][m]
              for i in range(m,N):
                 duration_mat[0][i]=0.0
              #print 'OFF_mat[0]['+str(d)+']=',OFF_mat[0][d]
           elif z_mat[0][m]<=x:
              ON_mat[0][e]=t_mat[0][N-1]-t_mat[0][m]    
              for i in range(m,N):
                 duration_mat[0][i]=1.0
              #print 'ON_mat[0]['+str(e)+']=',ON_mat[0][e]
        else:   
           if z_mat[0][n]<=x:
              if z_mat[0][n+1]>x: 
                 OFF_mat[0][d]=t_mat[0][n]-t_mat[0][m]
                 for i in range(m,n):
                    duration_mat[0][i]=0.0
                 #print 'OFF_mat[0]['+str(d)+']=',OFF_mat[0][d]
                 m=n
                 d=d+1
           elif z_mat[0][n]>=x:
              if z_mat[0][n+1]<x: 
                 ON_mat[0][e]=t_mat[0][n]-t_mat[0][m]
                 for i in range(m,n):
                    duration_mat[0][i]=1.0
                 #print 'ON_mat[0]['+str(e)+']=',ON_mat[0][e]
                 m=n
                 e=e+1
'''     

     #print'd=',d
     #print'e=',e
     #print'm=',m
     #print'ONtot_mat=',ON_mat

     #print'OFFtot_mat=',OFF_mat

     #sum_ON=0
     #sum_OFF=0

     #for i in range(ON):
        #sum_ON=sum_ON+ON_mat[0][i]
     #print'sum_ON=',sum_ON

     #for j in range(OFF):
        #sum_OFF=sum_OFF+OFF_mat[0][j]
     #print'sum_OFF=',sum_OFF

     #return ON_mat,OFF_mat,duration_mat
          
def duration_time_kai(N,ON_line,OFF_line,z_mat,t_mat):
     duration_mat=np.zeros(N)
     duration_mat=duration_mat.reshape(1,N)

     #x=20
     v=0; w=0; m=0
     for n in range(1,N):
         if m==0:
             if OFF_line<z_mat[0][m]<ON_line:
                 if z_mat[0][n]<=OFF_line or z_mat[0][n]>=ON_line:
                     m=n
             else:
                 if n==N-1:
                     if z_mat[0][m]<=OFF_line:
                         v=v+1
                     elif z_mat[0][m]>=ON_line:
                         w=w+1
                 else:
                     if z_mat[0][m]<=OFF_line and z_mat[0][n]>=ON_line:
                         v=v+1
                         m=n
                     elif z_mat[0][m]>=ON_line and z_mat[0][n]<=OFF_line:
                         w=w+1
                         m=n
         else:
             if n==N-1:
                 if z_mat[0][m]<=OFF_line:
                     v=v+1
                 elif z_mat[0][m]>=ON_line:
                     w=w+1
             else:
                 if z_mat[0][m]<=OFF_line and z_mat[0][n]>=ON_line:
                     v=v+1
                     m=n
                 elif z_mat[0][m]>=ON_line and z_mat[0][n]<=OFF_line:
                     w=w+1
                     m=n

     #print'v=',v
     #print'w=',w
     #print'm=',m
     #print'b=',b
     #print'c=',c
     #print'f=',f
     #print'g=',g
     #print'm=',m
     #print'A_mat[0][m]=',A_mat[0][m]


     #ON=ON+v
     #ON=ON+w
     #OFF=OFF+w
     #OFF=OFF+v
     #print 'ON=',ON
     #print 'OFF=',OFF

     ON_mat=np.zeros(w)
     OFF_mat=np.zeros(v)
     ON_mat=ON_mat.reshape(1,w)
     OFF_mat=OFF_mat.reshape(1,v)
     ON_ini_mat=np.zeros(w)
     ON_end_mat=np.zeros(w)
     OFF_ini_mat=np.zeros(v)
     OFF_end_mat=np.zeros(v)

     m=0; d=0; e=0;w_i=0;w_e=0;v_i=0;v_e=0
     for n in range(1,N):
         #print 'n=',n
         #print 'm=',m
         if m==0:
             if OFF_line<z_mat[0][m]<ON_line:
                 if z_mat[0][n]<=OFF_line or z_mat[0][n]>=ON_line:
                     m=n
                     if z_mat[0][n]<=OFF_line:
                         OFF_ini_mat[v_i]=t_mat[0][n]
                         v_i=v_i+1
                     elif z_mat[0][n]>=ON_line:
                         ON_ini_mat[w_i]=t_mat[0][n]
                         w_i=w_i+1
             else:
                 if n==N-1:
                     if z_mat[0][m]<=OFF_line:
                         OFF_mat[0][d]=t_mat[0][n]-t_mat[0][m]
                         for i in range(m,n+1):
                             duration_mat[0][i]=0.0
                         OFF_ini_mat[v_i]=t_mat[0][m]
                         OFF_end_mat[v_e]=t_mat[0][n]
                     elif z_mat[0][m]>=ON_line:
                         ON_mat[0][e]=t_mat[0][n]-t_mat[0][m]
                         for i in range(m,n+1):
                             duration_mat[0][i]=1.0
                         ON_ini_mat[w_i]=t_mat[0][m]
                         ON_end_mat[w_e]=t_mat[0][n]
                 else:
                     if z_mat[0][m]<=OFF_line and z_mat[0][n]>=ON_line:
                         OFF_mat[0][d]=t_mat[0][n]-t_mat[0][m]
                         for i in range(m,n+1):
                             duration_mat[0][i]=0.0
                             #print 'OFF_mat[0]['+str(d)+']=',OFF_mat[0][d]
                         OFF_ini_mat[v_i]=t_mat[0][m]
                         v_i=v_i+1
                         m=n
                         d=d+1
                         ON_ini_mat[w_i]=t_mat[0][n]
                         OFF_end_mat[v_e]=t_mat[0][n]
                         w_i=w_i+1
                         v_e=v_e+1
                     elif z_mat[0][m]>=ON_line and z_mat[0][n]<=OFF_line:
                         ON_mat[0][e]=t_mat[0][n]-t_mat[0][m]
                         for i in range(m,n+1):
                             duration_mat[0][i]=1.0
                             #print 'ON_mat[0]['+str(e)+']=',ON_mat[0][e]
                         #print m
                         ON_ini_mat[w_i]=t_mat[0][m]
                         w_i=w_i+1
                         m=n
                         e=e+1
                         ON_end_mat[w_e]=t_mat[0][n]
                         OFF_ini_mat[v_i]=t_mat[0][n]
                         w_e=w_e+1
                         v_i=v_i+1
         else:
             if n==N-1:
                 if z_mat[0][m]<=OFF_line:
                     OFF_mat[0][d]=t_mat[0][n]-t_mat[0][m]
                     for i in range(m,n+1):
                         duration_mat[0][i]=0.0
                     OFF_end_mat[v_e]=t_mat[0][n]
                 elif z_mat[0][m]>=ON_line:
                     ON_mat[0][e]=t_mat[0][n]-t_mat[0][m]
                     for i in range(m,n+1):
                         duration_mat[0][i]=1.0
                     ON_end_mat[w_e]=t_mat[0][n]
             else:
                 if z_mat[0][m]<=OFF_line and z_mat[0][n]>=ON_line:
                     OFF_mat[0][d]=t_mat[0][n]-t_mat[0][m]
                     for i in range(m,n+1):
                         duration_mat[0][i]=0.0
                         #print 'OFF_mat[0]['+str(d)+']=',OFF_mat[0][d]
                     m=n
                     d=d+1
                     ON_ini_mat[w_i]=t_mat[0][n]
                     OFF_end_mat[v_e]=t_mat[0][n]
                     w_i=w_i+1
                     v_e=v_e+1
                 elif z_mat[0][m]>=ON_line and z_mat[0][n]<=OFF_line:
                     ON_mat[0][e]=t_mat[0][n]-t_mat[0][m]
                     for i in range(m,n+1):
                         duration_mat[0][i]=1.0
                         #print 'ON_mat[0]['+str(e)+']=',ON_mat[0][e]
                     #print m
                     m=n
                     e=e+1
                     ON_end_mat[w_e]=t_mat[0][n]
                     OFF_ini_mat[v_i]=t_mat[0][n]
                     w_e=w_e+1
                     v_i=v_i+1
     
     #print'd=',d
     #print'e=',e
     #print'm=',m
     #print'ONtot_mat=',ON_mat

     #print'OFFtot_mat=',OFF_mat

     #sum_ON=0
     #sum_OFF=0

     #for i in range(ON):
        #sum_ON=sum_ON+ON_mat[0][i]
     #print'sum_ON=',sum_ON

     #for j in range(OFF):
        #sum_OFF=sum_OFF+OFF_mat[0][j]
     #print'sum_OFF=',sum_OFF

     return ON_mat,OFF_mat,duration_mat,ON_ini_mat,ON_end_mat,OFF_ini_mat,OFF_end_mat
