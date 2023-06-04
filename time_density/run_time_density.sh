#!/bin/bash

#cp -rf HMM-b 00
#cd 00
#qsub wsh.mpirun.murasaki 

#j='s_mat-ON-Gibbs.csv'
j='s_mat-ON-Block.csv'
#j='s_mat_eye_ON.csv'
k='T_no_noise.csv'
#j='s_mat-ON-itr40000-Gibbs.csv'
echo $j
echo $k
for i in {0..999}; do
    echo $i
    #cp ../$i/s_mat-ON.csv .
    cp ../$i/$j .
    #cp ../$i/eye/$j .
    #cp ../$i/$k .
    cp ../$i/$k .
    mv $j s_mat-ON-$i.csv
    mv $k T_no_noise-$i.csv
    #mv s_mat-ON.csv s_mat-ON-$i.csv
    #qsub wsh.mpirun.moegi	 
    #qsub wsh.mpirun.murasaki 
done 

#python ./main_duration.py --M 2000 --ON_line 0.9 --OFF_line 0.1 --event_error OFF --event_number 2 --T_ON_line 1.65 --T_OFF_line 1.1 --answer_check ON --T_sig 0.5 > LOG.duration
#python ./main_duration.py --M 1000 --ON_line 0.9 --OFF_line 0.1 
#python ./main_time_density.py --a 0.1 > LOG.time_density
#python ./main_time_density.py --a 0.01 > LOG.time_density
qsub wsh.mpirun.murasaki

#for i in 5 9 14 18 19 28 31 34 39 43; do
#    echo $i
#    cp -rf HMM-b $i
#    cd ${i}
#    pwd 
#    qsub wsh.mpirun.murasaki 
#    cd ..
#done 
