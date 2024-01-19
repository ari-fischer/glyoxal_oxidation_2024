#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import scipy as sci
from matplotlib import pyplot as plt
import os

#import cartesian coordinates from QChem input or output file in format:
#  "C        -0.7331289090   -0.6760306448    0.3672980418"
name_file = r'geom_out'
df = pd.read_csv(name_file+'.txt',skiprows=0,header=None,delim_whitespace=True)
df=df.rename(columns={0: 'atom',1: 'x',2:'y',3:'z'})

# find the mean x y z positions for the structure
COP = [df['x'].mean(),df['y'].mean(),df['z'].mean()]

#split into matrices consisting of coordinates of each element
df_H = df[df['atom']=='H']
df_C = df[df['atom']=='C']
df_O = df[df['atom']=='O']

#convert pandas df to numpy matrix and shift origin to the mean position
M_H = df_H[['x','y','z']].to_numpy()-COP
M_C = df_C[['x','y','z']].to_numpy()-COP
M_O = df_O[['x','y','z']].to_numpy()-COP

#calculate vdW V for a C-atom
vdW_R_C = 1.7 #radius, Ang
vdW_V_C = np.pi*4/3*vdW_R_C**3 #vol of single C
#H-atom
vdW_R_H = 1.2 #radius, Ang
vdW_V_H = np.pi*4/3*vdW_R_H**3 #vol of single H
#O-atom
vdW_R_O = 1.52 
vdW_V_O = np.pi*4/3*vdW_R_O**3 #vol of single C

##construct a 3D grid to enclose the molecule

#find the limits to be +- 3 angstroms beyond the min and max coordinates
c_mins = np.array([np.min(df['x']),np.min(df['y']),np.min(df['z'])])-3
c_maxs = np.array([np.max(df['x']),np.max(df['y']),np.max(df['z'])])+3
#gives the range of cartesian coordinates
c_ranges = c_maxs-c_mins

#dimensions of dividing cubes
d = 0.01 #discrete size
#number of elements in each dimension
M_size = np.ceil(c_ranges/d)
#make them all even
M_size = (np.ceil(M_size/2))*2 

#construct a matrix of zeroes to represent the grid with dimesions from M_size
M_fill = np.zeros([int(M_size[0]),int(M_size[1]),int(M_size[2])])

#find the grid coordinates of all H, C, and O atoms
pos_H = M_H/d+np.ceil(M_size/2)
pos_C = M_C/d+np.ceil(M_size/2)
pos_O = M_O/d+np.ceil(M_size/2)


#volume and coordinates of C
C_size = [int(2/d)*2,int(2/d)*2,int(2/d)*2]
C_fill = np.zeros(C_size)
R = vdW_R_C

zs = np.array([])

#build a volume element for the C-atom in coordinate/grid space
#loop over x
for i in range(int(np.ceil(R)/d)):
#loop over y
    z_slice = np.array([])
    for j in range(int(np.ceil(R)/d)):
        x = i
        y = j
        z_plus = +np.sqrt((R/d)**2-x**2-y**2)
        z_minus = -np.sqrt((R/d)**2-x**2-y**2)
        z_slice = np.append(z_slice,[z_plus])
    z_slice = np.append(np.flip(z_slice),z_slice)
    
    for k in range(C_size[0]):
        #print(i)
        z_max = np.ceil(z_slice[k])
        if z_max>0:
            C_fill[i+int(C_size[0]/2),k,int(C_size[0]/2-z_max-1):int(z_max-1+C_size[0]/2)]=1
            C_fill[int(C_size[0]/2)-i,k,int(C_size[0]/2-z_max-1):int(z_max-1+C_size[0]/2)]=1
            

#volume and coordinates of H
H_size = [int(2/d)*2,int(2/d)*2,int(2/d)*2]
H_fill = np.zeros(H_size)
R = vdW_R_H
zs = np.array([])

#build a volume element for the H-atom in coordinate/grid space
#loop over x
for i in range(int(np.ceil(R)/d)):
#loop over y
    z_slice = np.array([])
    for j in range(int(np.ceil(R)/d)):
        x = i
        y = j
        z_plus = +np.sqrt((R/d)**2-x**2-y**2)
        z_minus = -np.sqrt((R/d)**2-x**2-y**2)
        z_slice = np.append(z_slice,[z_plus])
    z_slice = np.append(np.flip(z_slice),z_slice)
    
    for k in range(H_size[0]):
        #print(i)
        z_max = np.ceil(z_slice[k])
        if z_max>0:
            H_fill[i+int(H_size[0]/2),k,int(H_size[0]/2-z_max-1):int(z_max-1+H_size[0]/2)]=1
            H_fill[int(H_size[0]/2)-i,k,int(H_size[0]/2-z_max-1):int(z_max-1+H_size[0]/2)]=1
    

#volume and coordinates of O
O_size = [int(2/d)*2,int(2/d)*2,int(2/d)*2]
O_fill = np.zeros(O_size)
R = vdW_R_O

zs = np.array([])

#build a volume element for the O-atom in coordinate/grid space
#loop over x
for i in range(int(np.ceil(R)/d)):
#loop over y
    z_slice = np.array([])
    for j in range(int(np.ceil(R)/d)):
        x = i
        y = j
        z_plus = +np.sqrt((R/d)**2-x**2-y**2)
        z_minus = -np.sqrt((R/d)**2-x**2-y**2)
        z_slice = np.append(z_slice,[z_plus])
    z_slice = np.append(np.flip(z_slice),z_slice)
    
    for k in range(O_size[0]):
        z_max = np.ceil(z_slice[k])
        if z_max>0:
            O_fill[i+int(O_size[0]/2),k,int(O_size[0]/2-z_max-1):int(z_max-1+O_size[0]/2)]=1
            O_fill[int(O_size[0]/2)-i,k,int(O_size[0]/2-z_max-1):int(z_max-1+O_size[0]/2)]=1
            

    
#superimpose the volume element for each element onto the molecular grid centered around their atomic coordinates
#take coordinate of each O, C, and H. Add matrix elements +/- dimensions of the O_atom to the M matrix.

i = 0
for i in range(np.shape(pos_O)[0]):
    M_fill[int(np.ceil(pos_O[i][0])-int(O_size[0]/2)):int(np.ceil(pos_O[i][0])+int(O_size[0]/2)),
           int(np.ceil(pos_O[i][1])-int(O_size[0]/2)): int(np.ceil(pos_O[i][1])+int(O_size[0]/2)),
           int(np.ceil(pos_O[i][2])-int(O_size[0]/2)): int(np.ceil(pos_O[i][2])+int(O_size[0]/2))] += O_fill

i = 0
for i in range(np.shape(pos_C)[0]):
    M_fill[int(np.ceil(pos_C[i][0])-int(C_size[0]/2)):int(np.ceil(pos_C[i][0])+int(C_size[0]/2)),
           int(np.ceil(pos_C[i][1])-int(C_size[0]/2)): int(np.ceil(pos_C[i][1])+int(C_size[0]/2)),
           int(np.ceil(pos_C[i][2])-int(C_size[0]/2)): int(np.ceil(pos_C[i][2])+int(C_size[0]/2))] += C_fill

i = 0  
for i in range(np.shape(pos_H)[0]):
    M_fill[int(np.ceil(pos_H[i][0])-int(H_size[0]/2)):int(np.ceil(pos_H[i][0])+int(H_size[0]/2)),
           int(np.ceil(pos_H[i][1])-int(H_size[0]/2)): int(np.ceil(pos_H[i][1])+int(H_size[0]/2)),
           int(np.ceil(pos_H[i][2])-int(H_size[0]/2)): int(np.ceil(pos_H[i][2])+int(H_size[0]/2))] += H_fill


#avoid double-counting by replacing values larger than 1 with 1
M_fill[M_fill>1]=1

#sum all of the ones in the volume element and multiple by d**3 to get the molecular volume.
vol =np.sum(M_fill)*d**3

s = 'volume: '+str(vol)

out = open(r'volume'+'.txt','w+')

out.write(name_file)
out.write('\n')
out.write(s)
out.write('\n')
out.write('end')

out.close()

