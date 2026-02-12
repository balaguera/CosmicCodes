#!/usr/bin/env python3
import numpy as np
import scipy.interpolate
import scipy.ndimage
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from PIL import Image
import sys
from scipy.ndimage.filters import gaussian_filter
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LogNorm
import csv
import math
import matplotlib.animation as animation
import time
import os.path
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import h5py
import pickle
# ============================================================================
# Input parameters
N=500
L=1500.
MAS = 0

MAS_NAME = ['ngp', 'cic','tsc']


if(MAS>0):
    Ncells_MAX = 3
else:
    Ncells_MAX = 1

# Derived parameters
delta = L/N
xmin=0
ymin=0
zmin=0
Ngrid = N*N*N

#THESE ARE THE DIFFERENT DATA SETS IN THE HDF5 FILES
ID="/PartType1/ParticleIDs"
COORDS="/PartType1/Coordinates"
VELS="/PartType1/Velocities"


# ============================================================================
# ============================================================================
# ============================================================================
def traverse_datasets(hdf_file):  
    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            path = '{0}/{1}'.format(prefix, key)
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from h5py_dataset_iterator(item, path)
                
    with h5py.File(hdf_file, 'r') as f:
        for path, _ in h5py_dataset_iterator(f):
            yield path


# ============================================================================
# ============================================================================
           
def index_3d(i,j,k,Nj,Nk):
    return k+Nk*j + Nk*Nj*i

# ============================================================================
# ============================================================================
# given the coordinate x, this function returns the index of its cell
# with two mirrorws ibb and iff. If the cell is sz0 then i->ibb =N, if the cell is N-1, then i->iff =-1
def get_bin (xi, ximin, delta_x, N):

    i = int((xi-ximin)/delta_x)
    if(i==N):
        i=i-1

    ibb=i
    if(i==0):
        ibb=N

    iff=i
    if(i==N-1):
        iff=-1
      
    return (ibb,i,iff)
   
# ============================================================================
# ============================================================================
def mas(mass, x):
    xf=np.abs(x)
    
    if(mass==2):
        w=0.
        if(xf<0.5):
            w=(0.75-x*x)
        elif(xf>=0.5 and xf>1.5):
            w=(0.5*pow(1.5-xf,2));
    
    elif(mass==1):
        w=0.
        if(xf<1.):
            w=1.-xf

    elif(mass==0):
        w=0
        if(xf<0.5):
            w=1.0
            
    return w

# ============================================================================
# ============================================================================
def get_MAS(x,y,z,step):
    i=int(x/delta)
    tx=mas(MAS,x-delta*(i+step))
    j=int(y/delta)
    ty=mas(MAS,y-delta*(j+step))
    k=int(z/delta)
    tz=mas(MAS,z-delta*(k+step))
    return (tx, ty, tz)

# ============================================================================
def get_MAS(mass_scheme,x):
    i=int(x/delta)
    txb=mas(mass_scheme,x-delta*(i-0.5))
    tx =mas(mass_scheme,x-delta*(i+0.5))
    txf=mas(mass_scheme,x-delta*(i+1.5))
    return (txb, tx, txf)
    
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
print("")
print("======================================")
print ("READING HDF5 FILES FROM MINERVA")
print ("Grid assignment")
print("======================================")
print("")


write_bin = "false"
if(write_bin == "true"):
    binfile_dm_x= 'minerva_dm_x.dat'
    fo_bin_x=open(binfile_dm_x,'wb')
    
    binfile_dm_y= 'minerva_dm_y.dat'
    fo_bin_y=open(binfile_dm_y,'wb')
    
    binfile_dm_z= 'minerva_dm_z.dat'
    fo_bin_z=open(binfile_dm_z,'wb')
    
    
    binfile_dm_vx= 'minerva_dm_velx.dat'
    fo_bin_vx=open(binfile_dm_vx,'wb')
    
    binfile_dm_vy= 'minerva_dm_vely.dat'
    fo_bin_vy=open(binfile_dm_vy,'wb')
    
    binfile_dm_vz= 'minerva_dm_velz.dat'
    fo_bin_vz=open(binfile_dm_vz,'wb')

    
dir_path = "/net/deimos/scratch1/balaguera/data/Numerics/Minerva_DM/z1/001/"
#os.path.dirname(os.path.realpath(__file__))
print("")
print ("Reading hdf5 file from directory", dir_path)
print("")


print("")
print("Defining grid-vectos")

vel_field_x=[0 for i in range (0,Ngrid)]
vel_field_y=[0 for i in range (0,Ngrid)]
vel_field_z=[0 for i in range (0,Ngrid)]
ncounts_ngp=[0 for i in range (0,Ngrid)]

if (MAS>0):
    ncounts =[0 for i in range (0,Ngrid)]


print("Done")


for filename in os.listdir(dir_path):
#for ij in range(0,1):

#    filename="snap_001.05.hdf5"
    if filename.endswith(".hdf5"):
    
        with h5py.File(filename, 'r') as f:


            print("")
            print("Reading from file", filename)
           
            print("")
            for dset in traverse_datasets(filename):
                print('Path:', dset)
                print('Shape:', f[dset].shape)
                print('Data type:', f[dset].dtype)
                if(dset==ID):
                    n_parts=len(f[dset])

            print("Done")
            print("")
            print("Number of dm particles = ", n_parts)

            print("")
            print("Extracting coordinates and velocities from hdf5 file:")

                
            x=[]
            x= f[COORDS][:,0]

            y=[]
            y= f[COORDS][:,1]

            z=[]
            z= f[COORDS][:,2]

            vx=[]
            vx= f[VELS][:,0]
            
            vy=[]
            vy= f[VELS][:,1]

            vz=[]
            vz= f[VELS][:,2]

            
            if(write_bin == "true"):           
                print("Writing x-coordinates to binary files")
                pickle.dump(x,file_bin_x)
                print("Writing y-coordinates to binary files")
                pickle.dump(y,file_bin_y)
                print("Writing z-coordinates to binary files")
                pickle.dump(z,file_bin_z)
                print("Writing x-velocities to binary files")
                pickle.dump(x,file_bin_vx)
                print("Writing y-velocities to binary files")
                pickle.dump(y,file_bin_vy)
                print("Writing z-velocities to binary files")
                pickle.dump(z,file_bin_vz)


            print("Done")
            print("")


            if(Ncells_MAX==1):
                print("Interpolating velocity field using MAS = ", MAS)                
                for kp in range (n_parts):
                    i_ind = get_bin(x[kp],xmin,delta, N)
                    j_ind = get_bin(y[kp],ymin,delta, N)
                    k_ind = get_bin(z[kp],zmin,delta, N)
                    index_c=index_3d(i_ind[1],j_ind[1],k_ind[1],N,N)
                    ncounts_ngp[index_c]+=1
                    vel_field_x[index_c]+=vx[kp]
                    vel_field_y[index_c]+=vy[kp]
                    vel_field_z[index_c]+=vz[kp]
                print("Done")                    

            else:
                print("Interpolating velocity field using MAS = ", MAS)                
                for kp in range (n_parts):
                    i_ind = get_bin(x[kp],xmin,delta, N)
                    i=i_ind[1]
                    j_ind = get_bin(y[kp],ymin,delta, N)
                    j=j_ind[1]
                    k_ind = get_bin(z[kp],zmin,delta, N)
                    k=k_ind[1]
                    ncounts_ngp[index_3d(i,j,k,N,N)]+=1
                    x_weight=get_MAS(MAS,x[kp])
                    y_weight=get_MAS(MAS,y[kp])
                    z_weight=get_MAS(MAS,z[kp])
                    
                    for ic in range (0,Ncells_MAX):
                        iaux=i_ind[ic]+ic-1
                        twx=x_weight[ic]
                        for jc in range (0,Ncells_MAX):
                            jaux=j_ind[jc]+jc-1
                            twy=y_weight[jc]
                            for kc in range (0,Ncells_MAX):
                                kaux=k_ind[kc]+kc-1
                                twz=z_weight[kc]
                                index=index_3d(iaux,jaux,kaux,N,N)
                                #                                print (iaux, jaux,kaux,Ngrid, index)
                                weight=twx*twy*twz
                                ncounts[index]+=weight
                                vel_field_x[index]+=vx[kp]*weight
                                vel_field_y[index]+=vy[kp]*weight
                                vel_field_z[index]+=vz[kp]*weight
                print("Done")                                


print("Done")
print("Getting mean velocities in cells")

for i in range(Ngrid):
    if(ncounts_ngp[i]!=0):
        vel_field_x[i]=vel_field_x[i]/ncounts_ngp[i]
        vel_field_y[i]=vel_field_y[i]/ncounts_ngp[i]
        vel_field_z[i]=vel_field_z[i]/ncounts_ngp[i]
    else:
        vel_field_x[i]=0.
        vel_field_y[i]=0.
        vel_field_z[i]=0.
print("Done")

print("Writting to binary files")


binfile_dm_nc_grid= 'minerva_dm_ncounts_'+str(MAS_NAME[MAS])+'.dat'
file_bin_ncounts=open(binfile_dm_nc_grid,'wb')
pickle.dump(ncounts, file_bin_ncounts)

if(MASS>0):
    binfile_dm_nc_grid_ngp= 'minerva_dm_ncounts_ngp.dat'
    file_bin_ncounts_ngp=open(binfile_dm_nc_grid_ngp,'wb')
    pickle.dump(ncoutns_ngp, file_bin_ncounts_ngp)


binfile_dm_vx_grid= 'minerva_dm_velx_grid_'+str(MAS_NAME[MAS])+'.dat'
file_bin_vx=open(binfile_dm_vx_grid,'wb')
pickle.dump(vel_field_x, file_bin_vx)

binfile_dm_vy_grid= 'minerva_dm_vely_grid_'+str(MAS_NAME[MAS])+'.dat'
file_bin_vy=open(binfile_dm_vy_grid,'wb')
pickle.dump(vel_field_y, file_bin_vy)


binfile_dm_vz_grid= 'minerva_dm_velz_grid_MAS_'+str(MAS_NAME[MAS])+'.dat'
file_bin_vz=open(binfile_dm_vz_grid,'wb')
pickle.dump(vel_field_z, file_bin_vz)


file_bin_vx.close()
file_bin_vy.close()
file_bin_vz.close()


if(write_bin == "true"):
    file_bin_x.close()
    file_bin_y.close()
    file_bin_z.close()

        
print("Done")
