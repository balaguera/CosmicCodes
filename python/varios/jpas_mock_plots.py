#!/usr/bin/env python3
import numpy as np
import scipy.interpolate
import scipy.ndimage
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.ticker import AutoMinorLocator
import os.path
import h5py
import pickle
# ============================================================================
#THESE ARE THE DIFFERENT DATA SETS IN THE HDF5 FILES
lwidth=1.0
# ============================================================================
# ============================================================================
grid = plt.GridSpec(4, 4, wspace=0.4, hspace=0.4)
X_min=0
X_max=1
Y_min=0
Y_max=2.5e5


fig=plt.subplot(grid[0, 0])
file_theo='/home/andres/data/Numerics/JPAS/compare_dndz_040_041.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= data[:, 0]
    y= data[:, 1]/500
    plt.plot(x,y,color='grey',linestyle='solid',linewidth=0.8, alpha=0.8,label=r'GalCat40')
    
file_theo='/home/andres/data/Numerics/JPAS/dndz_smooth_bspline.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= data[:, 0]
    y= data[:, 1]/500
    plt.plot(x,y,color='red',linestyle='solid',linewidth=1., alpha=0.3,label=r'dNdz smoothed')
    

for i in range(0,10):
    file_theo='/home/andres/data/Numerics/JPAS/dndz_random.txt_file'+str(i)
    if os.path.isfile(file_theo):
        data=np.loadtxt(file_theo)
        x= data[:, 0]
        y= data[:, 1]/50
        if(i==0):
            plt.plot(x,y,color='green',linestyle='solid',linewidth=0.3, alpha=0.2,label=r'Randoms')
        else:
            plt.plot(x,y,color='green',linestyle='solid',linewidth=0.3, alpha=0.2)
                
plt.legend(fontsize=2)

plt.xlabel(r"$z$", fontsize=4)
plt.ylabel(r"$dN/dz$", fontsize=4)
plt.xscale("linear")
plt.yscale("linear")
fig.yaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)
    
    
# ============================================================================
# ============================================================================
fig=plt.subplot(grid[1, 0])
X_min=-2
X_max=-0.2
Y_min=2.5
Y_max=5.5

file_theo='/home/andres/data/Numerics/JPAS/MOCKS/Pk_fkp_JpasMock_Nft512_CIC_galaxies.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= np.log10(data[:, 0])
    y= np.log10(data[:, 1])
    plt.plot(x,y,color='green',linestyle='solid',linewidth=1.2, alpha=0.8,label=r'CIC')


file_theo='/home/andres/data/Numerics/JPAS/MOCKS/Pk_fkp_JpasMock_Nft512_TSC_galaxies.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= np.log10(data[:, 0])
    y= np.log10(data[:, 1])
    plt.plot(x,y,color='red',linestyle='solid',linewidth=1.2, alpha=0.8,label=r'TSC')



file_theo='/home/andres/data/Numerics/JPAS/MOCKS/Pk_fkp_JpasMock_Nft512_TSC_window'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= np.log10(data[:, 0])
    y= np.log10(data[:, 1])
    plt.plot(x,y,color='black',linestyle='dotted',linewidth=1.2, alpha=0.8,label=r'$|W(k)|^{2}$')

plt.legend(fontsize=2)


    
plt.xlabel(r"$\log_{10}k$", fontsize=4)
plt.ylabel(r"$\log_{10} P_{0}(k)$", fontsize=4)
    
    
plt.xscale("linear")
plt.yscale("linear")
fig.yaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)
       
# ============================================================================
# ============================================================================
fig=plt.subplot(grid[2, 0])
X_min=-2
X_max=-0.2
Y_min=1
Y_max=5.5

file_theo='/home/andres/data/Numerics/JPAS/MOCKS/Pk_fkp_JpasMock_Nft512_TSC_galaxies.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= np.log10(data[:, 0])
    y= np.log10(data[:, 2])
    plt.plot(x,y,color='red',linestyle='solid',linewidth=1.2, alpha=0.8,label=r'$\ell = 2$')


file_theo='/home/andres/data/Numerics/JPAS/MOCKS/Pk_fkp_JpasMock_Nft512_TSC_galaxies.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= np.log10(data[:, 0])
    y= np.log10(data[:, 3])
    plt.plot(x,y,color='blue',linestyle='solid',linewidth=1.2, alpha=0.8,label=r'$\ell = 4$')

    
plt.xlabel(r"$\log_{10}k$", fontsize=4)
plt.ylabel(r"$\log_{10} P_{2}(k)$", fontsize=4)

plt.legend(fontsize=2)

    
plt.xscale("linear")
plt.yscale("linear")
fig.yaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)
    

# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
fig=plt.subplot(grid[0, 1])
X_min=11
X_max=15
Y_min=-700
Y_max=700


file_theo='/home/andres/data/Numerics/JPAS/MOCKS/reduced_gal_cat.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= np.log10(data[:, 1])
    y= (data[:, 4])
    plt.scatter(x,y,c=data[:,0],alpha=0.1, s=0.1, marker='o', linewidth=0.01)

   
plt.xlabel(r"$\log10 M_{h}$", fontsize=4)
plt.ylabel(r"Relative Bias", fontsize=4)
    
plt.xscale("linear")
plt.yscale("linear")
fig.yaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)

# ============================================================================
fig=plt.subplot(grid[1, 1])
X_min=-2
X_max=4
Y_min=-700
Y_max=700

data=np.loadtxt(file_theo)
x= (data[:, 2])
y= (data[:, 4])
plt.scatter(x,y,c=data[:,0],alpha=0.1, s=0.1, marker='o', linewidth=0.01)

   
plt.xlabel(r"Color", fontsize=4)
plt.ylabel(r"Relative Bias", fontsize=4)


    
plt.xscale("linear")
plt.yscale("linear")
fig.yaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)

    
# ============================================================================
fig=plt.subplot(grid[2, 1])
X_min=6
X_max=12
Y_min=-700
Y_max=700

data=np.loadtxt(file_theo)
x= (data[:, 3])
y= (data[:, 4])
plt.scatter(x,y,c=data[:,0],alpha=0.1, s=0.1, marker='o', linewidth=0.01)

    
plt.xlabel(r"$\log M_{\star}$", fontsize=4)
plt.ylabel(r"Relative Bias", fontsize=4)


    
plt.xscale("linear")
plt.yscale("linear")
fig.yaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)



plt.savefig("rand.png",dpi=600, bbox_inches='tight')    