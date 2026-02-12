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
from scipy.ndimage.filters import gaussian_filter

import os.path
import h5py
import pickle
# ============================================================================
#THESE ARE THE DIFFERENT DATA SETS IN THE HDF5 FILES
lwidth=1.0
# ============================================================================
# ============================================================================
grid = plt.GridSpec(2, 3, wspace=0.25, hspace=0.3)
interpol='spline16'
cmap_two="RdBu_r"
sigma_gaus=2.0
c_levels=[1, 10, 50, 80]
Nbins = 200




# ============================================================================
X_min=0
X_max=1
Y_min=0
Y_max=2.5e5
fig=plt.subplot(grid[1, 0])
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
                
plt.legend(fontsize=4)

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
X_min=-1
X_max=4
Y_min=0
Y_max=5e5
fig=plt.subplot(grid[1, 1])
file_theo='/home/andres/data/Numerics/JPAS/MOCKS/dNdc.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= data[:, 0]
    y= data[:, 1]
    plt.plot(x,y,color='grey',linestyle='dotted',linewidth=1.2, alpha=0.8,label=r'GalCat40')

for i in range(0,10):
    file_theo='/home/andres/data/Numerics/JPAS/MOCKS/dNdc.txt_ranfile'+str(i)
    if os.path.isfile(file_theo):
        data=np.loadtxt(file_theo)
        x= data[:, 0]
        y= data[:, 1]
        if(i==0):
            plt.plot(x,y,color='green',linestyle='solid',linewidth=0.8, alpha=0.2,label=r'Randoms')
        else:
            plt.plot(x,y,color='green',linestyle='solid',linewidth=0.8, alpha=0.2)
                
plt.legend(fontsize=4)




plt.xlabel(r"$Color$", fontsize=4)
plt.ylabel(r"$dN/dc$", fontsize=4)
fig.set_aspect('auto')
for axis in ['top','bottom','left','right']:
    fig.spines[axis].set_linewidth(lwidth)
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)
    
# ============================================================================
X_min=6
X_max=13
Y_min=0
Y_max=7e5
fig=plt.subplot(grid[1, 2])
file_theo='/home/andres/data/Numerics/JPAS/MOCKS/dNdMs.txt'
if os.path.isfile(file_theo):
    data=np.loadtxt(file_theo)
    x= data[:, 0]
    y= data[:, 1]
    plt.plot(x,y,color='grey',linestyle='dotted',linewidth=1.2, alpha=0.8,label=r'GalCat40')
    
    

for i in range(0,10):
    file_theo='/home/andres/data/Numerics/JPAS/MOCKS/dNdMs.txt_ranfile'+str(i)
    if os.path.isfile(file_theo):
        data=np.loadtxt(file_theo)
        x= data[:, 0]
        y= data[:, 1]
        if(i==0):
            plt.plot(x,y,color='green',linestyle='solid',linewidth=0.3, alpha=0.2,label=r'Randoms')
        else:
            plt.plot(x,y,color='green',linestyle='solid',linewidth=0.8, alpha=0.2)
                

plt.xlabel(r"$\log M_{\star}$", fontsize=4)
plt.ylabel(r"$dN/dM_{s}$", fontsize=4)
plt.xscale("linear")
plt.yscale("linear")
fig.set_aspect('auto')
for axis in ['top','bottom','left','right']:
    fig.spines[axis].set_linewidth(lwidth)
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_ticks_position('both')

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
fig=plt.subplot(grid[0, 0])
X_min=-1
X_max=4
Y_min=6
Y_max=13

z_bins_min=[0,0.0,0.2,0.4,0.6,0.8]
z_bins_max=[1,0.2,0.4,0.6,0.8,1.0]
col_gal=[]
ms_gal=[]


file_theo = "/home/andres/data/Numerics/JPAS/MOCKS/reduced_gal_cat.txt"
data=np.loadtxt(file_theo)
H, xedges,yedges=np.histogram2d(data[:, 2], data[:, 3],  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
H=H.T
H=data = gaussian_filter(H, sigma_gaus)
minH=np.min(H)
maxH=np.max(H)
im=plt.imshow(H, interpolation=interpol, cmap=cmap_two, origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, alpha=0.7, label="GalCat40")


delta_x_a=(X_max-X_min)/float(Nbins)
delta_y_a=(Y_max-Y_min)/float(Nbins)    
xbins_a=[X_min+(i+0.5)*delta_x_a for i in range(0,Nbins)]
ybins_a=[Y_min+(i+0.5)*delta_y_a for i in range(0,Nbins)]
file_ran = "/home/andres/data/Numerics/JPAS/RANDOM/xaa"
col_ran=[]
ms_ran=[]
data_ran=np.loadtxt(file_ran)



H, xedges,yedges=np.histogram2d(data_ran[:, 3], data_ran[:, 4],  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
H=H.T
H=gaussian_filter(H, sigma_gaus)
minH=np.min(H)
maxH=np.max(H)
plt.contour(xbins_a, ybins_a, H,levels=c_levels,alpha=0.5, colors="black",origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, linewidths=.6, inline=True,label="Random")
    

plt.ylabel(r"$\log M_{\star}$", fontsize=4)
plt.xlabel(r"Color", fontsize=4)


    
fig.set_aspect('auto')
for axis in ['top','bottom','left','right']:
    fig.spines[axis].set_linewidth(lwidth)
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)

# ============================================================================
fig=plt.subplot(grid[0, 1])
X_min=0
X_max=1
Y_min=6
Y_max=13
Nbins = 100

file_theo = "/home/andres/data/Numerics/JPAS/MOCKS/reduced_gal_cat.txt"
data=np.loadtxt(file_theo)
H, xedges,yedges=np.histogram2d(data[:, 0], data[:, 3],  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
H=H.T
H=gaussian_filter(H, sigma_gaus)
minH=np.min(H)
maxH=np.max(H)
im=plt.imshow(H, interpolation=interpol, cmap=cmap_two, origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, alpha=0.7)


delta_x_a=(X_max-X_min)/float(Nbins)
delta_y_a=(Y_max-Y_min)/float(Nbins)    
xbins_a=[X_min+(i+0.5)*delta_x_a for i in range(0,Nbins)]
ybins_a=[Y_min+(i+0.5)*delta_y_a for i in range(0,Nbins)]
#file_theo = "/home/andres/data/Numerics/JPAS/RANDOM/Random_cat_ranfile4"
#data=np.loadtxt(file_theo)
Ha, xedges,yedges=np.histogram2d(data_ran[:, 0], data_ran[:, 4],  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
Ha=Ha.T
Ha=gaussian_filter(Ha, sigma_gaus)
minH=np.min(Ha)
maxH=np.max(Ha)
plt.contour(xbins_a, ybins_a, Ha,levels=c_levels,alpha=0.5, colors="black",origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, linewidths=.6, inline=True)
   

plt.ylabel(r"$\log M_{\star}$", fontsize=4)
plt.xlabel(r"Redshift", fontsize=4)
fig.set_aspect('auto')
for axis in ['top','bottom','left','right']:
    fig.spines[axis].set_linewidth(lwidth)
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_ticks_position('both')
    
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)
        

# ============================================================================
fig=plt.subplot(grid[0, 2])
X_min=0
X_max=1
Y_min=-1
Y_max=4
Nbins = 100

file_theo = "/home/andres/data/Numerics/JPAS/MOCKS/reduced_gal_cat.txt"
data=np.loadtxt(file_theo)
H, xedges,yedges=np.histogram2d(data[:, 0], data[:, 2],  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
H=H.T
H=gaussian_filter(H, sigma_gaus)
minH=np.min(H)
maxH=np.max(H)
im=plt.imshow(H, interpolation=interpol, cmap=cmap_two, origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, alpha=0.7)


delta_x_a=(X_max-X_min)/float(Nbins)
delta_y_a=(Y_max-Y_min)/float(Nbins)    
xbins_a=[X_min+(i+0.5)*delta_x_a for i in range(0,Nbins)]
ybins_a=[Y_min+(i+0.5)*delta_y_a for i in range(0,Nbins)]
#file_theo = "/home/andres/data/Numerics/JPAS/RANDOM/xvi"
#data=np.loadtxt(file_theo)
H, xedges,yedges=np.histogram2d(data_ran[:, 0], data_ran[:, 3],  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
H=H.T
H=gaussian_filter(H, sigma_gaus)
minH=np.min(H)
maxH=np.max(H)
plt.contour(xbins_a, ybins_a, H,levels=c_levels,alpha=0.5, colors="black",origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, linewidths=.6, inline=True)
    

plt.ylabel(r"Color", fontsize=4)
plt.xlabel(r"Redshift", fontsize=4)


    
fig.set_aspect('auto')
for axis in ['top','bottom','left','right']:
    fig.spines[axis].set_linewidth(lwidth)
fig.yaxis.set_ticks_position('both')
fig.xaxis.set_ticks_position('both')
plt.axis([X_min, X_max, Y_min, Y_max])
plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=1)
plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)




plt.savefig("random_jpas.png",dpi=600, bbox_inches='tight')    