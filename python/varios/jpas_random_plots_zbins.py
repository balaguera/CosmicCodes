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
interpol='spline16'
cmap_two="RdBu_r"
sigma_gaus=2
c_levels=[1,10, 20]
Nbins = 200


z_bins_min=[0,0.0,0.3,0.6]
z_bins_max=[1,0.3,0.6,1.0]

grid = plt.GridSpec(len(z_bins_min), 5, wspace=0.1, hspace=0.1)

# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
Y_min=0
Y_max=3
X_min=8
X_max=12



file_theo = "/home/andres/data/Numerics/JPAS/MOCKS/reduced_gal_cat.txt"
data_gal=np.loadtxt(file_theo)
file_ran = "/home/andres/data/Numerics/JPAS/RANDOM/xaa"
data_ran=np.loadtxt(file_ran)


delta_x_a=(X_max-X_min)/float(Nbins)
delta_y_a=(Y_max-Y_min)/float(Nbins)    
xbins_a=[X_min+(i+0.5)*delta_x_a for i in range(0,Nbins)]
ybins_a=[Y_min+(i+0.5)*delta_y_a for i in range(0,Nbins)]


for iz in range(len(z_bins_max)):

    fig=plt.subplot(grid[iz, 0])

    col_gal=[]
    ms_gal=[]
    for i in range(len(data_gal)):
        if(data_gal[i,0]<z_bins_max[iz] and data_gal[i,0]>=z_bins_min[iz]):
            col_gal.append(data_gal[i,2])
            ms_gal.append(data_gal[i,3])


    plt.text(12.2, .5,  str(z_bins_min[iz])+"$\leq z <$"+str(z_bins_max[iz]), rotation=270, fontsize=4)


    #For the galcat we show contours and ima
    H, xedges,yedges=np.histogram2d(ms_gal, col_gal,  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
    H=H.T
    H=gaussian_filter(H, sigma_gaus)
    minH=np.min(H)
    maxH=np.max(H)
    im=plt.imshow(H, interpolation=interpol, cmap=cmap_two, origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, alpha=0.7)
    
    H, xedges,yedges=np.histogram2d(ms_gal, col_gal,  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
    H=H.T
    H=gaussian_filter(H, sigma_gaus)
    minH=np.min(H)
    maxH=np.max(H)
    plt.contour(xbins_a, ybins_a, H,levels=c_levels,alpha=0.3, colors="black",origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, linewidths=.3, inline=True,label="GalCat")

    
    #For the randoms we show contours:
    col_ran=[]
    ms_ran=[]
    for i in range(len(data_ran)):
        if(data_ran[i,0]<z_bins_max[iz] and data_ran[i,0]>=z_bins_min[iz]):
            col_ran.append(data_ran[i,3])
            ms_ran.append(data_ran[i,4])
        
    
    
    
    H, xedges,yedges=np.histogram2d(ms_ran, col_ran,  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
    H=H.T
    H=gaussian_filter(H, sigma_gaus)
    minH=np.min(H)
    maxH=np.max(H)
    plt.contour(xbins_a, ybins_a, H,levels=c_levels,alpha=0.5, colors="black",origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, linewidths=.6, inline=True,label="Random")
        
    if(iz==0):
        plt.legend(fontsize=4)    
    
    labelbottom=0
    if (iz==len(z_bins_min)-1):
        labelbottom=1

    if (labelbottom>0):
        plt.xlabel(r"$\log M_{\star}$", fontsize=5)
    plt.ylabel(r"Color", fontsize=5)

        
    fig.set_aspect('auto')
    for axis in ['top','bottom','left','right']:
        fig.spines[axis].set_linewidth(lwidth)
    fig.yaxis.set_ticks_position('both')
    fig.xaxis.set_ticks_position('both')
    plt.axis([X_min, X_max, Y_min, Y_max])
    plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth)
    plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
    plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth, labelbottom=labelbottom)
    plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=labelbottom)
    





plt.savefig("/home/andres/data/Numerics/JPAS/random_jpas_zbins.png",dpi=600, bbox_inches='tight')    