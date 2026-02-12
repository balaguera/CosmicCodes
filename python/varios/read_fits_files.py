#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 08:25:57 2024
Reading fits catalogs and converting to asciii to be read by Cosmicodes.
@author: andres balaguera
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

lwidth= 1.2

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def fits_to_ascii(path, infile):
    gal_file=path+infile
    uri=gal_file
    hdul = fits.open(uri+".fits")
    hdul.info()
    ngals=len(hdul[1].data)
    print ("Found ",ngals," lines")
    outfile = open(gal_file+".txt", 'w') 
    headers=hdul[0].header
    print (headers)


    print ("Writting to ascii")
    for i in range(ngals):
#        outfile.write(str(hdul[1].data[i][0])+" "+str(hdul[1].data[i][1])+" "+str(hdul[1].data[i][2])+" "+str(hdul[1].data[i][3])+" "+str(hdul[1].data[i][4])+ "\n") #radndom
#        outfile.write(str(hdul[1].data[i][5])+" "+str(hdul[1].data[i][1]/hdul[1].data[i][11])+ "\n")    #radndom dndz
        outfile.write(str(hdul[1].data[i][0])+" "+str(hdul[1].data[i][1])+" "+str(hdul[1].data[i][2])+" "+str(hdul[1].data[i][3])+" "+str(hdul[1].data[i][4])+" "+str(hdul[1].data[i][5])+" "+str(hdul[1].data[i][6])+" "+str(hdul[1].data[i][7])+ "\n") #power

    outfile.close()

    print ("Done")
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_power():
    gal_file="/home/andres/data/Numerics/Euclid/FlagShip/Power/PK_flux2_FS230degCircle_m3_smL5_obsz_z0.9-1.1_rebinned"
    uri=gal_file
    hdul = fits.open(uri+".fits")
    hdul.info()
    nlines=len(hdul[1].data)

    grid = plt.GridSpec(3,2, wspace=0.01, hspace=0.1)


    fig=plt.subplot(grid[:2, 0])

    X_min=0
    X_max=0.2
    
    Y_min=3
    Y_max=5
    Y_max_max=5
    tick_y=0.5
    min_locator_y=4
    
    kv=[]
    pk=[]
    pq=[]
    ph=[]
    for i in range(0,nlines):
        kv.append(hdul[1].data[i][0])
        pk.append(np.log10(hdul[1].data[i][2]))
        pq.append(np.log10(hdul[1].data[i][4]))
        ph.append(np.log10(hdul[1].data[i][6]))
        
    plt.plot(kv,pk,color='blue',linestyle='solid',linewidth=1.5, alpha=0.8,label=r'PK_flux2_FS230degCircle_m3_smL5_obsz_z0.9-1.1_rebinned, $\ell=0$')
    plt.plot(kv,pq,color='blue',linestyle='solid',linewidth=1.5, alpha=0.8,label=r'PK_flux2_FS230degCircle_m3_smL5_obsz_z0.9-1.1_rebinned, $\ell=2$')
    
    file_theo='/home/andres/data/Numerics/Euclid/FlagShip/Power/Pk_ys_FlagShip_Nft512_TSC_galaxies_z0.9_z1.1.txt'
    data=np.loadtxt(file_theo)
    x= data[:, 0]
    y= np.log10(data[:, 1])
    z=  np.log10(data[:, 2])
    plt.plot(x,y,color='red',linestyle='dashed',linewidth=1.5, alpha=0.7,label=r'ABA $\ell=0$')
    plt.plot(x,z,color='red',linestyle='dashed',linewidth=1.5, alpha=0.7,label=r'ABA $\ell=2$')
    
    fig.set_aspect('auto')
    for axis in ['top','bottom','left','right']:
        fig.spines[axis].set_linewidth(1)
    
    plt.legend(fontsize=3)
    
    plt.xscale("log")
    fig.xaxis.set_ticks_position('both')
    
    plt.yscale("linear")
    fig.yaxis.set_ticks_position('both')
    plt.axis([X_min, X_max, Y_min, Y_max])

    plt.grid(alpha=0.5)

    
#    plt.xlabel(r"k ", fontsize=8)
    plt.ylabel(r"$\log_{10} P(k)$", fontsize=8)
    
#    plt.xticks(np.arange(X_min ,X_max,tick_x))
    plt.yticks(np.arange(Y_min ,Y_max,tick_y))
    
#    fig.xaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
    fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
    
    plt.tick_params(direction='in',axis='y', which='major', labelsize=4, size=6, width=lwidth, labelbottom=0,labelleft=1)
    plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
    plt.tick_params(direction='in',axis='x', which='major', labelsize=4, size=6, width=lwidth, labelbottom=0)
    plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=0)



    #RATIOS
    fig=plt.subplot(grid[2, 0])

    
    Y_min=0.9
    Y_max=1.1
    Y_max_max=1
    tick_y=0.05
    min_locator_y=4


    
    kv=[]
    pk=[]
    pq=[]
    ph=[]
    for i in range(0,len(x)):
        kv.append(hdul[1].data[i][0])
        pk.append(hdul[1].data[i][2]/pow(10,y[i]))
        pq.append(hdul[1].data[i][4]/pow(10,z[i]))
        
    plt.plot(kv,pk,color='blue',linestyle='solid',linewidth=1.5, alpha=0.8,label=r'$\ell=0$')
    plt.plot(kv,pq,color='red',linestyle='solid',linewidth=1.5, alpha=0.8,label=r'$\ell=2$')
    
    
    fig.set_aspect('auto')
    for axis in ['top','bottom','left','right']:
        fig.spines[axis].set_linewidth(1)
    
    plt.legend(fontsize=4)
    
    plt.xscale("log")
    fig.xaxis.set_ticks_position('both')
    
    plt.yscale("linear")
    fig.yaxis.set_ticks_position('both')
    plt.axis([X_min, X_max, Y_min, Y_max])
    

    plt.grid(alpha=0.5)

    plt.xlabel(r"k ", fontsize=8)
    plt.ylabel(r"Ratio Official to ABA ", fontsize=6)
    
#    plt.xticks(np.arange(X_min ,X_max,tick_x))
    plt.yticks(np.arange(Y_min ,Y_max,tick_y))
    
#    fig.xaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
    fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
    
    plt.tick_params(direction='in',axis='y', which='major', labelsize=4, size=6, width=lwidth)
    plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth)
    plt.tick_params(direction='in',axis='x', which='major', labelsize=4, size=6, width=lwidth, labelbottom=1)
    plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth, labelbottom=1)

    
    plt.savefig("/home/andres/data/Numerics/Euclid/FlagShip/Power/comp.png",dpi=600, bbox_inches='tight')


#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


path="/home/andres/data/Numerics/Euclid/FlagShip/"
#fits_to_ascii(path,"dndz_flux2_flagship_FS2.1.10b_FS230degCircle_m3_obsz")
fits_to_ascii(path,"Power/PK_flux2_FS230degCircle_m3_smL5_obsz_z0.9-1.1_rebinned")
#fits_to_ascii(path,"data_flux2_flagship_FS2.1.10b_FS230degCircle_m3_obsz_z1.1-1.3")
#fits_to_ascii(path,"random_flux2_flagship_FS2.1.10b_FS230degCircle_m3_smL5_obsz_50x_z1.1-1.3")
plot_power()
