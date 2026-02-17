#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:23:29 2023
@author: Andres Balaguera-Antolínez

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm  
import os.path
import configparser
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
from matplotlib.ticker import AutoMinorLocator                                                                                                                                                                 
from matplotlib.colors import LogNorm                                                                                                                                                                          
import matplotlib                        
import json
import scipy.ndimage
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes                                                                                                                                                   
from mpl_toolkits.axes_grid1 import make_axes_locatable     
from pathlib import Path
import shutil
import subprocess
import matplotlib.colors as colors
from scipy.ndimage.filters import gaussian_filter


def get_title_coords_d(xmin, xmax, ymin, ymax):
    x_title=xmin+0.15*(xmax-xmin)
    y_title=ymin+0.2*(ymax-ymin)
    return x_title, y_title

def get_title_coords(xmin, xmax, ymin, ymax, fracx, fracy):
    x_title=xmin+fracx*(xmax-xmin)
    y_title=ymin+fracy*(ymax-ymin)
    return x_title, y_title

# ==========================================================================
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
    cmap(np.linspace(minval, maxval, n)))
    return new_cmap

# ==========================================================================
# ==========================================================================
# ==========================================================================

class COSMOLIB_PLOTS:

    def __init__(self,file):
        self.parameter_file=file
        print("Reading input file:",file)
        with open(self.parameter_file, 'r',encoding="utf-8") as file:
           self.data = json.load(file)
    
    def plot_mass_dependent_quantities(self, quantity):
        plt.figure()

        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)

        fig = plt.subplot(grid[0, 0]) 
        lwidth=2
        X_min=self.data["mmin"]
        X_max=self.data["mmax"]
        X_max_max=X_max+1
        tick_x=1  # this sets the label every 0.05 in x
        min_locator_x=1
        file_name=self.data["output_file_mass_function"]

        if(quantity=="mass_function"):
            Y_min=1e-29
            Y_max=1e-1
            Y_max_max=2e-1
            index=self.data["index_mass_function"]
            plt.ylabel(r"$n(M)$", fontsize=12)
            plt.title("Halo mass function")
            labels=[r"$n(M)$"]
            plt.yscale('log')
            xtit, ytit=get_title_coords_d(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))
            xtit=pow(10,xtit)
            ytit=pow(10,ytit)
        elif(quantity=="mass_dispersion"):
            Y_min=0
            Y_max=8.0
            Y_max_max=5.0
            index=self.data["index_mass_dispersion"]
            plt.ylabel(r"$\sigma(M)$", fontsize=12)
            plt.title("Mass dispersion")
            labels=[r"$\sigma(M)$"]
            plt.yscale('linear')
            xtit, ytit=get_title_coords_d(np.log10(X_min), np.log10(X_max), Y_min, Y_max)
            xtit=pow(10,xtit)
        elif(quantity=="halo_mass_bias"):
            Y_min=0
            Y_max=4.0
            Y_max_max=5.
            index=self.data["index_halo_mass_bias"]
            index2=self.data["index_effective_halo_mass_bias"]
            plt.ylabel(r"$b(M)$", fontsize=12)
            plt.title("Halo-mass bias")
            labels=[r"$b(M)$", r"$b_{eff}(M)$"]
            plt.yscale('linear')
            xtit, ytit=get_title_coords_d(np.log10(X_min), np.log10(X_max), Y_min, Y_max)
            xtit=pow(10,xtit)

        plt.text(xtit, ytit, "z  = "+str(self.data["redshift"]), fontsize=11)
   
        tick_y=0.5
        min_locator_y=1
        aspect='auto'

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]

        if(quantity=="mass_function" or quantity=="mass_dispersion"):
            data=np.genfromtxt(file_name)     
            x=(data[:,0])
            psc=(data[:,index])
            plt.plot(x,psc,color=colors[0],linestyle=lines_s[0],linewidth=3, alpha=0.6, label=labels[0])
        elif(quantity=="halo_mass_bias"):
            data=np.genfromtxt(file_name)     
            x=(data[:,0])
            psc=(data[:,index])
            plt.plot(x,psc,color=colors[0],linestyle=lines_s[0],linewidth=3, alpha=0.6, label=labels[0])
            psc=(data[:,index2])
            plt.plot(x,psc,color=colors[1],linestyle=lines_s[1],linewidth=3, alpha=0.6, label=labels[1])


        plt.legend(loc="upper right", fontsize=10)
        fig.set_aspect(aspect)

        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)

        plt.xscale('log')

        plt.xlabel(r"$\log M  \,[h  {\rm M_{\odot}}^{-1}]$", fontsize=12)
        plt.grid(alpha=0.2)
        plt.axis([X_min, X_max, Y_min, Y_max])
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        #fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=10, size=8,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=6, size=6,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=10, size=8,labelbottom=1, labeltop=0)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=6, size=6,labelbottom=1, labeltop=0)                                                                                      
   #     fig.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    def plot_density_profile(self, space):
        plt.figure()
        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)

        lwidth=2
        if space == "r" :
            X_min=self.data["rmin_dp"]
            X_max=self.data["rmax_dp"]
            file_name=self.data["output_file_density_profile_k"]

        elif space == "k" :
            X_min=self.data["kmin_dp"]
            X_max=self.data["kmax_dp"]
            file_name=self.data["output_file_density_profile_r"]

        X_max_max=X_max
        tick_x=0.5  # this sets the label every 0.05 in x
        min_locator_x=2
        Y_min=0.01
        Y_max=1e5
        Y_max_max=5
        tick_y=0.5
        min_locator_y=1
        aspect='auto'
        fig = plt.subplot(grid[0, 0]) 

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]
        if space == "r" :
            labels=[r"$\rho(r)$"]
        elif space == "k" :
            labels=[r"$\rho(k)$"]

        plt.title("Dark matter halo density profile")
        xtit, ytit=get_title_coords_d(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))
        plt.text(pow(10,xtit), pow(10,ytit), "z  = "+str(self.data["redshift"])+"    Halo mass " +"{:.3e}".format(self.data["Mass"]), fontsize=11)

        file_name=self.data["output_file_galaxy_power_spectrum"]
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            i=0
            x=(data[:,0])
            psc=(data[:,1])
            plt.plot(x,psc,color=colors[i],linestyle=lines_s[i],linewidth=3, alpha=0.6, label=labels[0])

        plt.legend(loc="upper right", fontsize=11)
       
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
            

        fig.set_xticks([0.01,0.03, 0.1,1])

        plt.xscale('log')
        plt.yscale('log')
        if space == "r" :
            plt.ylabel(r"$\rho(r)$", fontsize=12)
            plt.xlabel(r"$r  \,[{\rm Mpc}\, h^{-1}]$", fontsize=12)
        elif space == "k" :
            plt.ylabel(r"$\rho(k)$", fontsize=12)
            plt.xlabel(r"$k  \,[h  {\rm Mpc}^{-1}]$", fontsize=12)
    
        plt.grid(alpha=0.2)

        plt.axis([X_min, X_max, Y_min, Y_max])
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        #fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=10, size=8,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=6, size=6,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=10, size=8,labelbottom=1, labeltop=0)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=6, size=6,labelbottom=1, labeltop=0)                                                                                      



#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    def plot_power_spectrum_hmodel(self):
        plt.figure()
        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)

        lwidth=2

        X_min=self.data["kmin"]
        X_max=self.data["kmax"]
        X_max_max=X_max
        tick_x=1  # this sets the label every 0.05 in x
        min_locator_x=0.5
        Y_min=0.01
        Y_max=1e5
        Y_max_max=5
        tick_y=0.5
        min_locator_y=1
        aspect='auto'
        fig = plt.subplot(grid[0, 0]) 

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]
        labels=[r"$P_{ss}$", r"$P_{sc}$",r"$P_{2h}$",r"$P$",r"$P_{dm}$",r"$P_{dm,hf}$"]
        plt.title("Galaxy power spectrum")
        xtit, ytit=get_title_coords_d(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))
        plt.text(pow(10,xtit), pow(10,ytit), "z  = "+str(self.data["redshift"]), fontsize=11)

        file_name=self.data["output_file_galaxy_power_spectrum"]
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            for i in 1,2,3,4,5,6:
                x=(data[:,0])
                psc=(data[:,i])
                plt.plot(x,psc,color=colors[i],linestyle=lines_s[i],linewidth=3, alpha=0.6, label=labels[i-1])

        plt.legend(loc="upper right", fontsize=11)
       
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
            

        fig.set_xticks([0.01,0.03, 0.1,1])

        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel(r"$P(k)$", fontsize=12)
        plt.xlabel(r"$k  \,[h  {\rm Mpc}^{-1}]$", fontsize=12)
        plt.grid(alpha=0.2)

        plt.axis([X_min, X_max, Y_min, Y_max])
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        #fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=10, size=8,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=6, size=6,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=10, size=8,labelbottom=1, labeltop=0)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=6, size=6,labelbottom=1, labeltop=0)                                                                                      

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


    def plot_power_spectrum_measurement(self):
        plt.figure()
        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)

        lwidth=2

        X_min=self.data["kmin"]
        X_max=self.data["kmax"]
        X_max_max=X_max
        tick_x=1  # this sets the label every 0.05 in x
        min_locator_x=0.5
        Y_min=1
        Y_max=1e5
        Y_max_max=5
        tick_y=0.5
        min_locator_y=1
        aspect='auto'
        fig = plt.subplot(grid[0, 0]) 

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]
        labels=[r"$P_{0}(k)$", r"$P_{2}(k)$",r"$P_{4}(k)$",r"$P$",r"$P_{dm}$",r"$P_{dm,hf}$"]
        plt.title("Power spectrum in "+self.data["clustering_space"])

        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max), 0.15, 0.15)
        plt.text(pow(10,xtit), pow(10,ytit),str(self.data["sample"]), fontsize=11)
        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max), 0.15, 0.1)
        plt.text(pow(10,xtit), pow(10,ytit),"Cosmological redshift = "+str(self.data["redshift"]), fontsize=11)

        file_name=self.data["output_file"]
        print("Reading file", file_name)
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            for i in 1,2,3,:
                x=(data[:,0])
                psc=(data[:,i])
                plt.scatter(x,psc,color=colors[i],alpha=0.6, label=labels[i-1])

        plt.legend(loc="upper right", fontsize=11)
       
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
            

        fig.set_xticks([0.01,0.03, 0.1,1])

        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel(r"$P(k)$", fontsize=12)
        plt.xlabel(r"$k  \,[h  {\rm Mpc}^{-1}]$", fontsize=12)
        plt.grid(alpha=0.2)

        plt.axis([X_min, X_max, Y_min, Y_max])
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        #fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=10, size=8,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=6, size=6,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=10, size=8,labelbottom=1, labeltop=0)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=6, size=6,labelbottom=1, labeltop=0)                                                                                      


    def plot_power_spectrum_grf(self):
        plt.figure()
        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)

        lwidth=2

        X_min=self.data["kmin"]
        X_max=self.data["kmax"]
        X_max_max=X_max
        tick_x=1  # this sets the label every 0.05 in x
        min_locator_x=0.5
        Y_min=10
        Y_max=2e5
        Y_max_max=5
        tick_y=0.5
        min_locator_y=1
        aspect='auto'
        fig = plt.subplot(grid[0, 0]) 

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]
        labels=[r"$P_{grf}(k)$", r"$\hat{P}_{th}(k)$", r"$\hat{P}_{fa}(k)$", r"$\hat{P}_{mixed}(k)$",r"$P_{th}(k)$",r"$P$",r"$P_{dm}$",r"$P_{dm,hf}$"]
        plt.title(self.data["clustering_space"])

        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max), 0.15, 0.15)
        plt.text(pow(10,xtit), pow(10,ytit),str(self.data["sample"]), fontsize=11)
        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max), 0.15, 0.1)
        plt.text(pow(10,xtit), pow(10,ytit),"Cosmological redshift = "+str(self.data["redshift"]), fontsize=11)

        file_name=self.data["output_file"]
        print("Reading file", file_name)
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            x=(data[:,0])
            psc=data[:,1]
            plt.scatter(x,psc,color=colors[0],alpha=0.5, label=labels[0])
            psc=data[:,2]
            plt.scatter(x,psc,color=colors[1],alpha=0.4, label=labels[1])
            psc=data[:,3]
            plt.scatter(x,psc,color=colors[2],alpha=0.4, label=labels[2])
            psc=data[:,4]
            plt.scatter(x,psc,color=colors[3],alpha=0.4, label=labels[3])

        file_name=self.data["input_file"]
        print("Reading file", file_name)
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            x=(data[:,0])
            psc=(data[:,1])
            plt.plot(x,psc,linewidth=2, linestyle="solid", color=colors[4],alpha=0.7, label=labels[4])


        plt.legend(loc="upper right", fontsize=11)
       
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
            

        fig.set_xticks([0.01,0.03, 0.1,1])

        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel(r"$P(k)$", fontsize=12)
        plt.xlabel(r"$k  \,[h  {\rm Mpc}^{-1}]$", fontsize=12)
        plt.grid(alpha=0.2)

        plt.axis([X_min, X_max, Y_min, Y_max])
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        #fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=10, size=8,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=6, size=6,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=10, size=8,labelbottom=1, labeltop=0)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=6, size=6,labelbottom=1, labeltop=0)                                                                                      

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    def plot_cosmology_distances(self):
        plt.figure()
        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)

        lwidth=2

        X_min=self.data["zmin"]
        X_max=self.data["zmax"]
        X_max_max=X_max
        tick_x=1  # this sets the label every 0.05 in x
        min_locator_x=8
        Y_min=1
        Y_max=1e5
        Y_max_max=5
        tick_y=1.0
        min_locator_y=1
        aspect='auto'
        fig = plt.subplot(grid[0, 0]) 

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]
        labels=[r"Comoving distance", r"Comoving angular diameter distance","Comoving sound horizon","Luminosity distance",r"$P_{dm}$",r"$P_{dm,hf}$"]
        plt.title("Cosmological distances")
        xtit, ytit=get_title_coords_d(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))

        file_name=self.data["output_file"]
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            x= data[:,self.data["column_redshift"]]
            psc= data[:,self.data["column_comoving_distance"]]
            plt.plot(x,psc,color=colors[0],linestyle=lines_s[0],linewidth=3, alpha=0.6, label=labels[0])
            psc= data[:,self.data["column_comoving_angular_diameter_distance"]]
            plt.plot(x,psc,color=colors[1],linestyle=lines_s[1],linewidth=3, alpha=0.6, label=labels[1])
            psc= data[:,self.data["column_comoving_sound_horizon"]]
            plt.plot(x,psc,color=colors[2],linestyle=lines_s[2],linewidth=3, alpha=0.6, label=labels[2])
            psc= data[:,self.data["column_luminosity_distance"]]
            plt.plot(x,psc,color=colors[3],linestyle=lines_s[3],linewidth=3, alpha=0.6, label=labels[3])

        plt.legend(loc="upper left", fontsize=9)
       
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
            

        fig.set_xticks([0.01,0.03, 0.1,1])

        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel(r"Comoving distance Mpc$/h$", fontsize=12)
        plt.xlabel(r"Cosmological redshift", fontsize=12)
        plt.grid(alpha=0.2)

        plt.axis([X_min, X_max, Y_min, Y_max])
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        #fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=10, size=8,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=6, size=6,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=10, size=8,labelbottom=1, labeltop=0)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=6, size=6,labelbottom=1, labeltop=0)                                                                                      

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


    def plot_cosmology_growth(self):
        plt.figure()
        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)

        lwidth=2

        X_min=self.data["zmin"]
        X_max=self.data["zmax"]
        X_max_max=X_max
        tick_x=1  # this sets the label every 0.05 in x
        min_locator_x=8
        Y_min=0
        Y_max=1.
        Y_max_max=5
        tick_y=1.0
        min_locator_y=1
        aspect='auto'
        fig = plt.subplot(grid[0, 0]) 

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]
        labels=[r"Growth factor", r"Comoving angular diameter distance","Comoving sound horizon","Luminosity distance",r"$P_{dm}$",r"$P_{dm,hf}$"]
        plt.title("Cosmological linear growth factor")
        xtit, ytit=get_title_coords_d(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))

        file_name=self.data["output_file"]
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            x= data[:,self.data["column_redshift"]]
            psc= data[:,self.data["column_growth_factor"]]
            plt.plot(x,psc,color=colors[0],linestyle=lines_s[0],linewidth=3, alpha=0.6, label=labels[0])

        plt.legend(loc="upper left", fontsize=9)
       
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
            
        fig.set_xticks([0.01,0.03, 0.1,1])

        plt.xscale('log')
        plt.yscale('linear')
        plt.ylabel(r"Cosmological growth factor", fontsize=12)
        plt.xlabel(r"Cosmological redshift", fontsize=12)
        plt.grid(alpha=0.2)

        plt.axis([X_min, X_max, Y_min, Y_max])
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_ticks_position('both')
        #fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_y))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=10, size=8,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=6, size=6,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=10, size=8,labelbottom=1, labeltop=0)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=6, size=6,labelbottom=1, labeltop=0)                                                                                      


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    def plot_interpolated_field(self):

        plt.rc('text', usetex=True)  


        # set tick width
        mpl.rcParams['xtick.major.size'] = 5
        mpl.rcParams['xtick.major.width'] = 2
        mpl.rcParams['xtick.minor.size'] = 5
        mpl.rcParams['xtick.minor.width'] = 2
        mpl.rcParams['ytick.major.size'] = 5
        mpl.rcParams['ytick.major.width'] = 2
        mpl.rcParams['ytick.minor.size'] = 5
        mpl.rcParams['ytick.minor.width'] = 2
        mpl.rcParams['axes.linewidth'] = 2
        meanfiles = np.fromfile(self.data["output_file"],dtype=np.float32)
        nc = self.data["Nft"]
        L = self.data["Lbox"]
        name = self.data["sample"]
        start = self.data["Initial_slice"]
        num_slices = self.data["Nslices"] 
        end = start+num_slices
        b= meanfiles.reshape((nc,nc,nc))
        bmean=np.mean(b)
        bmin=b.min()
        bmax=b.max()
        type = 1

        if type==1:
            c=b[:,:,start:end].mean(axis=-1) #([:,:,num]= y:z, x:y)  ([:,num,:]= y:z, x:x) ([num,:,:]= y:y, x:x) 
        if type==2:
            c=b[:,start:end,:].mean(axis=1) #([:,:,num]= y:z, x:y)  ([:,num,:]= y:z, x:x) ([num,:,:]= y:y, x:x) 
        if type==3:
            c=b[start:end,:,:].mean(axis=0) #([:,:,num]= y:z, x:y)  ([:,num,:]= y:z, x:x) ([num,:,:]= y:y, x:x)
        c = scipy.ndimage.interpolation.zoom(c ,order=1, zoom=2)
        cmap = plt.get_cmap('nipy_spectral')
        cmap = truncate_colormap(cmap, 0.15, .96)
        cmap.set_under('w')
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        minc=c.min()
        #+0.01
        maxc=c.max() 
        #im = plt.imshow(c,norm=LogNorm(vmin=0.1,vmax=3), interpolation='bilinear',aspect='auto', extent=[0,2500,0,2500],cmap=cmap) #DM
        im = plt.imshow(c,vmin=minc,vmax=maxc, interpolation='bilinear',aspect='auto', extent=[0,L,0,L],cmap=cmap) #DM
        ax.set_aspect('equal')
        plt.xlabel(r'$z[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
        plt.ylabel(r'$x[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
        plt.text(0.25*L, 1.03*L, name, fontsize=20, color='black')

        plt.setp(ax.get_xticklabels(), fontsize=18)
        #plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), fontsize=18)

        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", "5%", pad="2%")
        cb = plt.colorbar(im,orientation='vertical',cax=cax)
        cb.ax.get_yaxis().labelpad = 7.
        cb.set_label(self.data["name"],size=20, rotation=270)
        cb.ax.yaxis.set_tick_params(labelsize=15)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


    def plot_interpolated_field_GRF(self):
        color='gist_earth'
        aspect='auto'
        lwidth = 1.5
        nc = self.data["Nft"]
        L = self.data["Lbox"]
        start = self.data["Initial_slice"]
        num_slices = self.data["Nslices"] 
        end = start+num_slices
        par_x_max=L
        par_x_min=0
        par_x_min_min=0
        par_x_max_max=L
        tick_x = L/5.0
        min_locator_x = 3
        name = self.data["sample"]+r" $k>$"+"{:.1e}".format(self.data["kmax_fa"])
        grid = plt.GridSpec(4,3, wspace=0.02, hspace=-0.1)

        fig=plt.subplot(grid[:2, 2])

        files = self.data["output_file_grf"]
        if os.path.isfile(files):  
            meanfiles = np.fromfile(files,dtype=np.float32)
            b= meanfiles.reshape((nc,nc,nc),order='C')
        sizeb=len(b)
        bmax=np.max(b)
        bmin=np.min(b)
        c=b[:,:,start:end].mean(axis=-1) 
        c = scipy.ndimage.interpolation.zoom(c ,order=2, zoom=1)
        cmap = plt.get_cmap(color)
        cmap = truncate_colormap(cmap, 0.15, .96)
        cmap.set_under('w')
        minc=c.min()
        maxc=c.max()
        im = plt.imshow(c,vmin=minc,vmax=maxc, interpolation='bilinear',aspect='auto', extent=[par_x_min,par_x_max,par_x_min,par_x_max],cmap=cmap) #DM
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
        plt.xscale('linear')
        plt.xticks(np.arange(par_x_min_min,par_x_max_max,tick_x))
        plt.yscale('linear')
        plt.yticks(np.arange(par_x_min_min,par_x_max_max,tick_x))
        plt.axis([par_x_min, par_x_max, par_x_min, par_x_max] )
        fig.xaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
        fig.yaxis.set_ticks_position('both')
        plt.xlabel(r'$x[\mathrm{Mpc}\,h^{-1}]$',fontsize=5)
        plt.Circle( (par_x_min+0.5*(par_x_max-par_x_min), par_x_min+0.5*(par_x_max-par_x_min)), 100, fill=False )
        plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth,labelbottom=1)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth,labelbottom=1)
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("top", "5%", pad="6%")
        cb = plt.colorbar(im,orientation='horizontal',cax=cax, fraction=0.1)
        cb.set_label("GRF",rotation=0, fontsize=7)
        cb.ax.get_yaxis().labelpad = 0
        cb.ax.get_xaxis().labelpad = -14
        cb.ax.yaxis.set_tick_params(labelsize=6)
        cb.ax.get_yaxis().direction = "in"
        cb.ax.get_xaxis().direction = "in"
        cb.ax.tick_params(size=1, labelsize=4, width=0.2, labelrotation=0, pad=-5)
    #------------------------------------

        fig=plt.subplot(grid[:2, 0])

        files = self.data["output_file_fa"]
        if os.path.isfile(files):  
            meanfiles = np.fromfile(files,dtype=np.float32)
            b= meanfiles.reshape((nc,nc,nc),order='C')
        sizeb=len(b)
        c=b[:,:,start:end].mean(axis=-1)
        c = scipy.ndimage.interpolation.zoom(c ,order=2, zoom=1)
        cmap = plt.get_cmap(color)
        cmap = truncate_colormap(cmap, 0.15, .96)
        cmap.set_under('w')
        minc=c.min()
        maxc=c.max()
        im = plt.imshow(c,vmin=minc,vmax=maxc, interpolation='bilinear',aspect='auto', extent=[par_x_min,par_x_max,par_x_min,par_x_max],cmap=cmap) #DM
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
        plt.xscale('linear')
        plt.xticks(np.arange(par_x_min_min,par_x_max_max,tick_x))
        plt.yscale('linear')
        plt.yticks(np.arange(par_x_min_min,par_x_max_max,tick_x))
        plt.axis([par_x_min, par_x_max, par_x_min, par_x_max] )
        fig.xaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
        fig.yaxis.set_ticks_position('both')
        plt.ylabel(r'$y[\mathrm{Mpc}\,h^{-1}]$',fontsize=5)
        plt.xlabel(r'$x[\mathrm{Mpc}\,h^{-1}]$',fontsize=5)
        plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth,labelleft=1)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth,labelleft=1)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth,labelbottom=1)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth,labelbottom=1)
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("top", "5%", pad="6%")
        cb = plt.colorbar(im,orientation='horizontal',cax=cax, fraction=0.1)
        cb.ax.get_xaxis().labelpad = -14
        cb.set_label("FA",rotation=0, fontsize=7)
        cb.ax.yaxis.set_tick_params(labelsize=6)
        cb.ax.get_yaxis().direction = "in"
        cb.ax.get_xaxis().direction = "in"
        cb.ax.tick_params(size=1, labelsize=4, width=0.2, labelrotation=0, pad=-5)
    #------------------------------------

        fig=plt.subplot(grid[:2, 1])
    
        files = self.data["output_file_mixed"]
        if os.path.isfile(files):  
            meanfiles = np.fromfile(files,dtype=np.float32)
            b= meanfiles.reshape((nc,nc,nc),order='C')
        sizeb=len(b)
        c=b[:,:,start:end].mean(axis=-1)
        c = scipy.ndimage.interpolation.zoom(c ,order=2, zoom=1)
        cmap = plt.get_cmap(color)
        cmap = truncate_colormap(cmap, 0.15, .96)
        cmap.set_under('w')
        minc=c.min()
        maxc=c.max()
        im = plt.imshow(c,vmin=minc,vmax=maxc, interpolation='bilinear',aspect='auto', extent=[par_x_min,par_x_max,par_x_min,par_x_max],cmap=cmap) #DM
        fig.set_aspect(aspect)
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
        plt.xscale('linear')
        plt.xticks(np.arange(par_x_min_min,par_x_max_max,tick_x))
        plt.yscale('linear')
        plt.yticks(np.arange(par_x_min_min,par_x_max_max,tick_x))
        plt.axis([par_x_min, par_x_max, par_x_min, par_x_max] )
        fig.xaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
        fig.xaxis.set_ticks_position('both')
        fig.yaxis.set_minor_locator(AutoMinorLocator(min_locator_x))
        fig.yaxis.set_ticks_position('both')
        plt.xlabel(r'$x[\mathrm{Mpc}\,h^{-1}]$',fontsize=5)
        plt.tick_params(direction='in',axis='y', which='major', labelsize=3, size=3, width=lwidth,labelleft=0)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=2, size=2, width=lwidth,labelleft=0)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=3, size=3, width=lwidth,labelbottom=1)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=2, size=2, width=lwidth,labelbottom=1)
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("top", "5%", pad="6%")
        cb = plt.colorbar(im,orientation='horizontal',cax=cax, fraction=0.1)
        cb.ax.get_xaxis().labelpad = -14
        cb.set_label(r"Mixed",rotation=0, fontsize=7)
        cb.ax.yaxis.set_tick_params(labelsize=6)
        cb.ax.get_yaxis().direction = "in"
        cb.ax.get_xaxis().direction = "in"
        cb.ax.tick_params(size=1, labelsize=4, width=0.2, labelrotation=0, pad=-5)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    def plot_halo_mass_plane(self):
        cmap_two="RdBu"
        grid = plt.GridSpec(1,1, wspace=0.05, hspace=0.05)
        fraction=1.0
        Nbins=120
        lwidth=2
        sigma_gaus=1.0
        interpol='spline16'
        fig = plt.subplot(grid[0,0]) 
        pathtofile=self.data["output_new_catalog"]
        print("Reading file ", pathtofile)
        file_mean=pathtofile
        data_z1=np.loadtxt(file_mean)
        c_levels=[0.5, 1.5, 2.5, 3.0, 3.5]

        Ntracers=len(data_z1)
        print("Number of tracers read", Ntracers)
        i_propy=self.data["column_bias"]  #Bias
        i_propx=self.data["column_mass"]  #Mass
        mass=data_z1[:,i_propx]
        bias=data_z1[:,i_propy]
        show_contours=1
        show_mean=1
        labeltop=0
        labelright=0
        labelleft=0
        labelbottom=1  
        labelleft=1
        counter=0

        Nbins_low=10
        X_min=min(mass)   
        X_max=max(mass)   
        Y_min=-15   
        Y_max=15   
  
        H, xedges,yedges=np.histogram2d(data_z1[:, i_propx], data_z1[:, i_propy],  range=[[X_min,X_max],[Y_min,Y_max]], bins=Nbins, density=False)
        H=H.T
        H=np.log10(H+1)
        minH=np.min(H)
        maxH=np.max(H)
        delta_x=(X_max-X_min)/float(Nbins_low)
        delta_y=(Y_max-Y_min)/float(Nbins_low)
        delta_x_a=(X_max-X_min)/float(Nbins)
        delta_y_a=(Y_max-Y_min)/float(Nbins)
   
        xbins=[X_min+(i+0.5)*delta_x for i in range(0,Nbins_low)]
        ybins=[Y_min+(i+0.5)*delta_y for i in range(0,Nbins_low)]
        
        
        xbins_a=[X_min+(i+0.5)*delta_x_a for i in range(0,Nbins)]
        ybins_a=[Y_min+(i+0.5)*delta_y_a for i in range(0,Nbins)]
        

        mean_prop=[0 for i in range(0,Nbins_low)]
        Ntracers_bin=[0 for i in range(0,Nbins_low)]
        for i in range(0,Ntracers):
            if(data_z1[i, i_propx]>=X_min and data_z1[i, i_propx]<X_max):
                if(data_z1[i, i_propx]>=X_min and data_z1[i, i_propx]<X_max):
                    inx=int(np.floor((data_z1[i, i_propx]-X_min)/delta_x))
                    mean_prop[inx]+=data_z1[i, i_propy]
                    Ntracers_bin[inx]+=1
        f_x=[]
        f_y=[]
        for i in range(0,len(mean_prop)):
            if(Ntracers_bin[i]>0):
                mean_prop[i]=mean_prop[i]/Ntracers_bin[i]
                f_x.append(xbins[i])
                f_y.append(mean_prop[i])

        var_prop=[0 for i in range(0,Nbins_low)]
        for i in range(0,Ntracers):
            if(data_z1[i, i_propx]>=X_min and data_z1[i, i_propx]<X_max):
                inx=int(np.floor((data_z1[i, i_propx]-X_min)/delta_x))
                var_prop[inx]+=pow(data_z1[i, i_propy]-mean_prop[inx],2)

        err_prop=[] 
        for i in range(0,len(mean_prop)):
            if(Ntracers_bin[i]>0):
                err_prop.append(fraction*np.sqrt(var_prop[i])/Ntracers_bin[i]) 

        label=r"UNITsim $\langle b_{hm}|m_{\rm vir}\rangle $"
        cmap_two="BuGn"

        plt.errorbar(f_x,f_y, yerr=err_prop, xerr=0.47*delta_x, fmt="none", ecolor="darkred", color="black", alpha=0.7, linewidth=0.9,capsize=2)    
        plt.scatter(f_x,f_y, color="black", alpha=1, linestyle="solid", linewidth=0.5, s=20, label="Mean halo-mass bias")    
        plt.errorbar(f_x,f_y, yerr=err_prop, xerr=0.47*delta_x,  fmt="bo", ecolor="black", linestyle="solid", alpha=0.7, linewidth=1.2,capsize=3,ms=5, marker="o")               


        #PLOT THE MEASUREMENTS FROM POWER SPECTRRU,
        plt.legend(loc="upper left", fontsize=7)

        
        if(show_contours):
            Hg=data = gaussian_filter(H, sigma_gaus)
            plt.contour(xbins_a, ybins_a, Hg,levels=c_levels,alpha=0.5, colors="grey",origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, linewidths=.6, inline=True)
        
        im=plt.imshow(H, interpolation=interpol, cmap=cmap_two, origin='lower', extent=[X_min,X_max, Y_min, Y_max],vmin=minH,vmax=maxH, alpha=0.7)
        cbar=plt.colorbar(im, fraction=0.3)
        cbar.ax.tick_params(labelsize=4)
        cbar.set_label(r"$log (N_{H})$", rotation=270, fontsize=7)
        cbar.ax.get_yaxis().labelpad = -27
        cbar.ax.get_yaxis().direction = "in"
        cbar.ax.tick_params(size=1, labelsize=5, width=lwidth, labelrotation=0)

        
        fig.set_aspect('auto')
        for axis in ['top','bottom','left','right']:
            fig.spines[axis].set_linewidth(lwidth)
        plt.xscale("linear")
        fig.xaxis.set_ticks_position('both')
        plt.yscale("linear")
        fig.yaxis.set_ticks_position('both')
        plt.axis([ 1.02*X_min, X_max , Y_min, Y_max])
        plt.xticks(np.arange(X_min ,X_max,1))
        plt.yticks(np.arange(Y_min ,Y_max,5))
        
        plt.xlabel("Mass $\log M$", fontsize=10)
        plt.ylabel("Halo bias", fontsize=10)
        
        fig.xaxis.set_minor_locator(AutoMinorLocator(4))
        fig.yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.tick_params(direction='in',axis='y', which='major', labelsize=7, size=6, width=lwidth,labelleft=labelleft, labelright=labelright)
        plt.tick_params(direction='in',axis='y', which='minor', labelsize=3, size=3, width=lwidth,labelleft=labelleft, labelright=labelright)
        plt.tick_params(direction='in',axis='x', which='major', labelsize=7, size=6, width=lwidth,labelbottom=1, labeltop=labeltop)
        plt.tick_params(direction='in',axis='x', which='minor', labelsize=3, size=3, width=lwidth,labelbottom=1, labeltop=labeltop)
               
        plt.tight_layout()


#========================================================================================
#========================================================================================
#========================================================================================
#========================================================================================
#========================================================================================
print("=========================================================")
print("Python script to plot functions computed in CosmicCodes")
print("=========================================================")



parser = argparse.ArgumentParser(description="Plots from CosmicCodes")
parser.add_argument(
    "json_input_file",
    help="Path to the json input file"
)
args = parser.parse_args()
plotCosmo = COSMOLIB_PLOTS(args.json_input_file)

if plotCosmo.data.get("show_mass_function_plot", False):
    print("Showing mass function")
    plotCosmo.plot_mass_dependent_quantities("mass_function")

if plotCosmo.data.get("show_mass_dispersion_plot", False):
    print("Showing mass dispersion")
    plotCosmo.plot_mass_dependent_quantities("mass_dispersion")

if plotCosmo.data.get("show_halo_mass_bias_plot", False):
    print("Showing halo-bias fit")
    plotCosmo.plot_mass_dependent_quantities("halo_mass_bias")

if plotCosmo.data.get("show_density_profile_r_plot", False):
    print("Showing density profile (r)")
    plotCosmo.plot_density_profile("r")

if plotCosmo.data.get("show_density_profile_k_plot", False):
    print("Showing density profile (k)")
    plotCosmo.plot_density_profile("k")

if plotCosmo.data.get("show_power_spectrum_plot", False):
    print("Showing halo model galaxy power spectrum")
    plotCosmo.plot_power_spectrum_hmodel()

if plotCosmo.data.get("cosmology", False):
    print("Showing cosmology")
    plotCosmo.plot_cosmology_distances()

if plotCosmo.data.get("growth", False):
    print("Showing growth factor")
    plotCosmo.plot_cosmology_growth()

if plotCosmo.data.get("power_spectrum_measurement", False):
    print("Showing measurement of power spectrum")
    plotCosmo.plot_power_spectrum_measurement()

if plotCosmo.data.get("show_interpolated_field", False):
    print("Showing intepolated field")
    plotCosmo.plot_interpolated_field()

if plotCosmo.data.get("power_spectrum_grf", False):
    print("Showing measurement of power spectrum")
    plotCosmo.plot_power_spectrum_grf()

if plotCosmo.data.get("grf", False):
    print("Showing intepolated field")
    plotCosmo.plot_interpolated_field()

if plotCosmo.data.get("fa", False):
    print("Showing intepolated field")
    plotCosmo.plot_interpolated_field_GRF()

if plotCosmo.data.get("show_individual_halo_mass_bias", False):
    print("Showing indivdual halo mass bias relation")
    plotCosmo.plot_halo_mass_plane()

if plotCosmo.data.get("show_bias_field", False):
    print("Showing interpolated bias")
    plotCosmo.plot_interpolated_field()


plt.show()