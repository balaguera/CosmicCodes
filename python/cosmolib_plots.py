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
import argparse
from matplotlib.ticker import AutoMinorLocator                                                                                                                                                                 
from matplotlib.colors import LogNorm                                                                                                                                                                          
import matplotlib                        
import json
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes                                                                                                                                                   
from mpl_toolkits.axes_grid1 import make_axes_locatable     
from pathlib import Path
import shutil
import subprocess

def get_title_coords(xmin, xmax, ymin, ymax):
    x_title=xmin+0.15*(xmax-xmin)
    y_title=ymin+0.2*(ymax-ymin)
    return x_title, y_title


# ==========================================================================
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
        if(quantity=="mass_function"):
            Y_min=1e-29
            Y_max=1e-1
            Y_max_max=2e-1
            file_name=self.data["output_file_mass_function"]
            index=1
            plt.ylabel(r"$n(M)$", fontsize=12)
            plt.title("Halo mass function")
            labels=[r"$n(M)$"]
            plt.yscale('log')
            xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))
            xtit=pow(10,xtit)
            ytit=pow(10,ytit)
        elif(quantity=="mass_dispersion"):
            Y_min=0
            Y_max=8.0
            Y_max_max=5.0
            file_name=self.data["output_file_mass_function"]
            index=2
            plt.ylabel(r"$\sigma(M)$", fontsize=12)
            plt.title("Mass dispersion")
            labels=[r"$\sigma(M)$"]
            plt.yscale('linear')
            xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), Y_min, Y_max)
            xtit=pow(10,xtit)
        elif(quantity=="halo_mass_bias"):
            Y_min=0
            Y_max=4.0
            Y_max_max=5.
            index=1
            file_name=self.data["output_file_halo_mass_bias"]
            plt.ylabel(r"$b(M)$", fontsize=12)
            plt.title("Halo-mass bias")
            labels=[r"$b(M)$"]
            plt.yscale('linear')
            xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), Y_min, Y_max)
            xtit=pow(10,xtit)

        plt.text(xtit, ytit, "z  = "+str(self.data["redshift"]), fontsize=11)
   
        tick_y=0.5
        min_locator_y=1
        aspect='auto'

        colors=["red", "blue", "green", "orange","cyan", "black", "red"]
        lines_s=["solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot","solid", "dotted", "dashed", "dashdot"]

       
        if os.path.isfile(file_name):
            data=np.genfromtxt(file_name)     
            x=(data[:,0])
            psc=(data[:,index])
            plt.plot(x,psc,color=colors[0],linestyle=lines_s[0],linewidth=3, alpha=0.6, label=labels[0])

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
        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))
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

    def plot_power_spectrum(self):
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
        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))
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
        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))

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
        xtit, ytit=get_title_coords(np.log10(X_min), np.log10(X_max), np.log10(Y_min), np.log10(Y_max))

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
    plotCosmo.plot_power_spectrum()


if plotCosmo.data.get("cosmology", False):
    print("Showing cosmology")
    plotCosmo.plot_cosmology_distances()

if plotCosmo.data.get("growth", False):
    print("Showing growth factor")
    plotCosmo.plot_cosmology_growth()


plt.show()