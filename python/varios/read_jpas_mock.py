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
# ============================================================================
# ============================================================================
    
dir_path = "/home/andres/data/Numerics/JPAS/MOCKS/"

#filename=dir_path+"galcat_v0p1_pinocchio_miriam040_ID.hdf5"
#ofilename=dir_path+"galcat_v0p1_pinocchio_miriam040_ID.txt"
filename=dir_path+"galaxy_catalogue_pinocchio_v0.2_Miriam040_AllMass.hdf5"
ofilename=dir_path+"galaxy_catalogue_pinocchio_v0.2_Miriam040_AllMass.txt_b"
outfile = open(ofilename, 'w') 



"""
delimchar = '\t'
colnames = ['/ra', '/dec', '/zobs', '/app_mag', '/col',
            '/halo_mass', '/log_stell_mass']
fmts = ['%f', '%f', '%f', '%f', '%f', '%f', '%f']
np.savetxt(ofilename, h5py.File(filename)['/home/andres/data/Numerics/JPAS/MOCKS/galaxy_catalogue_pinocchio_v0.2_Miriam040_AllMass'],
           fmt=fmts, delimiter=delimchar, header=delimchar.join(colnames))

"""

with h5py.File(filename, 'r') as f:
    print("Reading from file", filename)
    for dset in traverse_datasets(filename):
        print('Path:', dset)
        print('Shape:', f[dset].shape)
        print('Data type:', f[dset].dtype)

        n_parts=len(f[dset])
    print("Writting")
    for i  in range(21216874, n_parts):
        outfile.write(str(round(f["/ra"][i],6))+" "+str(round(f["/dec"][i],6))
                      +" "+str(round(f["/zobs"][i],6))+" "+str(round(f["/app_mag"][i],6))+" "+str(round(f["/col"][i],6))+" "+str(round(f["/halo_mass"][i],6))+" "+str(round(f["/log_stell_mass"][i],6)))
        outfile.write("\n")     
    print("Done")




"""
    centrals=[0 for i in range(n_parts)]
    print("Assigning c/s")

    for i  in range(n_parts):
        if(f["/is_cen"][i]==True):
            centrals[i]
 
    print("Done")

    for i  in range(n_parts):
        for dset in traverse_datasets(filename):
             if(f["/is_cen"][i]==True):
                 outfile.write(str(1)+" ")
             elif (f[dset][i]==False)  :  
                outfile.write(str(0)+" ")
             else:
                outfile.write(str(round(f[dset][i],6))+" ")
                 
        outfile.write("\n")     

"""
    




