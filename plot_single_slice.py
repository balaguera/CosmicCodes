import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.ndimage
from matplotlib.colors import LogNorm
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from scipy import ndimage
import matplotlib.colors as colors
import sys

plt.rc('text', usetex=True)  

if (len(sys.argv)!=8):
  sys.exit('number of arguments \n <name> <number of cells><number of adding slices><1,2,3 for coordinates><starting slice><Name><lenght of box-side>')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
#mpl.rcParams.update({'axes.titlesize': 'large'})
#mpl.rcParams.update({'axes.labelsize': 'large'})
#axis_font = {'size':'20'}

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

files = str(sys.argv[1])
nc = int(sys.argv[2])
L = int(sys.argv[7])

meanfiles = np.fromfile(files,dtype=np.float32)

#with open(files,'r') as f:
#      header = np.fromfile(f, dtype=np.float64, count=2)
#      header1 = np.fromfile(f, dtype=np.int32, count=1)
#      header2 = np.fromfile(f, dtype=np.float64, count=3)
#      meanfiles = np.fromfile(f, dtype=np.float64)
#      f.close()
          

b= meanfiles.reshape((nc,nc,nc),order='C')
bmean=np.mean(b)


#b=np.log10((b/bmean)+1)
#b=b/bmean-1
bmin=b.min()
bmax=b.max()

print (bmean, bmin, bmax)
start = int(sys.argv[5])
name = str(sys.argv[6])
num_slices = int(sys.argv[3]) 
end = start+num_slices
type = int(sys.argv[4])

if type==1: #slice at a fixd x coordinate
  c=b[:,:,start:end].mean(axis=-1) #([:,:,num]= y:z, x:y)  ([:,num,:]= y:z, x:x) ([num,:,:]= y:y, x:x) 
if type==2:
  c=b[:,start:end,:].mean(axis=1) #([:,:,num]= y:z, x:y)  ([:,num,:]= y:z, x:x) ([num,:,:]= y:y, x:x) 
if type==3:
  c=b[start:end,:,:].mean(axis=0) #([:,:,num]= y:z, x:y)  ([:,nuxm,:]= y:z, x:x) ([num,:,:]= y:y, x:x)

c = scipy.ndimage.interpolation.zoom(c ,order=1, zoom=2)

cmap = plt.get_cmap('terrain')
cmap = truncate_colormap(cmap, 0.15, .96)
cmap.set_under('w')

fig = plt.figure(figsize=(10, 8))


ax = fig.add_subplot(111)
minc=c.min()
maxc=c.max() 


im = plt.imshow(c,vmin=minc,vmax=maxc, interpolation='bilinear',aspect='auto', extent=[0,L,0,L],cmap=cmap) #DM

ax.set_aspect('equal')
if type==1:
    plt.xlabel(r'$y[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
    plt.ylabel(r'$z[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
if type==2:
    plt.xlabel(r'$x[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
    plt.ylabel(r'$z[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
if type==3:
    plt.xlabel(r'$x[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
    plt.ylabel(r'$y[h^{-1}\,\mathrm{Mpc}]$',fontsize=25)
plt.text(0.3*L, 1.03*L, name, fontsize=20, color='black')


plt.setp(ax.get_xticklabels(), fontsize=18)
#plt.setp(ax.get_yticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), fontsize=18)

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "5%", pad="2%")
cb = plt.colorbar(im,orientation='vertical',cax=cax)
cb.ax.get_yaxis().labelpad = 0.
#cb.set_label(r'$1+\delta$',size=20)
cb.ax.yaxis.set_tick_params(labelsize=18)




plt.show()
#fig.savefig("slice_IC_LOS980_Nres180_MAS1_original.pdf",dpi=200,bbox_inches='tight')
