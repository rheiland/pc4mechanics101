import sys
import os
import glob
import numpy as np
from pyMCDS_cells import pyMCDS_cells
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle
from matplotlib.collections import PatchCollection

argc = len(sys.argv)-1
print("# args=",argc)

#data_dir = 'output'
if (argc < 2):
#  data_dir = int(sys.argv[kdx])
  print("Usage: provide output subdir and idx")
  sys.exit(-1)

kdx = 1
data_dir = sys.argv[kdx]
print('data_dir = ',data_dir)
kdx += 1
iframe = int(sys.argv[kdx])
print('iframe = ',iframe)

os.chdir(data_dir)
xml_files = glob.glob('output*.xml')
os.chdir('..')
xml_files.sort()
#print('xml_files = ',xml_files)

ds_count = len(xml_files)
print("----- ds_count = ",ds_count)
mcds = [pyMCDS_cells(xml_files[i], data_dir) for i in range(ds_count)]

tval = np.linspace(0, mcds[-1].get_time(), ds_count)
print('tval= ',tval)

#-----------------------------------------------------
def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    See https://gist.github.com/syrte/592a062c562cd2a98a83 

    Make a scatter plot of circles. 
    Similar to plt.scatter, but the size of circles are in data scale.
    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.
    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`
    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    plt.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.

    zipped = np.broadcast(x, y, s)
    patches = [Circle((x_, y_), s_)
               for x_, y_, s_ in zipped]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        c = np.broadcast_to(c, zipped.shape).ravel()
        collection.set_array(c)
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    plt.draw_if_interactive()
    if c is not None:
        plt.sci(collection)
    return collection

#-----------------------------------
# default cell type
#yval1 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 0) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100.) == True)) for idx in range(ds_count)] )

# plt.plot(tval, yval1, label='Mac', linewidth=1, color='red')
#plt.legend(loc='center left', prop={'size': 15})
#plt.legend(loc='upper left', prop={'size': 10})

#idx = 0
xvals = mcds[iframe].data['discrete_cells']['position_x']
yvals = mcds[iframe].data['discrete_cells']['position_y']
total_vol = mcds[iframe].data['discrete_cells']['total_volume']
r = total_vol
# print(r)
r = r / 3.1415
# print(r)
r = np.sqrt(r)
# print(r)
# plt.plot(xvals,yvals,'o')
#circles(xvals,yvals, s=rvals, color=rgbs)
circles(xvals,yvals, s=r)

plt.title(data_dir + ": frame= " + str(iframe))
#plt.savefig(data_dir + '.png')
plt.show()
