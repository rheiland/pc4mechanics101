#
# create_subcells_matrix.py - parse an annotated .svg file containing cells of interest and extract 
#        polyline approximations to the cells' boundaries.
#        Fill the boundaries with a dense packing of PhysiCell (sub)cells and create a .csv file.
#        Also, read in previously saved "matrix" agents and plot.
#
#
# Author: Randy Heiland
#
__author__ = "Randy Heiland"

import sys
import os
import xml.etree.ElementTree as ET
import math
from svg.path import parse_path, Path, Move
from scipy.spatial import ConvexHull
import itertools

try:
  import matplotlib
  from matplotlib.patches import Circle, Ellipse, Rectangle
  from matplotlib.collections import PatchCollection
except:
  print("\n---Error: cannot import matplotlib")
  print("---Try: python -m pip install matplotlib")
  print(join_our_list)
#  print("---Consider installing Anaconda's Python 3 distribution.\n")
  raise
try:
  import numpy as np  # if mpl was installed, numpy should have been too.
except:
  print("\n---Error: cannot import numpy")
  print("---Try: python -m pip install numpy\n")
  print(join_our_list)
  raise

import matplotlib.pyplot as plt
from collections import deque
from operator import add


print("# args=",len(sys.argv)-1)
fname = "rwh1b.svg"
fname = sys.argv[1]


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

#--------------------------------

fig = plt.figure(figsize=(7,7))
ax = fig.gca()
#ax.set_aspect("equal")
#plt.ion()

xlist = deque()
ylist = deque()
rlist = deque()
rgb_list = deque()

print('\n---- ' + fname + ':')
tree = ET.parse(fname)
root = tree.getroot()
#  print('--- root.tag ---')
#  print(root.tag)
#  print('--- root.attrib ---')
#  print(root.attrib)

#  print('--- child.tag, child.attrib ---')
numChildren = 0
bx = [0,1,2,3]
by = [0,2,1,4]
#plt.plot(bx,by)
by2 = list(map(add,by,[1,1,1,1]))
#plt.plot(bx,by2)
#plt.show()

xoff = -300
yoff = 250
scale_factor = 5.0

plot_pieces_flag = False

def perp( a ):
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def seg_intersect(a1,a2, b1,b2):
    # print("seg_intersect: a1,a2,b1,b2=", a1,a2,b1,b2)
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom)*db + b1

    # if abs(denom) < 1.e-6:
    #     # print("denom=",denom,"  < 1.e-6")
    #     # sys.exit()
    #     return np.array([999,999])
    # else:
    #     return (num / denom)*db + b1

#-------------------
# cell_radius = 8.412710547954228  # PhysiCell default
cell_radius = 1.  # ~2 micron spacing
cell_radius = 2.5 # ~5 micron spacing of subcells
cell_diam = cell_radius*2

#yc = -1.0
y_idx = -1
# hex packing constants
x_spacing = cell_radius*2
y_spacing = cell_radius*np.sqrt(3)

guess_cell_max_x = 50.0  # somewhat arbitrary

subcells_flag = False
mapped_geom_flag = False
subcells_flag2 = True
# subcells_flag2 = False
convexhull_flag = True
debug_print = False

my_data_name = ''
# proximal tubule (S1) 0, distal tubule 1, macula densa 2, endothelial cells 3, parietal endothelial cells 4, 
# juxtaglomerular granular cell 5, vascular smooth muscle cell 7
#id_color = ['Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Indigo', 'Violet']
#id_color = ['Red', 'Cyan', 'Blue', 'Pink', 'springgreen', 'yellow', 'red']
id_color = ['red', 'skyblue', 'dodgerblue', 'lightpink', 'palegreen', 'red', 'green','green']
#id_name_color_dict = ['proximal_tubule', 'proximal_tubule_S1_', 'basement_membranes', 'parietal_basement_membrane', 'glomerular_basement_membrane', 'distal_tubule', 'macula_densa', 'endothelial_cells', 'parietal_epithelial_cells', 'juxtaglomerular_granular_cell', 'renin_granules', 'vascular_smooth_muscle_cell', 'extra_labels']


# from https://github.com/PhysiCell-Models/Kidney_FTU/blob/main/PhysiCell/config/mymodel.xml 
    #   <cell_definition name="endothelial" ID="0">
    #   <cell_definition name="vascular_smooth_muscle" ID="1">
    #   <cell_definition name="parietal_epithelial" ID="2">
    #   <cell_definition name="podocyte" ID="3">
    #   <cell_definition name="parietal_basement_membrane" ID="4">
    #   <cell_definition name="proximal_tube_epithelial" ID="5">
    #   <cell_definition name="glomerular_basement_membrane" ID="6">
    #   <cell_definition name="glomerular_capillary_endothelial" ID="7">
    #   <cell_definition name="glomerular_mesangial" ID="8">
    #   <cell_definition name="extraglomerular_mesangial" ID="9">
    #   <cell_definition name="macula_densa" ID="10">
    #   <cell_definition name="juxtaglomerular_granule" ID="11">
    #   <cell_definition name="mesangial_matrix" ID="12">
    #   <cell_definition name="capillary" ID="13">

# cf. results of this script: ~/git/Kidney_Importer/python$ grep data-name foo.out
# data-name =  proximal tubule
# data-name =  proximal tubule (S1)
# data-name =  distal tubule
# data-name =  macula densa
# data-name =  endothelial cells
# data-name =  parietal epithelial cells
# data-name =  juxtaglomerular granular cell
# data-name =  vascular smooth muscle cell
# data-name =  extra labels

cell_type_dict = {'proximal':5, 'macula':10, 'distal':10, 'endothelial':0, 'endo-blib':0, 'parietal':4, 'juxtag':11, 'vascular':1, 'podocyte':3, 'capillary':7}


id_count = -2  # for some weird reason, use -2 instead of -1 (it's being incremented twice somehow before being used, maybe proximal_tubule)
pec_flag = False
mycolor = "black"
cell_id = -1
# cells_file = 'ftu_cells.csv'
# filep = open(cells_file, 'w')  # 'w', 'a' to append?
# filep.close()

def parse_children(parent):
    global id_count,pec_flag,id_color,my_data_name, mycolor, cell_id,filep

    path_count = 0
    for child in parent:
        # print(child.tag, child.attrib)
        # print("keys = ",child.attrib.keys())
        if 'id' in child.attrib.keys():
            if debug_print:
                print("--> id = ",child.attrib['id'])
            if "parietal" in child.attrib['id']:  # just for matrix agents testing
                continue
            if "basement_membranes" in child.attrib['id']:
                continue
            if "Bowman" in child.attrib['id']:
                continue
            if "extraglomerular" in child.attrib['id']:
                continue
            if "renin_granules" in child.attrib['id']:
                continue
            if "macula" in child.attrib['id']: # let's not include the "wheel" of cells on left
                continue
            if "distal" in child.attrib['id']: # let's not include the "wheel" of cells on left
                continue
            # if "podocytes" in child.attrib['id']:
            if "podocyte" in child.attrib['id']:
                mycolor = "lightsteelblue"
                # my_data_name = 'podocyte'
                my_data_name = child.attrib['id']
            if "endo-blib" in child.attrib['id']:
                mycolor = "cyan"
                my_data_name = child.attrib['id']
            # if "juxtaglomerular_granular_cell" in child.attrib['id']:
                # continue
            # if "extraglomerular_mesangium" in child.attrib['id']:
            # if "mesangium" == child.attrib['id']:  # avoid large rectangular blocks on left
                # continue
        if 'data-name' in child.attrib.keys():
            my_data_name = child.attrib['data-name']
            # print("keys = ",child.attrib.keys())
            print("data-name = ",child.attrib['data-name'])
            if "parietal epithelial cells" in child.attrib['data-name']:
                pec_flag = True
            else:
                pec_flag = False

            if 'proximal_tubule' in child.attrib['id']:
                # pass
                id_count += 1
            else:
                id_count += 1

        if "path" in child.tag:
        # if pec_flag and "path" in child.tag:
            # print("found 'path' in child.tag = ",child.tag,"\n")
            path_count += 1
            cell_id += 1
            # print("found 'path' in child.tag = ",child.tag,", cell_id=",cell_id,"\n")
            if debug_print:
                print("found 'path' in child.tag = ",child.tag,", path_count = ",path_count,", cell_id=",cell_id,"\n")
            d_str = child.attrib['d']
            # print("\n--- d_str = ",d_str,"\n")
            cpath = parse_path(d_str)   # from svg.path
            # print(cpath)
            if debug_print:
                print("-> id (child)= ",child.attrib['id'])
            xv = np.array([])
            yv = np.array([])
            xpts_path = np.array([])
            ypts_path = np.array([])
            # xy = np.array([],[])
            xy = np.array([[],[]], np.float64)
            for idx in range(0,len(cpath)):
                # print("--> ",cpath[idx])
                if 'CubicBezier' in str(cpath[idx]): # NB: "CubicBezier" is NOT in the .svg; it's from svg.path:parse_path
                    if debug_print:
                        print("  --> id (child with bezier) = ",child.attrib['id'])
                    cbez = cpath[idx]
                    # print("--> ",cpath[idx])
                    # bx = [cbez.point(0).real, cbez.point(1).real, cbez.point(2).real, cbez.point(3).real ]
                    # by = [cbez.point(0).imag, cbez.point(1).imag, cbez.point(2).imag, cbez.point(3).imag ]
                    bx = [cbez.start.real, cbez.control1.real, cbez.control2.real, cbez.end.real ]
                    nbx = np.array(bx)
                    # xv = np.append(xv,[cbez.start.real, cbez.control1.real, cbez.control2.real, cbez.end.real ])
                    by = [cbez.start.imag, cbez.control1.imag, cbez.control2.imag, cbez.end.imag ]
                    nby = np.array(by)
                    nby *= -1
                    # yv = np.append(yv,[cbez.start.imag, cbez.control1.imag, cbez.control2.imag, cbez.end.imag ])
                    # yv *= -1
                    # print("nbx= ",nbx)
                    # print("nby= ",nby)
                    nbx += xoff
                    nby += yoff
                    nbx *= scale_factor
                    nby *= scale_factor
                    xpts_path = np.append(xpts_path, nbx)
                    ypts_path = np.append(ypts_path, nby)
                    if plot_pieces_flag:
                        plt.plot(nbx,nby)
                    # plt.plot(xv,yv)  # wth?
                    # xy_pts = np.array(list(itertools.product(xv,yv)))
                    # xy_pts = np.array(list(itertools.product(nbx,nby)))
                    # plt.plot(xy_pts[:,0], xy_pts[:,1], 'o')
                    # plt.plot(xy_pts[:,0], xy_pts[:,1])
                    # break
                elif 'Line' in str(cpath[idx]):
                    if debug_print:
                        print("\n------ found: ",str(cpath[idx]))
                    line = cpath[idx]
                    bx = [line.start.real, line.end.real]
                    nbx = np.array(bx)
                    by = [line.start.imag, line.end.imag]
                    nby = np.array(by)
                    nby *= -1
                    nbx += xoff
                    nby += yoff
                    nbx *= scale_factor
                    nby *= scale_factor
                    xpts_path = np.append(xpts_path, nbx)
                    ypts_path = np.append(ypts_path, nby)
                    if plot_pieces_flag:
                        plt.plot(nbx,nby)
                elif 'Arc' in str(cpath[idx]):
                    if debug_print:
                        print("\n------ found: ",str(cpath[idx]))
                    arc = cpath[idx]
                    bx = [arc.start.real, arc.end.real]
                    nbx = np.array(bx)
                    by = [arc.start.imag, arc.end.imag]
                    nby = np.array(by)
                    nby *= -1   # invert in y for svg
                    nbx += xoff
                    nby += yoff
                    nbx *= scale_factor
                    nby *= scale_factor
                    xpts_path = np.append(xpts_path, nbx)
                    ypts_path = np.append(ypts_path, nby)
                    if plot_pieces_flag:
                        plt.plot(nbx,nby)

        # if "path" in child.tag:
            # for idx in range(0,len(cpath)):
            # obtain the slope (best fit line) of the cell's boundary points.
            if len(xpts_path) == 0:
                print("=========> break, due to len(xpts_path) = 0")
                break
            m, b = np.polyfit(xpts_path, ypts_path, 1)
            theta_rad = np.arctan(m)
            # print("m, theta = ",m,theta_rad)
            # theta_rad = -(np.pi - theta_rad)
            theta_rad += np.pi/2.0  # make vertical

            # translate the points to the origin.
            x_center = xpts_path.min() + (xpts_path.max() - xpts_path.min())/2.0
            y_center = ypts_path.min() + (ypts_path.max() - ypts_path.min())/2.0
            xpts_path -= x_center
            ypts_path -= y_center
            
            # rotate the points (centered at the origin) to lie in the ~vertical polyline boundary.
            xpts_new = xpts_path * np.cos(theta_rad) + ypts_path * np.sin(theta_rad) 
            ypts_new = xpts_path * -np.sin(theta_rad) + ypts_path * np.cos(theta_rad) 
            y_min = ypts_new.min()
            y_max = ypts_new.max()
            # print("y_min, y_max = ",y_min, y_max)

            # Find the convex hull of these transformed points
            # xy_hull = np.array(list(itertools.product(xpts_new,ypts_new)))
            # hull = ConvexHull(xy_hull)
            # # print("hull = ",hull)  # <scipy.spatial.qhull.ConvexHull object ...>
            # if path_count % 2 == 0:  # Warning: Assumes that the (smaller) nucleus <path> follows the (larger) cell <path>. We skip over it.
            #     hull_nucleus_x = xy_hull[hull.vertices,0]
            #     hull_nucleus_y = xy_hull[hull.vertices,1]
            #     # print("hull_nucleus_x = ",hull_nucleus_x)
            #     # print("hull_nucleus_y = ",hull_nucleus_y)
            #     # Does the nucleus lie within the cell? (approx)
            #     if (hull_nucleus_x.min() > hull_cell_x.min()) and (hull_nucleus_x.max() < hull_cell_x.max()) and \
            #        (hull_nucleus_y.min() > hull_cell_y.min()) and (hull_nucleus_y.max() < hull_cell_y.max()):
            #         continue
            # else:
            #     hull_cell_x = xy_hull[hull.vertices,0]
            #     hull_cell_y = xy_hull[hull.vertices,1]
            #     # print("hull_cell_x = ",hull_cell_x)
            #     # print("hull_cell_y = ",hull_cell_y)

            # # plt.plot(xy_hull[hull.vertices,0], xy_hull[hull.vertices,1], 'r--', lw=2)
            # if convexhull_flag:
            #     plt.plot(xy_hull[hull.vertices,0], xy_hull[hull.vertices,1], 'g', lw=2)

            if mapped_geom_flag:
                plt.plot(xpts_new, ypts_new, 'ro')
                plt.plot(xpts_new, ypts_new)

            # plt.title(child.attrib['id'])
            # plt.title("subcells for parietal epithelial cells")
            # plt.title("cell w/ nucleus; convex hulls")

            hline0 = np.array( [-guess_cell_max_x, 0.0] )  # artificially chosen x range for horiz line.
            hline1 = np.array( [guess_cell_max_x, 0.0] )
            lseg_p0 = np.array( [0.0, 0.0] )
            lseg_p1 = np.array( [0.0, 0.0] )

            cells_x = np.array([])
            cells_y = np.array([])

            if debug_print:
                print("\n\n-- bruteforce check all polyline line segments:")
            y_idx = 0
            # for yval in np.arange(y_min,y_max, cell_diam):
            for yval in np.arange(y_min,y_max, y_spacing):
                xvals = []
                y_idx += 1
                # print("\n--- yval, len(xpts_new) = ",yval, len(xpts_new))
                hline0[1] = yval
                hline1[1] = yval
                for idx in range(len(xpts_new)-1):
                    # print(idx,") ",xpts_new[idx],ypts_new[idx], " -> ", xpts_new[idx+1],ypts_new[idx+1])
                    lseg_p0[0] = xpts_new[idx]
                    lseg_p0[1] = ypts_new[idx]
                    lseg_p1[0] = xpts_new[idx+1]
                    lseg_p1[1] = ypts_new[idx+1]

                    ptint = seg_intersect( hline0,hline1, lseg_p0,lseg_p1)
                    # if ptint[0] >= 999:
                        # continue
                    xmin = min(lseg_p0[0], lseg_p1[0])
                    xmax = max(lseg_p0[0], lseg_p1[0])
                    if ptint[0] >= xmin and ptint[0] <= xmax:
                        # print("------------ ptint = ",ptint)
                        # print("------------ xmin,xmax = ",xmin,xmax)
                        xvals.append(ptint[0])
                        # print("--> ",ptint[0],ptint[1])
                    # else:
                        # print("-- no intersection.")

                lseg_p0[0] = xpts_new[idx+1]
                lseg_p0[1] = ypts_new[idx+1]
                lseg_p1[0] = xpts_new[0]
                lseg_p1[1] = ypts_new[0]
                # print("yval=",yval,", idx=", idx,") last -> 0th: ",lseg_p0[0],lseg_p0[1], " -> ", lseg_p1[0],lseg_p1[1])
                ptint = seg_intersect( hline0,hline1, lseg_p0,lseg_p1)
                xmin = min(lseg_p0[0], lseg_p1[0])
                xmax = max(lseg_p0[0], lseg_p1[0])
                if ptint[0] >= xmin and ptint[0] <= xmax:
                    xvals.append(ptint[0])
                    # print("--> ",ptint[0],ptint[1])

                # print("(presorted) xvals = ",end='')
                # for kdx in range(len(xvals)):
                #     print(xvals[kdx],',',end='')
                # print()

                xvals.sort()
                # print("(sorted) xvals = ",end='')
                # for kdx in range(len(xvals)):
                #     print(xvals[kdx],',',end='')
                # print()

                if len(xvals) == 1:
                    pass
                else:
                    for xval in np.arange(-guess_cell_max_x,guess_cell_max_x, x_spacing):
                        xval_offset = xval + (y_idx%2) * cell_radius
                        for kdx in range(0,len(xvals),2):
                            if (xval >= xvals[kdx]) and (xval <= xvals[kdx+1]):
                                cells_x = np.append(cells_x, xval_offset)
                                cells_y = np.append(cells_y, yval)

                # plt.plot(cells_x,cells_y,'go')
                # circles(cells_x,cells_y, s=cell_radius, color=rgbs, alpha=0.6, ed='black', linewidth=0.5)

            if subcells_flag:
                # circles(cells_x,cells_y, s=cell_radius, ec='black', linewidth=0.1)
                # if id_count < len(id_color):
                    # mycolor = id_color[id_count]
                # else:
                    # mycolor = 'b'
                # mycolor = id_color[id_count % len(id_color)]
                # print(" ----------->>>>>>>>>>>  mycolor = ",mycolor)
                circles(cells_x,cells_y, s=cell_radius, c='b', ec='black', linewidth=0.1)

            # Now do the inverse transformations to put back in the original cell's location

            # rotate the points back to their original
            theta_rad = -theta_rad
            xpts_new2 = cells_x * np.cos(theta_rad) + cells_y * np.sin(theta_rad) 
            ypts_new2 = cells_x * -np.sin(theta_rad) + cells_y * np.cos(theta_rad) 
            # y_min = ypts_new.min()
            # y_max = ypts_new.max()
            # # print("y_min, y_max = ",y_min, y_max)

            # # translate the points away the origin.
            xpts_new2 += x_center
            ypts_new2 += y_center
                
            if subcells_flag2:
                # if id_count < len(id_color):
                #     mycolor = id_color[id_count]
                # else:
                #     mycolor = 'b'
                if mycolor == "lightsteelblue":
                    pass
                else:
                    mycolor = id_color[id_count % len(id_color)]
                # print("id_count, len(id_color) = ", id_count, len(id_color) )
                # print(" ----------->>  my_data_name, id_count, mycolor = ",my_data_name, id_count, mycolor)
                if 'podocyte' in my_data_name:
                    tmp_cell_id = cell_id
                    if '1' in my_data_name:
                        cell_id = 901
                    elif '2' in my_data_name:
                        cell_id = 902
                    elif '3' in my_data_name:
                        cell_id = 903
                    elif '4' in my_data_name:
                        cell_id = 904
                    elif '5' in my_data_name:
                        cell_id = 905
                if 'endo-blib' in my_data_name:
                    tmp_cell_id = cell_id
                    if '1' in my_data_name:
                        cell_id = 801
                    elif '2' in my_data_name:
                        cell_id = 802
                    elif '3' in my_data_name:
                        cell_id = 803
                    elif '4' in my_data_name:
                        cell_id = 804
                print(" ----------->>  my_data_name, id_count, mycolor, # pts, cell_id = ",my_data_name, id_count, mycolor, len(xpts_new2),cell_id)
                # circles(xpts_new2,ypts_new2, s=cell_radius, ec='k', linewidth=0.1)
                # circles(xpts_new2,ypts_new2, c=mycolor, s=cell_radius, ec='k', linewidth=0.1)
                circles(xpts_new2,ypts_new2, c=mycolor, s=cell_radius, linewidth=0.1)

                found_k = False
                for k in cell_type_dict.keys():
                    if k in my_data_name:
                        # print(" ~~~~~   cell_id = ",cell_id)
                        found_k = True
                        # f = open(k+'.csv', 'a')  # 'w', 'a' to append?
                        # filep = open(cells_file, 'a')  # 'w', 'a' to append?
                        # if (len(xpts_new2) == 0):
                        #     print("----- skipping over 0-length subcells!")
                        # for ipt in range(len(xpts_new2)):
                        #     filep.write(f"{xpts_new2[ipt]},{ypts_new2[ipt]}, 0.0, {cell_type_dict[k]}, {cell_id}\n")
                        # filep.close()
                if not found_k:
                    print("\n ======= oops, ",k," not in ",my_data_name)
                        # f.close()
                if 'podocyte' in my_data_name: # restore saved 
                    cell_id = tmp_cell_id

                # if (id_count >= 1):
                    # break


        parse_children(child)

parse_children(root)

# filep.close()
# print("\n-------> ",cells_file)

# xv_mtx = np.load("../data/matrix_x.npy")
# xv_mtx *= 0.98
# xv_mtx -= 140
# yv_mtx = np.load("../data/matrix_y.npy")
# yv_mtx *= 0.55
# yv_mtx += 20
# circles(xv_mtx,yv_mtx, c='black', s=cell_radius, linewidth=0.1)

# match cell type in .xml: <cell_definition name="mesangial_matrix" ID="12">
# cells_file = "mtx_cells2.csv"
# filep = open(cells_file, 'w')
# matrix_cell_type = 12
# for ipt in range(len(xv_mtx)):
#     filep.write(f"{xv_mtx[ipt]},{yv_mtx[ipt]}, 0.0, {matrix_cell_type}, {601}\n")
# filep.close()
# print("\n-------> ",cells_file)

plt.show()
