"""Plot CAPE 1D and 3D CAPE-diff, relative error and histograms"""

import pandas
import numpy
import matplotlib.pyplot as plt
import scipy.interpolate

def force_aspect(axi,aspect):
    """ Function to set aspect ratio of figure """
    imag = axi.get_images()
    extent =  imag[0].get_extent()
    axi.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

# read and pre-calculate data
conv1d = pandas.read_csv('1d_convective_params.csv', header=[0])
conv1d = conv1d[['Lon', 'Lat', 'CAPE']]
conv1d.rename(columns={'CAPE': 'CAPE1D'}, inplace=True)

conv3d = pandas.read_csv('3d_convective_params.csv', header=[0])
conv3d = conv3d[['Lon', 'Lat', 'CAPE']]
conv3d.rename(columns={'CAPE': 'CAPE3D'}, inplace=True)

conv = pandas.merge(conv1d, conv3d)
conv['diff'] = conv['CAPE3D'] - conv['CAPE1D']

conv['rel_err'] = (conv['diff'].abs() / conv['CAPE1D']) * 100
conv['rel_err'] = conv['rel_err'].fillna(0)
conv['rel_err'].replace(numpy.inf, 100, inplace=True)

# plot histograms
conv['diff'].hist(bins=50, range=(-50,50)).get_figure().savefig("diffhist.png")
plt.clf()

conv['rel_err'].hist(bins=50, range=(0,10)).get_figure().savefig("errhist.png")
plt.clf()

# plot CAPE 1D
grid_x, grid_y = numpy.mgrid[14.1:23.9:500j, 49.1:54.9:500j]
coords = numpy.stack((conv["Lon"].values, conv["Lat"].values), axis=1)
values = conv["CAPE1D"].values
grid_z = scipy.interpolate.griddata(coords, values, (grid_x, grid_y), method='linear')

fig = plt.figure()
ax = fig.add_subplot(111)
img = plt.imshow(grid_z.T, extent=(14.25,23.75,49.25,54.75), origin='lower')
force_aspect(ax,aspect=1.0)
plt.colorbar()
fig.savefig("cape1d.png", dpi=300)
plt.clf()

# plot CAPE 3D
grid_x, grid_y = numpy.mgrid[14.1:23.9:500j, 49.1:54.9:500j]
coords = numpy.stack((conv["Lon"].values, conv["Lat"].values), axis=1)
values = conv["CAPE3D"].values
grid_z = scipy.interpolate.griddata(coords, values, (grid_x, grid_y), method='linear')

fig = plt.figure()
ax = fig.add_subplot(111)
img = plt.imshow(grid_z.T, extent=(14.25,23.75,49.25,54.75), origin='lower')
force_aspect(ax,aspect=1.0)
plt.colorbar()
fig.savefig("cape3d.png", dpi=300)
plt.clf()

# plot diff
grid_x, grid_y = numpy.mgrid[14.1:23.9:500j, 49.1:54.9:500j]
coords = numpy.stack((conv["Lon"].values, conv["Lat"].values), axis=1)
values = conv["diff"].values
grid_z = scipy.interpolate.griddata(coords, values, (grid_x, grid_y), method='linear')

fig = plt.figure()
ax = fig.add_subplot(111)
img = plt.imshow(grid_z.T, extent=(14.25,23.75,49.25,54.75), origin='lower', vmin=-25, vmax=25)
force_aspect(ax,aspect=1.0)
plt.colorbar()
fig.savefig("capediff.png", dpi=300)
plt.clf()

# plot err
grid_x, grid_y = numpy.mgrid[14.1:23.9:500j, 49.1:54.9:500j]
coords = numpy.stack((conv["Lon"].values, conv["Lat"].values), axis=1)
values = conv["rel_err"].values
grid_z = scipy.interpolate.griddata(coords, values, (grid_x, grid_y), method='linear')

fig = plt.figure()
ax = fig.add_subplot(111)
img = plt.imshow(grid_z.T, extent=(14.25,23.75,49.25,54.75), origin='lower', vmin=0, vmax=5)
force_aspect(ax,aspect=1.0)
plt.colorbar()
fig.savefig("capeerr.png", dpi=300)
plt.clf()
