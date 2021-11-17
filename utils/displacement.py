"""Plot CAPE 1D and 3D CAPE-diff, relative error and histograms"""

import pandas
import numpy
import matplotlib.pyplot as plt
import scipy.interpolate


def force_aspect(axi, aspect):
    """ Function to set aspect ratio of figure """
    imag = axi.get_images()
    extent = imag[0].get_extent()
    axi.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


# read and pre-calculate data
conv3d = pandas.read_csv('3d_convective_params.csv', header=[0])
conv3d = conv3d[['Lon', 'Lat', 'xHorizontalDisplacement',
                 'yHorizontalDisplacement']]

# plot histograms
conv3d['xHorizontalDisplacement'].hist(
    bins=50).get_figure().savefig("xdishist.png")
plt.clf()

conv3d['yHorizontalDisplacement'].hist(
    bins=50).get_figure().savefig("ydishist.png")
plt.clf()

# plot x displacement
grid_x, grid_y = numpy.mgrid[14.1:23.9:500j, 49.1:54.9:500j]
coords = numpy.stack((conv3d["Lon"].values, conv3d["Lat"].values), axis=1)
values = conv3d["xHorizontalDisplacement"].values
grid_z = scipy.interpolate.griddata(
    coords, values, (grid_x, grid_y), method='linear')

fig = plt.figure()
ax = fig.add_subplot(111)
scale_max = conv3d["xHorizontalDisplacement"].quantile(0.99)
scale_min = conv3d["xHorizontalDisplacement"].quantile(0.01)
img = plt.imshow(grid_z.T, extent=(14.25, 23.75, 49.25, 54.75),
                 origin='lower', vmin=scale_min, vmax=scale_max)
force_aspect(ax, aspect=1.0)
plt.colorbar()
fig.savefig("xdisplac.png", dpi=300)
plt.clf()

# plot y displacement
grid_x, grid_y = numpy.mgrid[14.1:23.9:500j, 49.1:54.9:500j]
coords = numpy.stack((conv3d["Lon"].values, conv3d["Lat"].values), axis=1)
values = conv3d["yHorizontalDisplacement"].values
grid_z = scipy.interpolate.griddata(
    coords, values, (grid_x, grid_y), method='linear')

fig = plt.figure()
ax = fig.add_subplot(111)
scale_max = conv3d["yHorizontalDisplacement"].quantile(0.99)
scale_min = conv3d["yHorizontalDisplacement"].quantile(0.01)
img = plt.imshow(grid_z.T, extent=(14.25, 23.75, 49.25, 54.75),
                 origin='lower', vmin=scale_min, vmax=scale_max)
force_aspect(ax, aspect=1.0)
plt.colorbar()
fig.savefig("ydisplac.png", dpi=300)
plt.clf()
