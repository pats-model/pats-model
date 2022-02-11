import pandas
import numpy
import scipy.interpolate
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.io import netcdf_file
from colormaps import get_radar_cmap, get_pivotal_cmap
from radar_reader import get_radar_data
from carto_plot import configure_carto_plot


def main():
    # read pre-calculate and prepare data
    conv1d = pandas.read_csv("1d_convective_params.csv", header=[0])
    conv1d = conv1d[["Lon", "Lat", "CAPE"]]
    conv1d.rename(columns={"CAPE": "CAPE1D"}, inplace=True)

    conv3d = pandas.read_csv("3d_convective_params.csv", header=[0])
    conv3d = conv3d[["Lon", "Lat", "CAPE"]]
    conv3d.rename(columns={"CAPE": "CAPE3D"}, inplace=True)

    conv = pandas.merge(conv1d, conv3d)

    # compute data
    grid_x, grid_y = numpy.mgrid[14.1:23.9:491j, 49.1:54.9:291j]
    coords = numpy.stack((conv["Lon"].values, conv["Lat"].values), axis=1)

    grid_z_1d = scipy.interpolate.griddata(
        coords, conv["CAPE1D"].values, (grid_x, grid_y), method="linear"
    )

    grid_z_3d = scipy.interpolate.griddata(
        coords, conv["CAPE3D"].values, (grid_x, grid_y), method="linear"
    )

    # read ERA5 data
    dataset = netcdf_file("era5_cape.nc", maskandscale=True, mmap=False)
    cape_era = dataset.variables["cape"][0, :, :]
    lats = dataset.variables["latitude"][:]
    lons = dataset.variables["longitude"][:]

    # read radar data
    radar_data = get_radar_data("imgw_radar.cmax")

    # plot
    fig = plt.figure(figsize=(12, 12), dpi=300)

    ##########
    ax0 = fig.add_subplot(2, 2, 3, projection=ccrs.Mercator())
    ax0.title.set_text("CAPE ERA5")

    ax0 = configure_carto_plot(ax0)

    plt.contourf(
        lons,
        lats,
        cape_era,
        numpy.linspace(0.0, 4000.0, 101),
        cmap=get_pivotal_cmap(),
        extend='max',
        transform=ccrs.PlateCarree(),
    )

    ############
    ax1 = fig.add_subplot(2, 2, 1, projection=ccrs.Mercator())
    ax1.title.set_text("CAPE 1D")

    ax1 = configure_carto_plot(ax1)

    ax1.pcolormesh(
        numpy.arange(14.09, 23.91, 0.02),
        numpy.arange(49.09, 54.91, 0.02),
        grid_z_1d.T,
        vmin=0.0,
        vmax=4000.0,
        cmap=get_pivotal_cmap(),
        transform=ccrs.PlateCarree(),
    )

    ###############
    ax2 = fig.add_subplot(2, 2, 2, projection=ccrs.Mercator())
    ax2.title.set_text("CAPE 3D")

    ax2 = configure_carto_plot(ax2)
    ax2.yaxis.tick_right()

    im = ax2.pcolormesh(
        numpy.arange(14.09, 23.91, 0.02),
        numpy.arange(49.09, 54.91, 0.02),
        grid_z_3d.T,
        vmin=0.0,
        vmax=4000.0,
        cmap=get_pivotal_cmap(),
        transform=ccrs.PlateCarree(),
    )

    ###############
    ax3 = fig.add_subplot(2, 2, 4, projection=ccrs.Mercator())
    ax3.title.set_text("RADAR")

    ax3 = configure_carto_plot(ax3)
    ax3.yaxis.tick_right()

    ax3.pcolormesh(
        numpy.arange(-450000.0, 450000.0, 1000.0),
        numpy.arange(-450000.0, 450000, 1000.0),
        radar_data,
        cmap=get_radar_cmap(),
        vmin=0.0,
        vmax=65.0,
        transform=ccrs.AzimuthalEquidistant(
            central_latitude=52.3469, central_longitude=19.0926
        ),
    )

    # add colorbar
    # position of the <given> edge of the subplots, as a fraction of the figure height
    fig.subplots_adjust(top=0.85)
    # [left, bottom, width, height]
    cbar_ax = fig.add_axes([0.1, 0.9, 0.8, 0.05])
    fig.colorbar(im, cax=cbar_ax, orientation="horizontal")

    print("Saving...")

    fig.savefig("cape_maps.png", bbox_inches="tight")


if __name__ == "__main__":
    main()
