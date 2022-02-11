import wradlib as wrl
import numpy


def get_radar_data(filename):
    radar_dict = wrl.io.read_rainbow(filename)["product"]["data"]["radarpicture"]

    datamin = float(radar_dict["@min"])
    datamax = float(radar_dict["@max"])
    datadepth = float(radar_dict["datamap"]["@depth"])

    rawdata = numpy.array(radar_dict["datamap"]["data"])
    data = datamin + (rawdata * ((datamax - datamin) / (2**datadepth)))

    return data
