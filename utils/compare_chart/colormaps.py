import numpy as np
import matplotlib.colors as mcolors
from pyparsing import col


def hex_to_rgb(value):
    """
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values"""
    value = value.strip("#")  # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i : i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    """
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values"""
    return [v / 256 for v in value]


def get_continuous_cmap(hex_list, float_list, alpha_list):
    """creates and returns a color map that can be used in heat map figures.
    If float_list is not provided, colour map graduates linearly between each color in hex_list.
    If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

    Parameters
    ----------
    hex_list: list of hex code strings
    float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

    Returns
    ----------
    colour map"""
    rgba_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]

    for i, col in enumerate(rgba_list):
        col.append(alpha_list[i])

    cdict = dict()
    for num, col in enumerate(["red", "green", "blue", "alpha"]):
        col_list = [
            [float_list[i], rgba_list[i][num], rgba_list[i][num]]
            for i in range(len(float_list))
        ]
        cdict[col] = col_list

    cmp = mcolors.LinearSegmentedColormap("my_cmp", segmentdata=cdict, N=256)

    return cmp


def get_radar_cmap():

    hex_list = [
        "#ffffff",
        "#0707c8",
        "#1261ff",
        "#00b4ff",
        "#3ccdff",
        "#78e6ff",
        "#bbf2ff",
        "#ffffff",
        "#fff799",
        "#ffee33",
        "#ffcb1a",
        "#ff7300",
        "#ff1e00",
        "#c80000",
        "#af0014",
        "#960028",
        "#be0078",
        "#d6119b",
        "#ee23be",
        "#f646d2",
        "#ff69e6",
        "#ffffff",
    ]

    values = list(
        np.array(
            [
                0.0,
                5.0,
                8.0,
                11.0,
                14.0,
                17.0,
                20.0,
                23.0,
                28.0,
                31.0,
                33.0,
                34.0,
                37.0,
                40.0,
                44.0,
                47.0,
                50.0,
                53.0,
                56.0,
                59.0,
                62.0,
                65.0,
            ]
        )
        / 65.0
    )

    alpha = [0.0] + [1.0] * 21

    colmap = get_continuous_cmap(hex_list, float_list=values, alpha_list=alpha)

    return colmap


def get_meteologix_cmap():
    hex_list = [
        "#1e7800",
        "#449900",
        "#68ac06",
        "#8cc00d",
        "#b1d414",
        "#d5e81b",
        "#fafc22",
        "#fad024",
        "#faa427",
        "#fb7929",
        "#fb4d2c",
        "#fc222f",
        "#fc2256",
        "#fc227e",
        "#fc22a6",
        "#fc22ce",
        "#fc22f6",
        "#fc41f7",
        "#fc61f8",
        "#fd80f9",
        "#fda0fb",
        "#febffc",
        "#fedffd",
    ]

    values = list(
        np.array(
            [
                0.0,
                20.0,
                40.0,
                100.0,
                200.0,
                300.0,
                400.0,
                600.0,
                800.0,
                1000.0,
                1200.0,
                1400.0,
                1600.0,
                1800.0,
                2000.0,
                2200.0,
                2400.0,
                2800.0,
                3200.0,
                3600.0,
                4000.0,
                4500.0,
                5000.0,
            ]
        )
        / 5000.0
    )

    alpha = [1.0] * 23

    colmap = get_continuous_cmap(hex_list, float_list=values, alpha_list=alpha)

    return colmap


def get_pivotal_cmap():
    hex_list = [
        "#ffffff",
        "#696969",
        "#37536a",
        "#a7c8ce",
        "#e9dd96",
        "#e16f02",
        "#dc4110",
        "#8b0950",
        "#73088a",
    ]

    values = list(
        np.array([0.0, 900.0, 1000.0, 1900.0, 2000.0, 2900.0, 3000.0, 3900.0, 4000.0])
        / 4000.0
    )

    alpha = [1.0] * 9

    colmap = get_continuous_cmap(hex_list, float_list=values, alpha_list=alpha)

    return colmap
