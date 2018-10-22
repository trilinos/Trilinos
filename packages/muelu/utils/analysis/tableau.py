#!/usr/bin/env python

def rgb2float(colors):
    """Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts"""
    for i in range(len(colors)):
        r, g, b = colors[i]
        colors[i] = (r / 255., g / 255., b / 255.)
    return colors

def tableau10():
    """'Tableau 10' colors as RGB"""
    colors = [
        ( 31, 119, 180), (255, 127,  14), ( 44, 160,  44), (214,  39,  40),
        (148, 103, 189), (140,  86,  75), (227, 119, 194), (127, 127, 127),
        (188, 189,  34), ( 23, 190, 207)
    ]
    return rgb2float(colors)

def tableau10_light():
    """'Tableau 10 Light' colors as RGB"""
    colors = [
        (174, 199, 232), (255, 187, 120), (152, 223, 138), (255, 152, 150),
        (197, 176, 213), (196, 156, 148), (247, 182, 210), (199, 199, 199),
        (219, 219, 141), (158, 218, 229)
    ]
    return rgb2float(colors)

def tableau10_medium():
    """'Tableau 10 Medium' colors as RGB"""
    colors = [
        (114, 158, 206), (255, 158,  74), (103, 191,  92), (237, 102,  93),
        (173, 139, 201), (168, 120, 110), (237, 151, 202), (162, 162, 162),
        (204, 204,  93), (109, 204, 218)
    ]
    return rgb2float(colors)

def tableau20():
    """'Tableau 20' colors as RGB"""
    colors = [
        ( 31, 119, 180), (174, 199, 232), (255, 127,  14), (255, 187, 120),
        ( 44, 160,  44), (152, 223, 138), (214,  39,  40), (255, 152, 150),
        (148, 103, 189), (197, 176, 213), (140,  86,  75), (196, 156, 148),
        (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
        (188, 189,  34), (219, 219, 141), ( 23, 190, 207), (158, 218, 229)
    ]
    return rgb2float(colors)
