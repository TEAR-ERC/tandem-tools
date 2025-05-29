#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:06:21 2024

@author: bar
"""
'''
Vths: threshold velocity defining coseismic phase
Vp: plate velocity
V0: reference velocity
vmin,vmax: min./max. velocity for the colormap
Based on Jeena's code!
'''
import matplotlib as mpl
import numpy as np

#%%



def ReturnCmap(vmin, vmax, Vths=1e-4, Vp=1e-9, V0=1e-6):
    '''
    Vths: threshold velocity defining coseismic phase
    Vp: plate velocity
    V0: reference velocity
    vmin, vmax: min./max. velocity for the colormap
    '''

    if vmin >= vmax:
        raise ValueError("vmin is larger than vmax")

    if Vp >= Vths:
        raise ValueError("Vp is larger than Vths")

    cm = mpl.colormaps['RdYlBu_r']
    col_list = [cm(i)[0:3] for i in [0.15, 0.5, 0.8, 0.9]]
    col_list = [cm(0.15)[0:3], mpl.colormaps['jet'](0.67), cm(0.8)[0:3], cm(0.9)[0:3]]
    col_list.insert(0, (0, 0, 0))
    col_list.insert(3, mpl.colors.to_rgb('w'))
    col_list.append(mpl.colormaps['turbo'](0.95))

    # Calculate float_list with values in strictly increasing order
    float_list = [0, mpl.colors.LogNorm(vmin, vmax)(Vp), mpl.colors.LogNorm(vmin, vmax)(V0)]
    float_list += list(np.linspace(mpl.colors.LogNorm(vmin, vmax)(Vths), 1, 4))
    float_list = sorted(set(float_list))  # Ensure uniqueness and sort

    # Ensure float_list is strictly increasing
    if any(float_list[i] >= float_list[i + 1] for i in range(len(float_list) - 1)):
        raise ValueError("float_list values must be strictly increasing")

    cmap_n = get_continuous_cmap(col_list, input_hex=False, float_list=float_list)

    return cmap_n



def get_continuous_cmap(col_list,input_hex=False,float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in col_list.
        If float_list is provided, each color in col_list is mapped to the respective location in float_list.

        Parameters
        ----------
        col_list: list of color code strings
        float_list: list of floats between 0 and 1, same length as col_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    if input_hex:
        rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in col_list]
    else:
        rgb_list = col_list.copy()

    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mpl.colors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]
