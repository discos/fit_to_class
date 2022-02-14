#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:22:33 2020

@author: debbio
"""


import argparse
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
from os import path
import sys


if __name__ == "__main__" :
    l_parser= argparse.ArgumentParser( \
                    description= "classfit offset CDELT2 CDELT3 plot to arcsec"\
                    )
    l_parser.add_argument("-f", "--file", type= str, help= "Input file")    
    l_parser.add_argument("-u", "--unit", type= str, help= "Input data unit", \
                          choices= ['deg', 'rad'], default= 'deg')
    args= l_parser.parse_args()
    # Check input file
    if not path.exists(args.file):
        print("Input file %s not found", args.file)
        sys.exit(-1)
    # Fits open
    try:
        l_fits= fits.open(args.file)
    except Exception as e:
        print(e)
        sys.exit(-1)
    # convert data
    l_cdelt_2= l_fits[1].data['cdelt2']    
    l_cdelt_3= l_fits[1].data['cdelt3']
    if args.unit == 'deg':
        l_angles_2= Angle(l_cdelt_2* u.deg).arcsec
        l_angles_3= Angle(l_cdelt_3* u.deg).arcsec
    elif args.unit == 'rad':
        l_angles_2= Angle(l_cdelt_2* u.rad).arcsec
        l_angles_3= Angle(l_cdelt_3* u.rad).arcsec
    # Plot
    fig, axs= plt.subplots(2)
    fig.suptitle('Ra Dec offsets from object coordinates (arcsec)')
    axs[0].plot(l_angles_2)
    axs[1].plot(l_angles_3)
    plt.show()
    
    