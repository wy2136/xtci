#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Sat Aug 24 15:02:58 EDT 2019
#import os.path#, os, sys
#import xarray as xr, numpy as np#, pandas as pd
#import matplotlib.pyplot as plt
from xtci.ERA5 import do_tci

    
if __name__ == '__main__':
    #year = 1979
    #do_tci(year, odir='/tigress/wenchang/data/era5/analysis/TCI/')
    years = range(1979, 2019)
    odir='/tigress/wenchang/data/era5/analysis/TCI/'
    for year in years:
        do_tci(year, odir=odir)
