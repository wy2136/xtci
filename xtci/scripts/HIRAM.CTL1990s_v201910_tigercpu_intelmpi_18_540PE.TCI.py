#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Mon Feb 17 11:50:04 EST 2020
import sys, os.path, os, datetime
#import xarray as xr, numpy as np, pandas as pd
#import matplotlib.pyplot as plt
from xtci.AMx import do_tci
print()
model = 'HIRAM'
expname = 'CTL1990s_v201910_tigercpu_intelmpi_18_540PE'
years = range(101,150+1)

exp_user = 'wenchang'
analysis_user = os.environ['USER']

 
if __name__ == '__main__':
    tformat = '%Y-%m-%d %H:%M:%S'
    t0 = datetime.datetime.now()
    print('[start]:', t0.strftime(tformat))
    
    idir = f'/tigress/{exp_user}/MODEL_OUT/{model}/{expname}'
    odir = f'/tigress/{analysis_user}/MODEL_OUT/{model}/{expname}/analysis_wy/TCI'
    if not os.path.exists(odir):
        os.makedirs(odir)
        print('[dir made]:', odir)
    for year in years:
        print(f'year = {year:04d}')
        ifile = os.path.join(idir, 'POSTP', f'{year:04d}0101.atmos_month.nc')
        do_tci(ifile, odir)

    t1 = datetime.datetime.now()
    print('[total time used]:', f'{(t1-t0).seconds:,} seconds')
    print('[end]:', t1.strftime(tformat))
    print()

