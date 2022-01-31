#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Fri Aug 23 12:52:55 EDT 2019
#readme: update from potential_intensity.py: input p units: hPa -> Pa
if __name__ == '__main__':
    import sys
    from misc.timer import Timer
    tt = Timer('strat ' + ' '.join(sys.argv))
#import os, os.path, sys
#import xarray as xr, numpy as np, pandas as pd
#import matplotlib.pyplot as plt
import numpy as np, xarray as xr

#from xtci.shared.pcmin import pcmin3
from tcpypi.pi import pi as pi1d
if __name__ == '__main__':
    tt.check('end import')

def potential_intensity(sst, slp, p, T, q, dim_z, CKCD=0.9,ascent_flag=0,diss_flag=1,V_reduc=0.8,ptop=50,miss_handle=1):
    '''xarray-wrapper of the tcpypi calculation of PI.
    sst: sea surface temperature in K;
    slp: seal level pressure in Pa;
    p: pressure levels in Pa;
    T: temperature in K;
    q: specific humidity in kg/kg;
    dim_z: dim name along the z/p direction e.g. "pfull" for FLOR/AM2.5.
    kwargs used in tcpypi: CKCD=0.9,ascent_flag=0,diss_flag=1,V_reduc=0.8,ptop=50,miss_handle=1.
    '''
    sstc = sst - 273.15 #units convert for tcpypi: K->C
    slp = slp/100 #units convert for tcpypi: Pa->hPa
    p = p/100 #units convert for tcpypi: Pa->hPa
    Tc = T - 273.15 #units convert for tcpypi: K->C
    r_v = q/(1-q) # specific humidity to mixing ratio
    r_v = r_v * 1000 # units convert from kg/kg to g/kg for tcpypi
    vmax, pmin, iflag, To, p_o = xr.apply_ufunc(pi1d,
        sstc, slp, p, Tc, r_v,
        kwargs=dict(CKCD=CKCD,ascent_flag=ascent_flag,diss_flag=diss_flag,V_reduc=V_reduc,ptop=ptop,miss_handle=miss_handle),
        input_core_dims=[[], [], [dim_z], [dim_z], [dim_z]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True)
    pmin.attrs['long_name'] = 'mininum central pressure'
    pmin.attrs['units'] = 'hPa'
    vmax.attrs['long_name'] = 'maximum surface wind speed'
    vmax.attrs['units'] = 'm/s'
    iflag = iflag.astype('int32')
    iflag.attrs['long_name'] = '1: OK; 0: no convergence; 2: CAPE routine failed; 3: CAPE missing data'
    To.attrs['long_name'] = 'outflow temperature'
    To.attrs['units'] = 'K'
    p_o.attrs['long_name'] = 'outflow temperature level'
    p_o.attrs['units'] = 'hPa'
    PI = xr.Dataset(dict(pmin=pmin, vmax=vmax, iflag=iflag, To=To, p_o=p_o))
    
    return PI

if __name__ == '__main__':
    ifile = '/tigress/wenchang/MODEL_OUT/CTL1860_noleap_tigercpu_intelmpi_18_576PE/POSTP/10000101.atmos_month.nc'
    ds = xr.open_dataset(ifile)#.sel(pfull=slice(100, None))
    is_ocean = ds.land_mask.load() < 0.1
    sst = ds.t_surf.where(is_ocean)
    slp = ds.slp.where(is_ocean)*100
    i_reversed = slice(-1, None, -1)
    p = ds.pfull.isel(pfull=i_reversed)*100
    T = ds.temp.where(is_ocean).isel(pfull=i_reversed).load()
    q = ds.sphum.where(is_ocean).isel(pfull=i_reversed).load()
    print('calculating...')
    PI = potential_intensity(sst, slp, p, T, q, dim_z='pfull')
    print('saving...')
    PI.to_netcdf('PI_tcpypi.nc', encoding={dname:{'zlib': True, 'complevel': 1} for dname in list(PI.data_vars)})

    tt.check('**done**')
