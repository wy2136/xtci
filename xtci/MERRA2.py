#!/usr/bin/env python
# Wenchang Yang (wenchang@princeton.edu)
# Sat Aug 24 15:02:58 EDT 2019
#readme: update from old version: potential_intensity -> potential_intensity_tcpypi; entropy_deficit also updated with new arg forGPI2010
if __name__ == '__main__':
    from misc.timer import Timer
    tt = Timer(f'start {__file__}')
import os.path, os, sys
import xarray as xr, numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import absolute, exp, log

from xtci.shared.entropy_deficit import entropy_deficit, relative_humidity
#from xtci.shared.potential_intensity import potential_intensity
from xtci.shared.potential_intensity_tcpypi import potential_intensity
from xtci.shared.wind_shear import wind_shear
from xtci.shared.absolute_vorticity import absolute_vorticity

#if True:
#    year = 1980
#    odir = None
def do_tci(year, odir=None):
    '''calculate TC indices (e.g. GPI, VI) and related variables given MERRA2 monthly reanalysis. Modified from the ERA5 script.'''
    print('[year]:', year)
    if odir is None:
        odir = '.'
    ibasename = f'merra2.monthly.{year}.nc' 

    # sst and ocean mask
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/sea_surface_temperature/era5.sea_surface_temperature.monthly.{year:04d}.nc'
    #sst = xr.open_dataarray(ifile)# units K
    #is_ocean = sst.isel(time=0).drop('time').pipe(lambda x: x*0==0)
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/tavgM_2d_ocn_Nx/TSKINWTR/MERRA2_*.tavgM_2d_ocn_Nx.{year:04d}??.nc4.nc4'
    sst = xr.open_mfdataset(ifile)['TSKINWTR'].load().resample(time='MS').mean()# units K
    is_ocean = sst.pipe(lambda x: x*0==0)
    ##test
    #is_ocean.isel(time=0).plot()
    #plt.show()
    #sys.exit()


    # slp
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/mean_sea_level_pressure/era5.mean_sea_level_pressure.monthly.{year:04d}.nc'
    #slp = xr.open_dataarray(ifile) # units Pa
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_2d_asm_Nx/SLP/MERRA2_*.instM_2d_asm_Nx.{year:04d}??.nc4.nc4'
    slp = xr.open_mfdataset(ifile)['SLP'].load().resample(time='MS').mean() # units Pa
    
    # t2m
    #ifile = f'/tigress/wenchang/data/era5/analysis/2m_temperature/monthly/era5.2m_temperature.monthly.{year:04d}.nc'
    #t2m = xr.open_dataarray(ifile) # units K
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_2d_asm_Nx/T2M/MERRA2_*.instM_2d_asm_Nx.{year:04d}??.nc4.nc4'
    t2m = xr.open_mfdataset(ifile)['T2M'].load().resample(time='MS').mean() # units K 

    # RH2m
    #ifile = f'/tigress/wenchang/data/era5/analysis/2m_relative_humidity/monthly/era5.2m_relative_humidity.monthly.{year:04d}.nc'
    #RH2m = xr.open_dataarray(ifile) # units %
    #QV2M
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_2d_asm_Nx/QV2M/MERRA2_*.instM_2d_asm_Nx.{year:04d}??.nc4.nc4'
    qv2m = xr.open_mfdataset(ifile)['QV2M'].load().resample(time='MS').mean() # units kg/kg
    #PS
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_2d_asm_Nx/PS/MERRA2_*.instM_2d_asm_Nx.{year:04d}??.nc4.nc4'
    ps = xr.open_mfdataset(ifile)['PS'].load().resample(time='MS').mean() # units Pa
    RH2m = relative_humidity(q=qv2m, p=ps, T=t2m).clip(min=0, max=1) * 100 # units %
    
    # Ta
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/temperature/era5.temperature.monthly.{year:04d}.nc'
    #Ta = xr.open_dataarray(ifile) # in K
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_3d_asm_Np/T/MERRA2_*.instM_3d_asm_Np.{year:04d}??.nc4.nc4'
    Ta = xr.open_mfdataset(ifile)['T'].load().resample(time='MS').mean() # in K

    # q
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/specific_humidity/era5.specific_humidity.monthly.{year:04d}.nc'
    #q = xr.open_dataarray(ifile) # in kg/kg
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_3d_asm_Np/QV/MERRA2_*.instM_3d_asm_Np.{year:04d}??.nc4.nc4'
    q = xr.open_mfdataset(ifile)['QV'].load().resample(time='MS').mean() # in kg/kg

    # RH
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/relative_humidity/era5.relative_humidity.monthly.{year:04d}.nc'
    #RH = xr.open_dataarray(ifile) # in %
    p = q.lev*100 #hPa -> Pa
    RH = relative_humidity(q=q, p=p, T=Ta).clip(min=0, max=1) * 100 # units %

    # u
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/u_component_of_wind/era5.u_component_of_wind.monthly.{year:04d}.nc'
    #u = xr.open_dataarray(ifile) # in m/s
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_3d_asm_Np/U/MERRA2_*.instM_3d_asm_Np.{year:04d}??.nc4.nc4'
    u = xr.open_mfdataset(ifile)['U'].load().resample(time='MS').mean() # in m/s 

    # v 
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/v_component_of_wind/era5.v_component_of_wind.monthly.{year:04d}.nc'
    #v = xr.open_dataarray(ifile) # in m/s
    ifile = f'/tigress/wenchang/data/merra2/raw/monthly/instM_3d_asm_Np/V/MERRA2_*.instM_3d_asm_Np.{year:04d}??.nc4.nc4'
    v = xr.open_mfdataset(ifile)['V'].load().resample(time='MS').mean() # in m/s 

    # vorticity
    #ifile = f'/tigress/wenchang/data/era5/raw/monthly/plevels/vorticity/era5.vorticity.monthly.{year:04d}.nc'
    #vort = xr.open_dataarray(ifile) # in s**-1
    ifile = f'/tigress/wenchang/data/merra2/analysis/instM_3d_asm_Np/VORT/MERRA2.vort.monthly.{year:04d}.nc'
    vort = xr.open_dataarray(ifile).load().resample(time='MS').mean() # units s**-1
    ##test
    #vort.isel(time=0, lev=0).plot()
    #plt.show()
    #sys.exit()


    # entropy deficit: (s_m_star - s_m)/(s_sst_star - s_b)
    print('entropy deficit')
    dname = 'chi'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        chi = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        p_m = 600*100 # Pa
        chi = entropy_deficit(
            sst=sst,
            slp=slp,
            Tb=t2m,
            RHb=RH2m/100,
            p_m=p_m,
            Tm=Ta.sel(lev=p_m/100).drop('lev'),
            RHm=RH.sel(lev=p_m/100).drop('lev')/100
            ).where(is_ocean)
        chi.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # entropy deficit for GPI2010: (s_b - s_m)/(s_sst_star - s_b)
    print('entropy deficit for GPI2010')
    dname = 'chi_sb'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        chi = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        p_m = 600*100 # Pa
        chi_sb = entropy_deficit(
            sst=sst,
            slp=slp,
            Tb=t2m,
            RHb=RH2m/100,
            p_m=p_m,
            Tm=Ta.sel(lev=p_m/100).drop('lev'),
            RHm=RH.sel(lev=p_m/100).drop('lev')/100
            ).where(is_ocean)
        chi_sb.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # potential intensity
    print('potential intensity')
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.PI.nc') )
    if os.path.exists(ofile):
        PI = xr.open_dataset(ofile)
        print('[opened]:', ofile)
    else:
        """
        reverse_plevels = lambda x: x.isel(level=slice(-1, None, -1)) # no need to reverse plevels for MERRA2 (lev = 1000, 975, ..., 0.1)
        PI = potential_intensity(
            sst=sst,
            slp=slp.where(is_ocean),
            p=Ta.lev.pipe(reverse_plevels),
            T=Ta.pipe(reverse_plevels).where(is_ocean),
            q=q.pipe(reverse_plevels).where(is_ocean),
            dim_x='longitude', dim_y='latitude', dim_z='level'
            )
        """
        """
        #old version
        PI = potential_intensity(
            sst=sst,
            slp=slp.where(is_ocean),
            p=Ta.lev,
            T=Ta.where(is_ocean),
            q=q.where(is_ocean),
            dim_x='lon', dim_y='lat', dim_z='lev'
            )
        """
        PI = potential_intensity(
            sst=sst,
            slp=slp.where(is_ocean),
            p=Ta.lev*100,
            T=Ta.where(is_ocean),
            q=q.where(is_ocean),
            dim_z='lev'
            )
        encoding = {dname:{'dtype': 'float32', 'zlib': True, 'complevel': 1} 
            for dname in ('pmin', 'vmax')}
        encoding['iflag'] = {'dtype': 'int32'}
        PI.to_netcdf(ofile, encoding=encoding, unlimited_dims='time')
        print('[saved]:', ofile)

    # wind shear: ( (u200-u850)**2 + (v200-v850)**2 )**0.5
    print('wind shear')
    dname = 'Vshear'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        Vshear = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        Vshear = wind_shear(
            u850=u.sel(lev=850),
            v850=v.sel(lev=850),
            u200=u.sel(lev=200),
            v200=v.sel(lev=200)
            )
        Vshear.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # ventilation index: Vshear * chi_m /V_PI
    print('ventilation index')
    dname = 'VI'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        VI = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        VI = Vshear*chi/PI.vmax.pipe(lambda x: x.where(x>0))
        VI.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # absolute vorticity at 850hPa
    print('absolute vorticity')
    dname = 'eta'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        eta = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        eta = absolute_vorticity(
            vort850=vort.sel(lev=850),
            lat=vort.lat
            )
        eta.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)


    # relative humidity at 600hPa in %
    print('relative humidity in %')
    dname = 'H'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        H = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        H = RH.sel(lev=600).drop('lev')
        H.attrs['long_name'] = '600hPa relative humidity'
        H.attrs['units'] = '%'
        H.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)
    
    # GPI (Emanuel and Nolan 2004): |10**5\eta|**(3/2) * (H/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)
    print('GPI')
    dname = 'GPI'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        GPI = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        GPI = (1e5 * absolute(eta) )**(3/2) \
            * (H/50)**3 \
            * (PI.vmax/70)**3 \
            * (1+0.1*Vshear)**(-2)
        GPI.attrs['long_name'] = 'Genesis Potential Index'
        GPI.attrs['history'] = '|10**5\eta|**(3/2) * (H/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)'
        GPI.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # GPI2010 (Emanuel 2010): |\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)
    print('GPI2010')
    dname = 'GPI2010'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        GPI2010 = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        GPI2010 = absolute(eta)**3 \
            * chi_sb.where(chi_sb>0)**(-4/3) \
            * (PI.vmax - 35).clip(min=0)**2 \
            * (25 + Vshear)**(-4)
        GPI2010.attrs['long_name'] = 'Genesis Potential Index of Emanuel2010'
        GPI2010.attrs['history'] = '|\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)'
        GPI2010.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)
    
if __name__ == '__main__':
    #year = 1979
    #do_tci(year, odir='/tigress/wenchang/data/era5/analysis/TCI/')
    odir = '/tigress/wenchang/data/merra2/analysis/TCI'
    if len(sys.argv)>1:
        year = int(sys.argv[1])
        do_tci(year, odir=odir)
    else:
    years = range(1980, 2021)
        for year in years:
            #do_tci(year, odir='/tigress/wenchang/data/era5/analysis/TCI/')
            do_tci(year, odir=odir)

    tt.check('**done**')
