#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This creates and a single NetCDF file containing bottom pressure
and other fields, for some time range.

"""

from datetime import datetime, timedelta
start_time = datetime.now()
import netCDF4 as nc

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
import zrfun
import zfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

n_layer = 0 # 0 = deepest
do_press = True

# model_type = 'LiveOcean'
model_type = 'Kurapov'

if model_type == 'LiveOcean':
    gridname = 'cascadia1'
    tag = 'base'
    ex_name = 'lobio5'
    list_type = 'low_passed'
    #
    Ldir = Lfun.Lstart(gridname, tag)
    Ldir['ex_name'] = ex_name
    Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']
    tag = Ldir['gtagex']
    #
    if Ldir['lo_env'] == 'pm_mac':
        dt0 = datetime(2017,8,5)
        dt1 = datetime(2017,8,9)
    elif Ldir['lo_env'] == 'pm_fjord':
        dt0 = datetime(2017,1,1)
        dt1 = datetime(2017,12,31)
    #
    if list_type == 'low_passed':
        fn_list = []
        dt = dt0
        while dt <= dt1:
            date_string = dt.strftime(format='%Y.%m.%d')
            Ldir['date_string'] = date_string
            f_string = 'f' + Ldir['date_string']
            in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
            if 'low_passed.nc' in os.listdir(in_dir):
                fn_list.append(in_dir + 'low_passed.nc')
            else:
                print('Missing file for ' + date_string)
            dt = dt + timedelta(days=1)
    #
    # make some things
    fn = fn_list[0]
    ds = nc.Dataset(fn)
    G, S, T = zrfun.get_basic_info(fn)
    h = ds['h'][:]
    z = zrfun.get_z(h, 0*h, S, only_rho=True)
    z0 = z[n_layer,:,:].squeeze()
    ds.close()

elif model_type == 'Kurapov':
    Ldir = Lfun.Lstart() #checks for _user directory and changes paths accordingly
    in_dir = Ldir['parent'] + 'ptools_data/Kurapov/' #location of .nc files
    tag = 'kurapov'

    # get grid info
    fng = in_dir + 'zts_Monterey_Exp41_3000.nc' #all files contain grid info
    dsg = nc.Dataset(fng) #gets metadata
    if False: #probably some kind of check
        print('\nGRID INFO')
        zfun.ncd(dsg)
    lon = dsg['lon_rho'][:,:] #probably taking a cutout from full model domain
    lat = dsg['lat_rho'][:,:] #[1030-1:1521, 375-1:615]
    dsg.close()

    # make file list
    fn_list = []
    if Ldir['lo_env'] == 'pm_mac': #indicates testing subset?
        #frange = range(3000,3013+1)
        frange = range(3000,3697+1)
    elif Ldir['lo_env'] == 'pm_fjord': #indicates full dataset?
        # we have 3000-3697, so use range(3000,3697+1)
        frange = range(3000,3697+1)
        #frange = range(1,20+1)
    for ii in frange:  #building file names
        nn = ('0000' + str(ii))[-4:]
        fn_list.append(in_dir + 'zts_Monterey_Exp41_' + nn + '.nc')

    # make some things
    S_info_dict = {'VTRANSFORM':2, 'VSTRETCHING':4, 'THETA_S':7,
        'THETA_B':4, 'TCLINE':50, 'N':40} #updated to reflect new files, but what is N?
    S = zrfun.get_S(S_info_dict)
    G = zrfun.get_basic_info(fng, only_G=True) #zrfun had to be modified to reflect variables in .nc files. But what are 'pm' and 'pn'?
    fn = fn_list[0]
    ds = nc.Dataset(fn)
    h = ds['h'][:]
    z = zrfun.get_z(h, 0*h, S, only_rho=True)
    z0 = z[n_layer,:,:].squeeze()
    ds.close()

if Ldir['lo_env'] == 'pm_mac':
    etag = '_mac'
elif Ldir['lo_env'] == 'pm_fjord':
    etag = ''

NT = len(fn_list)

# prepare a directory for results
outdir0 = Ldir['parent'] + 'ptools_output/slow_slip/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + 'bottom_pressure_extractions/'
Lfun.make_dir(outdir, clean=False)
# output file
out_name = 'pressure_' + tag + etag + '.nc'
out_fn = outdir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# make bottom pressure
if do_press:
    g = 9.81
    for tt in range(NT): #for each .nc file
        fn = fn_list[tt]
        ds1 = nc.Dataset(fn)
        if np.mod(tt,10)==0: #print to console every 10th file
            print('tt = ' + str(tt) + '/' + str(NT) + ' ' + str(datetime.now()))
            sys.stdout.flush()
        if model_type == 'LiveOcean':
            zeta = ds1['zeta'][0,:,:].squeeze()
            if tt == 0:
                bp_arr = (0*zeta) * np.ones((NT,1,1))
                DA = G['DX'] * G['DY']
                DAm = np.ma.masked_where(zeta.mask, DA)
            rho = ds1['rho'][0,:,:,:].squeeze() + 1000.
            # note that rho appears to be in situ density, not potential density
            z_w = zrfun.get_z(G['h'], zeta, S, only_w=True)
            DZ = np.diff(z_w, axis=0)
            bp_arr[tt,:,:] = (g * rho * DZ).sum(axis=0)
        elif model_type == 'Kurapov':
            import seawater
            zeta = ds1['zeta'][0,:,:].squeeze()
            if tt == 0:
                bp_arr = (0*zeta) * np.ones((NT,1,1))
                DA = G['DX'] * G['DY'] #are these grid step size? used to calculate area?
                DA = DA[:,:] #matches the (cutout?) range used earlier
                DAm = np.ma.masked_where(zeta.mask, DA)
            salt = ds1['salt'][0,:,:,:].squeeze()
            ptemp = ds1['temp'][0,:,:,:].squeeze() # potential temperature
            z_r = zrfun.get_z(G['h'][:,:], 0* zeta, S, only_rho=True)
            z_w = zrfun.get_z(G['h'][:,:], zeta, S, only_w=True)
            p = seawater.pres(-z_r, G['lat_rho'][0,0])
            temp = seawater.temp(salt, ptemp, p) # in situ temperature
            # for some reason seawater.dens throws errors if we don't do this
            sd = salt.data
            td = temp.data
            prd = p.data
            sd[salt.mask] = np.nan
            td[salt.mask] = np.nan
            prd[salt.mask] = np.nan
            prho = seawater.pden(sd, td, prd) # potential density
            prho = np.ma.masked_where(salt.mask, prho)
            rho = seawater.dens(sd, td, prd) # in-situ density from salinity, temperature, and pressure
            rho = np.ma.masked_where(salt.mask, rho)
            DZ = np.diff(z_w, axis=0)
            bp_arr[tt,:,:] = (g * rho * DZ).sum(axis=0) # sums vertically
        ds1.close()
    bp_mean = np.mean(bp_arr, axis=0)
    bp_anom = bp_arr - bp_mean

# initialize output Dataset
ds2 = nc.Dataset(out_fn, 'w')

# lists of variables to process
if model_type == 'LiveOcean':
    dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
elif model_type == 'Kurapov':
    dlist = ['xi_rho', 'eta_rho', 'ocean_time']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'mask_rho', 'h']
vn_list2t = ['zeta', 'ocean_time']
vn_list3t = ['salt', 'temp']

ds1 = nc.Dataset(fn_list[0])

# Copy dimensions
for dname, the_dim in ds1.dimensions.items():
    if dname in dlist:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# Copy variables
if model_type == 'Kurapov':
    dsg = nc.Dataset(fng)
for vn in vn_list2:
    if model_type == 'LiveOcean':
        varin = ds1[vn]
    elif model_type == 'Kurapov':
        varin = dsg[vn] #copy in coordinates and mask
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    if model_type == 'LiveOcean':
        vv[:] = ds1[vn][:]
    elif model_type == 'Kurapov':
        vv[:] = dsg[vn][:,:] #cutout again?
if model_type == 'Kurapov':
    dsg.close()
#
for vn in vn_list2t:
    varin = ds1[vn] #zeta and time
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.long_name = varin.long_name
    vv.units = varin.units
    try:
        vv.time = varin.time
    except AttributeError:
        # ocean_time has no time
        pass

# copy data
tt = 0
for fn in fn_list: #for each data file
    ds = nc.Dataset(fn)
    ds2['ocean_time'][tt] = ds['ocean_time'][0].squeeze()
    ds2['zeta'][tt,:,:] = ds['zeta'][0, :, :].squeeze()
    tt += 1
    ds.close()
#
for vn in vn_list3t:
    do_var = True
    try:
        varin = ds1[vn] #salinity and temperature
    except IndexError:
        # designed so that it still works when we don't have a variable from this list
        # e.g. when there is no bio or carbon
        do_var = False
        print(' - Variable not found: ' + vn)
    if do_var==True:
        dd = tuple([d for d in varin.dimensions if d != 's_rho'])
        vv = ds2.createVariable(vn, varin.dtype, dd)
        vv.long_name = varin.long_name
        try:
            vv.units = varin.units
        except AttributeError:
            # salt has no units
            pass
        vv.time = varin.time
        for tt in range(NT):
            fn = fn_list[tt]
            ds = nc.Dataset(fn)
            vv[tt,:,:] = ds[vn][0, n_layer, :, :].squeeze()
            ds.close()

# Add derived variables
vv = ds2.createVariable('z', float, ('eta_rho', 'xi_rho'))
vv.long_name = 'z position closest to free surface for 3D variables '
vv.units = 'meter'
vv[:] = z0
#
if do_press: #copy bottom pressure and anomaly
    vv = ds2.createVariable('bp', float, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'bottom pressure'
    vv.units = 'Pa'
    vv.time = 'ocean_time'
    vv[:] = bp_arr
    #
    vv = ds2.createVariable('bpa', float, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'bottom pressure anomaly'
    vv.units = 'Pa'
    vv.time = 'ocean_time'
    vv[:] = bp_anom

ds1.close()
ds2.close()

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

if True:
    zfun.ncd(out_fn)

if Ldir['lo_env'] == 'pm_mac': #make plots if working on subset?
    import matplotlib.pyplot as plt
    plt.close('all')
    fs = 12 # primary fontsize
    lw = 1 # primary linewidth

    # confirm profiles look right
    fig0 = plt.figure(figsize=(11,8.5))
    ax0 = fig0.add_subplot(131)
    ax0.plot(ptemp[:,15,5],z_r[:,15,5],'b--',lw=lw,label='potential')
    ax0.plot(temp[:,15,5],z_r[:,15,5],'b',lw=lw,label='in-situ')
    #ax0.set_yticks(np.arange(0,-2000,-250)) can't get this to work
    #ax0.set_yticklabels(np.arange(0,2000,250))
    ax0.legend(fontsize=fs-2, loc='upper left')
    ax0.set_xlabel('(C)',fontsize=fs)
    ax0.set_ylabel('Depth (m)',fontsize=fs)
    ax0.set_title('Temperature',fontsize=fs+2)
    ax1 = fig0.add_subplot(132)
    ax1.plot(prho[:,15,5],z_r[:,15,5],'b--',lw=lw,label='potential')
    ax1.plot(rho[:,15,5],z_r[:,15,5],'b',lw=lw,label='in-situ')
    ax1.legend(fontsize=fs, loc='upper right')
    ax1.set_xlabel('(g/mL)',fontsize=fs)
    ax1.set_title('Density',fontsize=fs+2)
    ax2 = fig0.add_subplot(133)
    ax2.plot(sd[:,15,5],z_r[:,15,5],'b',lw=lw,label='in-situ')
    ax2.legend(fontsize=fs, loc='upper right')
    ax2.set_xlabel('(psu)',fontsize=fs)
    ax2.set_title('Salinity',fontsize=fs+2)
    fig0.savefig(outdir + 'profiles.png')

    fig = plt.figure(figsize=(8.5,11))
    ax = fig.add_subplot(111)
    ds = nc.Dataset(out_fn)
    vn = 'bp'
    pch = ax.pcolormesh(ds['lon_rho'][:], ds['lat_rho'][:], ds[vn][-1,:,:].squeeze()/10065)
    ax.plot(ds['lon_rho'][15,5], ds['lat_rho'][15,5],'xk',markersize=10,markeredgewidth=2)
    fig.colorbar(pch)
    ax.set_title(ds[vn].long_name + ' (dbar)')
    pfun.dar(ax)
    fig.savefig(outdir + 'mapview.png')

    # put data into a DataFrame
    import pandas as pd
    bpam = ds['bpa'][:,15,5].squeeze()/100.65
    tm = ds['ocean_time'][:]
    dtm_list = []
    for t in tm:
        dtm_list.append(Lfun.modtime_to_datetime(t))
    dti = pd.to_datetime(dtm_list)
    dti = dti.tz_localize('UTC')
    df = pd.DataFrame(data={'bpa':bpam}, index = dti)
    df.index.name = 'Date'
    fig1 = plt.figure(figsize=(11,8.5))
    axa = fig1.add_subplot(111)
    df.plot(ax=axa)
    axa.set_ylabel('Bottom Pressure (cm)')
    fig1.savefig(outdir + 'timeseries.png')

    plt.show()
    ds.close()
