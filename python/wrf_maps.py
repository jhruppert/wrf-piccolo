# %% [markdown]
# #### Jupyter notebook to create snapshot maps from WRF output.
# 
# James Ruppert  
# jruppert@ou.edu  
# 11/16/23

# %% [markdown]
# #### Main Settings

# %%
from netCDF4 import Dataset
import numpy as np
from matplotlib import ticker, colors, rc#, cm
import matplotlib.pyplot as plt
from wrf import getvar#, get_cartopy
# import cartopy.crs as crs
# import cartopy.feature as cfeature
# from metpy.plots import ctables
from read_wrf_piccolo import *
import xarray as xr
# import matplotlib.animation as animation
import os

# %%
#### Directories and model output specs
scdir = "/glade/derecho/scratch/ruppert/piccolo/"
figdir = "/glade/work/ruppert/wrf-piccolo/python/figures/"

case = "sept1-4/"
memb_tag = "memb_01"
test_name = "ctl/"
wrf_domain="wrf_fine/"

wrfdir = scdir+case+memb_tag+'/'+test_name
wrffiles = get_wrf_file_list(wrfdir, wrf_domain+"wrfout")
lat, lon, nx1, nx2, nz, npd = wrf_dims(wrffiles[0])
nfiles = len(wrffiles)

# %% [markdown]
# ---
# ### Plotting routines

# %%
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 15}

rc('font', **font)

# %% [markdown]
# ##### Plot functions

# %%
# Variable settings

def figure_var_specs(plot_name):

    # Switches (default settings)
    i2d=True      # switch on if the data read-in needs to be done in 2D
    dosym=True    # switch off to specify min colorbar setting
    dolog=False   # switch for logarithmic color scale
    extend='both' # extend color fill beyond bounds
    scale=1.      # scale variable by x
    vartag=plot_name

    if plot_name == 'OLR':
        # OLR
        unittag='W m$^{-2}$'
        cmin=75
        cmax=320
        dosym=False
        cmap='RdGy'
    elif plot_name == 'U10':
        # OLR
        unittag='m/s'
        cmin=-10
        cmax=10
        dosym=True
        cmap='RdBu_r'
    elif plot_name == 'MRef':
        # Composite reflectivity (column-max value)
        vartag='mdbz'
        unittag='dBZ'
        cmin=-25
        cmax=60
        dosym=False
        extend='neither'
        # cmap = ctables.registry.get_colortable('NWSReflectivity')#'NWSReflectivityExpanded'
        cmap='Spectral_r'
    elif plot_name == 'refl_10cm':
        # Composite reflectivity (column-max value)
        vartag='REFL_10CM'
        unittag='dBZ'
        cmin=-25
        cmax=60
        dosym=False
        extend='neither'
        # cmap = ctables.registry.get_colortable('NWSReflectivity')#'NWSReflectivityExpanded'
        cmap='Spectral_r'

    return vartag, unittag, cmin, cmax, extend, cmap, i2d, dosym, dolog, scale

# %%
# wind barbs
def plot_wind(ax, u, v, lon, lat, skip, length=None):#, transform
    spacing=skip #barbspacing (smaller if zoomed in)
    mps_to_kts=1.94384 # conversion factor from m/s to knots for barbs
    uplt = u * mps_to_kts
    vplt = v * mps_to_kts
    # ax.barbs(lon[::spacing,::spacing], lat[::spacing,::spacing], 
    #          uplt[::spacing,::spacing], vplt[::spacing,::spacing], 
    ax.barbs(lon[::spacing], lat[::spacing], 
             uplt[::spacing,::spacing], vplt[::spacing,::spacing], 
             zorder=2, color='black', length=length,
             linewidth=0.8)
            #  transform=crs.PlateCarree(), linewidth=1.5)

# %%
def read_pltvar(wrffile, it_read, plot_name, vartag):
    wrffil_read = Dataset(wrffile.strip())
    # time = getvar(wrffil_read, 'Times')
    time = wrffil_read.variables['Times']
    time = [''.join(row.astype(str)) for row in time]
    time = time[it_read]
    u10 = getvar(wrffil_read, "U10", timeidx=it_read) # m/s
    v10 = getvar(wrffil_read, "V10", timeidx=it_read) # m/s
    if plot_name == "refl_10cm":
        pltvar = wrffil_read.variables[vartag][it_read,0,:,:]
        pltvar = xr.DataArray(pltvar, coords=u10.coords, attrs=u10.attrs)
    else:
        pltvar = getvar(wrffil_read, vartag, timeidx=it_read)
    wrffil_read.close()
    return pltvar, u10, v10, time

# %%
    # Add map features
    # states = cfeature.NaturalEarthFeature(category="cultural", scale="10m",
    #                          facecolor="none",
    #                          name="admin_1_states_provinces_lines")
    # countries = cfeature.NaturalEarthFeature(category="cultural", scale="10m",
    #                          facecolor="none",
    #                          name="admin_0_countries_lakes")
    # featurewidth=0.5
    # featurecol="black"
    # ax.add_feature(states, linewidth=featurewidth, edgecolor=featurecol, zorder=1)
    # ax.add_feature(countries, linewidth=featurewidth, edgecolor=featurecol, zorder=1)

# %%
def run_plot(plot_name, wrffile, it_read, plt_area=[lon[0], lon[-1], lat[0], lat[-1]]):

    vartag, unittag, cmin, cmax, extend, cmap, i2d, dosym, dolog, scale = figure_var_specs(plot_name)
    pltvar, u10, v10, time = read_pltvar(wrffile, it_read, plot_name, vartag)

    hr_tag = str(time)[0:10]+', '+str(time)[11:16]+' UTC'
    title_extra=''
    # if i2d:
    title = plot_name+title_extra+', '+hr_tag
    # else:
    #     title = plot_name+title_extra+', '+hr_tag+',  k-level='+str(ikread)+' (p = '+str(int(pres[ikread]))+' hPa)'

    # Color scale
    nlevs=31#71
    if dosym:
        delta=2*cmax/nlevs
        clevs = np.arange(-1*cmax,cmax+delta,delta)
    else:
        delta=(cmax-cmin)/nlevs
        clevs = np.arange(cmin,cmax+delta,delta)

    # create figure
    fig = plt.figure(figsize=(12,10))
    # proj = get_cartopy(pltvar)
    ax = fig.add_subplot(111)#,projection=proj)

    ax.set_title(title, fontsize=20)
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')

    # fill contour
    if dolog:
        im = ax.contourf(lon, lat, pltvar*scale, cmap=cmap, alpha=0.9,
                            extend='both', zorder=0, norm=colors.LogNorm(vmin=cmin, vmax=cmax))
                            # transform=crs.PlateCarree())
        ticks=ticker.LogLocator()
    else:
        im = ax.contourf(lon, lat, pltvar*scale, clevs, cmap=cmap, alpha=0.9,
                            extend=extend, zorder=0)#, transform=crs.PlateCarree())
        ticks=ticker.AutoLocator()
    cbar = plt.colorbar(im, ax=ax, shrink=0.45, ticks=ticks)
    cbar.ax.set_ylabel(unittag)

    # Set plot area
    # plt.set_extent(plt_area)
    ax.set_xlim(plt_area[0], plt_area[1])
    ax.set_ylim(plt_area[2], plt_area[3])
    ax.set_aspect('equal')

    # Add wind barbs
    # if plot_name == "2mTemp":
    skip = int(u10.shape[0]/30)
    plot_wind(ax, u10, v10, lon, lat, skip, length=5)#, crs.PlateCarree()

    plt.tight_layout()
    # format="pdf"
    format="png"
    hr_tag_fig = str(time)[0:16]
    plt.savefig(figdir+plot_name+"_"+memb_tag+"_"+hr_tag_fig+"."+format, format=format, bbox_inches="tight")
    # plt.show()
    plt.close('all')

# %% [markdown]
# #### Create plots

# %%
# plt_area=[lon[0,0], lon[0,-1], lat[0,0], lat[-1,0]] # W,E,S,N

# for ifile in range(0,nfiles):
for ifile in range(0,1):
    # for it_read in range(0,npd,9):
    for it_read in range(0,npd):
        # for var_name in ["OLR", "MRef", "900-600Thick"]:
        for var_name in ["refl_10cm"]:
            run_plot(var_name, wrffiles[ifile], it_read)#, plt_area=plt_area)

# %%
# Create animation
# var_name="refl_10cm"
# anim_file = figdir+var_name+"_"+memb_tag+".gif"
# convert_str = "convert -delay 50 -loop 0 "+figdir+"/"+var_name+"*.png "+anim_file
# os.system(convert_str)


