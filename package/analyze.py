import warnings
warnings.filterwarnings('ignore')
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from marvin.tools.maps import Maps
from marvin.tools import Cube
from marvin.tools import Image
import copy
import marvin.utils.plot.map as mapplot
import os
from os.path import join
from marvin.tools.query import Query
import pandas as pd




# ------------------------------------------------------------------------------------------------
def plot_image(object_name):
    image = Image(object_name)
    # overlay the IFU fibers
    ax = image.plot()

    image.overlay_fibers(ax)
	 
def load_maps(object_name):
    maps = Maps(object_name)
    # get stellar velocity and sigma map
    s_vel = maps['stellar_vel']
    s_sig = maps['stellar_sigma']
    g_vel = maps['emline_gvel_ha_6564']
    g_sig = maps['emline_gsigma_ha_6564']
    
    # get differences
    v_diff = s_vel - g_vel
    s_diff = s_sig - g_sig
    
    maps = [s_vel, g_vel, s_sig, g_sig, v_diff, s_diff]
    
    return maps
	 
def plot_maps(maps):
    fig, axes = plt.subplots(3, 2, figsize=[7,9])
    # cosmetics
    titles = ['Stellar Velocity','H$\\alpha$ Velocity',
              'Stellar Velocity\nDispersion','H$\\alpha$ Velocity\nDispersion',
              'Stars - Gas Vel.', 'Stars - Gas Sig.']
    ### color bars
    axes = axes.ravel()
    for i, (ax, map_) in enumerate(zip(axes, maps)):
        if i%2==0:
            lo = np.percentile(np.hstack([maps[i].value,maps[i+1].value])[:],5)
            hi = np.percentile(np.hstack([maps[i].value,maps[i+1].value])[:],95)
        if i > 3:
            mapplot.plot(dapmap=map_, fig=fig, ax=ax, cmap = 'RdBu',cbrange=(lo,hi),cb_kws={'shrink':0.4},symmetric=True)
        else:
            mapplot.plot(dapmap=map_, fig=fig, ax=ax,cbrange=(lo,hi),cb_kws={'shrink':0.4})

    for ax, title in zip(axes, titles):
        ax.set_title(title)
    plt.tight_layout()
    
    
def stack_spectra(cube, ind):
    # get a list of spaxels by passing a list of x-indices and a list of y-indices
    spaxels = cube.getSpaxel(x=ind[0], y=ind[1], xyorig='lower')

    # copy the first spaxel so that the stack is a Marvin Spaxel object
    stack = copy.deepcopy(spaxels[0])
    wl = stack.flux.wavelength

    # overwrite the flux with the mean flux of the spectra in the stack (but keep the other meta-data)
    stack.flux = np.mean([sp.flux for sp in spaxels], axis=0)

    return wl, stack.flux

def get_stacked_spectra(object_name):
    cube = Cube(object_name)
    
    # get positive and negative indices
    _,_,_,_,v_diff,_ = load_maps(object_name)
    ind_neg = np.where((v_diff.value > -70) & (v_diff.value < -40))
    ind_pos = np.where((v_diff.value > 40) & (v_diff.value < 70))
    
    # get wavelength and stack
    wavelength, pos_stack = stack_spectra(cube, ind_pos)
    wavelength, neg_stack = stack_spectra(cube, ind_neg)
    
    return (wavelength, pos_stack, neg_stack)
    
def get_redshift(object_name):
    drpall = fits.open('package/drpall-v2_1_2.fits')
    tbdata = drpall[1].data
    ind = np.where(tbdata['plateifu'] == object_name)
    drpall.close()
    if len(tbdata['nsa_z'][ind]) == 0:
        print('redshift not found :(')
        print('guessing...(beware of incorrect plot!)')
        return 0.05, True
    else:
        return tbdata['nsa_z'][ind][0], False
	 
    


def plot_features(object_name):
    # get the stacked spectra
    w, p, n = get_stacked_spectra(object_name)
    
    # rescale to be about the same amplitude
    n *= np.nanmean(p/n)
    
    # get redshift
    z, flag = get_redshift(object_name)
    
    # correct redshift
    w *= 1/(1+z)
    
    # plot
    fig, axes = plt.subplots(nrows=2,ncols=1,dpi=150,figsize=[8,6])
    cmap = plt.get_cmap('bwr')
    # h alpha
    axes[0].plot(w,p,color=cmap(0.9),label='+V')
    axes[0].plot(w,n,color=cmap(0.1),label='-V')
    sub_lo = np.min(p[np.argwhere((w.value>6473) & (w.value<6673))])-0.075
    sub_hi = np.max(p[np.argwhere((w.value>6473) & (w.value<6673))])+0.075
    axes[0].set_ylim(sub_lo,sub_hi)
    axes[0].set_xlim(6473,6673)
    axes[0].legend(loc='upper left',fontsize=15)
    axes[0].set_title('H$\\alpha$ Emission Lines',fontsize=20)
    axes[0].tick_params(labelsize=14)
    axes[0].axvspan(6563-0.25,6563+0.25,color='gray')

    # calcium h and k
    axes[1].plot(w,p,color=cmap(0.9),label='+V')
    axes[1].plot(w,n,color=cmap(0.1),label='-V')
	 # get min between x lims
    sub_lo = np.min(p[np.argwhere((w.value>3850) & (w.value<4050))])-0.075
    sub_hi = np.max(p[np.argwhere((w.value>3850) & (w.value<4050))])+0.075
    axes[1].set_ylim(sub_lo,sub_hi)
    axes[1].set_xlim(3850,4050)
    axes[1].legend(loc='upper left',fontsize=15)
    axes[1].set_title('Calcium H + K \nAbsorption Lines',fontsize=20)
    axes[1].tick_params(labelsize=14)
    axes[1].axvspan(3969-0.25,3969+0.25,color='gray')
    axes[1].axvspan(3934-0.25,3934+0.25,color='gray')


    for ax in axes.flat:
        ax.set_xlabel('Wavelength [angstroms]',fontsize=16)
        ax.set_ylabel('Flux [arb. scaling]',fontsize=16)

    # global
    plt.tight_layout()
	 
def get_query():
	myfilter = 'nsa.z < 0.1 and nsa.z > 0.05'
	query = Query(search_filter=myfilter, mode='remote', limit=1000)
	results = query.run()
	random_number = np.random.randint(low=0,high=1000)
	string = results.results[random_number].plateifu
	
	return string
	
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

def load_maps2(object_name):
    maps = Maps(object_name)
    return maps
	 
def plot_bpt(maps,save_dir):
	 a = maps.get_bpt(return_figure=True,show_plot=True)
	 a[1].savefig(save_dir + '/' + maps.plateifu + '_bpt.pdf')

def plot_maps_side_by_side(maps):
    OH = get_metallicity(maps, star_forming=False)
    OH_masked, mask = get_metallicity(maps, star_forming=True)
    mapps = [OH, OH_masked]
    masks = [None, mask]
    titles=['All Spaxels','Star Forming\nSpaxels']

    fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
    for ax, mask, title, m in zip(axes, masks, titles, mapps):
        m.plot(fig=fig, ax=ax, mask=mask,cblabel='12+log(O/H)',title=title)
    
    fig.tight_layout()
	 
def get_spaxel_size(maps):
    # get size in arcsec
    spaxel_size = float(maps.getCube().header['CD2_2']) * 3600
    # get distance
    redshift = get_redshift(maps.plateifu)
    c = 299792  # speed of light [km/s]
    H0 = 70  # [km s^-1 Mpc^-1]
    D = c * redshift[0] / H0  # approx. distance to galaxy [Mpc]
    scale = 1 / 206265 * D * 1e6  # 1 radian = 206265 arcsec [pc / arcsec]
    spaxel_area = (scale * spaxel_size)**2  # [pc^2]
    
    return spaxel_area

def get_surface_density(maps):
    write_data_file_from_firefly(maps.plateifu)
    mstar = pd.read_csv('data/manga-{}_mstar.csv'.format(maps.plateifu))
    spaxel_area = get_spaxel_size(maps)
    
    sigma_star = np.log10(10**mstar / spaxel_area)  # [Msun / pc^2]
    
    return sigma_star

	 
def write_data_file_from_firefly(object_name):
    
    '''
    https://github.com/sdss/marvin/issues/689 
    '''
    fin = fits.open('data/manga_firefly-v2_1_2-STELLARPOP.fits')
    plateifus = fin['GALAXY_INFO'].data['PLATEIFU']

    # Ngal x Nybin x Nxbin x
    #     (binid,
    #      xmin [arcsec],
    #      xmax [arcsec],
    #      ymin [arcsec],
    #      ymax [arcsec],
    #      image size [units of spaxel number])
    spaxel_binid = fin['SPAXEL_BINID'].data

    # Ngal x Nbin x (Mstar, Mstar_err)
    mstar_all = fin['STELLAR_MASS_VORONOI'].data

    # Select galaxy and binids
    ind1 = np.where(plateifus == object_name)[0][0]
    ind_binid = spaxel_binid[ind1, :, :, 0].astype(int)

    # Create 2D stellar mass array
    mstar = np.ones(ind_binid.shape) * np.nan
    for row, inds in enumerate(ind_binid):
        ind_nans = np.where(inds == -99)
        mstar[row] = mstar_all[ind1, inds, 0]
        mstar[row][ind_nans] = np.nan

    # trim mstar to match size of DAP maps and write to csv
    cube = Cube(object_name)
    len_x = int(cube.header['NAXIS1'])

    df = pd.DataFrame(mstar[:len_x, :len_x])
    fout = join('data/', 'manga-{}_mstar.csv'.format(object_name))
    df.to_csv(fout, index=False)


    fin.close()
	 
def get_metallicity(maps, star_forming=True):
    # get nii and ha 
    nii = maps['emline_gflux_nii_6585']
    ha  = maps['emline_gflux_ha_6564']
    ratio = nii/ha
    log_ratio = np.log10(ratio)
    # use Pettini + Pagel (2004) calibration
    OH = 8.90+0.57*log_ratio
    
    
    if star_forming:
        # get the star forming, good data masks
        masks_bpt, __, __ = maps.get_bpt(show_plot=False)
        mask_non_sf = ~masks_bpt['sf']['global'] * ratio.pixmask.labels_to_value('DONOTUSE')
        mask_bad_data = ratio.pixmask.get_mask(['NOCOV', 'UNRELIABLE', 'DONOTUSE'])
        # signal to noise threshold
        min_snr = 3.
        mask_nii_low_snr = (np.abs(nii.value * np.sqrt(nii.ivar)) < min_snr)
        mask_ha_low_snr = (np.abs(ha.value * np.sqrt(ha.ivar)) < min_snr)
        mask = mask_non_sf | mask_bad_data | mask_nii_low_snr | mask_ha_low_snr
        return OH, mask
        
    else:
        return OH
		  
def plot_mzr(maps,save_dir):
    # get data
    sigma_star = get_surface_density(maps)
    OH,mask = get_metallicity(maps,True)
    mstar = pd.read_csv('data/manga-{}_mstar.csv'.format(maps.plateifu))
    
    # fitting formula Barrera-Ballesteros et al. (2016)
    aa = 8.55
    bb = 0.014
    cc = 3.14
    xx = np.linspace(1, 3, 1000)
    yy = aa + bb * (xx - cc) * np.exp(-(xx - cc))
    
    # plot
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize=[9,6],dpi=125)
    plt.subplot(121)
    im = plt.imshow(mstar, origin='lower')
    plt.xlabel('spaxel')
    plt.ylabel('spaxel')
    plt.colorbar(im,label='log(Stellar Mass) [M$_\odot$]',fraction=0.046, pad=0.04)
    
    plt.subplot(122)
    im = plt.imshow(sigma_star, origin='lower',cmap='Spectral_r')
    plt.xlabel('spaxel')
    plt.ylabel('spaxel')
    plt.colorbar(im,label='log(Stellar Surface Density) [M$_\odot$ / pc$^{2}$]',fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    plt.savefig(save_dir + '/' + maps.plateifu + '_mzr.pdf')
	 
    plt.figure(figsize=(4,4),dpi=150)
    plt.scatter(sigma_star.values[mask == 0], OH.value[mask == 0], alpha=0.15,color='k',label='Spaxels')
    plt.plot(xx, yy,label='Barrera-Ballesteros et al. (2016) fit')
    plt.xlabel('log(Stellar Surface Density) [M$_\odot$ / pc$^{2}$]')
    plt.ylabel('12+log(O/H)')
    plt.axis([0, 4, 8.0, 8.8])
    plt.legend()
    plt.savefig(save_dir + '/' + maps.plateifu + '_sptially_resolved_mzr.pdf')
    

# ------------------------------------------------------------------------------------------------
	 
def run_kinematics(object_name,save_dir):
	if object_name == None:
		object_name = get_query()
	# check if dir exists
	if not os.path.isdir(save_dir):
		os.mkdir(save_dir)
	plot_image(object_name)
	plt.savefig(save_dir + '/' + object_name + '_image.pdf')
	maps = load_maps(object_name)
	plot_maps(maps)
	plt.savefig(save_dir + '/' + object_name + '_maps.pdf')
	plot_features(object_name)
	plt.savefig(save_dir + '/' + object_name + '_features.pdf')
	
def run_metallicity(object_name,save_dir):
	if object_name == None:
		object_name = get_query()
	# check if dir exists
	if not os.path.isdir(save_dir):
		os.mkdir(save_dir)
		
	maps = load_maps2(object_name)
	
	plot_image(object_name)
	plt.savefig(save_dir + '/' + object_name + '_image.pdf')
	plot_bpt(maps,save_dir)
	plot_maps_side_by_side(maps)
	plt.savefig(save_dir + '/' + object_name + '_star_forming_spaxels.pdf')
	plot_mzr(maps,save_dir)