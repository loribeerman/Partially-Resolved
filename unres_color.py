import numpy as np
import matplotlib.pyplot as plt
import pyfits
plt.ion()


def get_phot(phot_file):
    '''read in photometry file for cluster'''
    phot = pyfits.open(phot_file)
    phot_data = phot[1].data
    phot_275 = phot_data.field('F275W_VEGA')
    phot_336 = phot_data.field('F336W_VEGA')
    phot_475 = phot_data.field('F475W_VEGA')
    phot_814 = phot_data.field('F814W_VEGA')
    phot_110 = phot_data.field('F110W_VEGA')
    phot_160 = phot_data.field('F160W_VEGA')
    return phot_275, phot_336, phot_475, phot_814, phot_110, phot_160


def mag_to_flux(mag):
    '''convert mag to flux'''
    flux = 10.**((mag)/(-2.5))
    return flux


def flux_to_mag(flux):
    '''convert flux to mag'''
    mag = -2.5*np.log10(flux)
    return mag


#read in table containing integrated phot
clust_file = pyfits.open('apdata_cluster_v2.fits')
clust = clust_file[1].data
clust_id = clust.field('id')
int_mag_275 = clust.field('mag275')
int_mag_336 = clust.field('mag336')
int_mag_475 = clust.field('mag475')
int_mag_814 = clust.field('mag814')
int_mag_110 = clust.field('mag110')
int_mag_160 = clust.field('mag160')


#define cutoffmag
cutoffmag = 24.47-3.0

#read in phot for individual stars
ap_id = np.genfromtxt('done_list_1215.txt', dtype=None)

#define arrays
unres_color_48 = np.zeros(len(ap_id))
unres_color_38 = np.zeros(len(ap_id))
unres_color_34 = np.zeros(len(ap_id))
int_color_38 = np.zeros(len(ap_id))
int_color_48 = np.zeros(len(ap_id))
int_color_34 = np.zeros(len(ap_id))
unres_mag_275 = np.zeros(len(ap_id))
unres_mag_336 = np.zeros(len(ap_id))
unres_mag_475 = np.zeros(len(ap_id))
unres_mag_814 = np.zeros(len(ap_id))
unres_mag_110 = np.zeros(len(ap_id))
unres_mag_160 = np.zeros(len(ap_id))
int_mag_275_keep = np.zeros(len(ap_id))
int_mag_336_keep = np.zeros(len(ap_id))
int_mag_475_keep = np.zeros(len(ap_id))
int_mag_814_keep = np.zeros(len(ap_id))
int_mag_110_keep = np.zeros(len(ap_id))
int_mag_160_keep = np.zeros(len(ap_id))
int_id = np.zeros(len(ap_id))


#for i in range(1,2):
for i in range(len(ap_id)):

    name = ap_id[i]
    int_id[i] = int(name.split('ap')[1])

    try:
	#read 6 band photometry
        phot_275, phot_336, phot_475, phot_814, phot_110, phot_160 = get_phot('../April2014_runs/' +name+ '/' +name+ '_phot.fits')

	#find stars brighter than mag limit
	bright = np.where(phot_814 <= cutoffmag)

	#find total flux brighter than mag limit
	bright_flux_275 = np.sum(mag_to_flux(phot_275[bright]))
	bright_flux_336 = np.sum(mag_to_flux(phot_336[bright]))
	bright_flux_475 = np.sum(mag_to_flux(phot_475[bright]))
	bright_flux_814 = np.sum(mag_to_flux(phot_814[bright]))
	bright_flux_110 = np.sum(mag_to_flux(phot_110[bright]))
	bright_flux_160 = np.sum(mag_to_flux(phot_160[bright]))

	#calc int flux
	clust_place = np.where(clust_id == int_id[i])
	int_flux_275 = mag_to_flux(int_mag_275[clust_place])
	int_flux_336 = mag_to_flux(int_mag_336[clust_place])
	int_flux_475 = mag_to_flux(int_mag_475[clust_place])
	int_flux_814 = mag_to_flux(int_mag_814[clust_place])
	int_flux_110 = mag_to_flux(int_mag_110[clust_place])
	int_flux_160 = mag_to_flux(int_mag_160[clust_place])

	#subtract bright flux from int flux to get 'unres' flux
	if int_flux_275 >= bright_flux_275:  
            unres_flux_275 = int_flux_275 - bright_flux_275
	else:
	    unres_flux_275 = -1
	if int_flux_336 >= bright_flux_336:  
	    unres_flux_336 = int_flux_336 - bright_flux_336
	else:
	    unres_flux_336 = -1
	if int_flux_475 >= bright_flux_475:  
	    unres_flux_475 = int_flux_475 - bright_flux_475
	else:
	    unres_flux_475 = -1
	if int_flux_814 >= bright_flux_814:  
	    unres_flux_814 = int_flux_814 - bright_flux_814
	else:
	    unres_flux_814 = -1
	if int_flux_110 >= bright_flux_110:  
	    unres_flux_110 = int_flux_110 - bright_flux_110
	else:
	    unres_flux_110 = -1
	if int_flux_160 >= bright_flux_160:  
	    unres_flux_160 = int_flux_160 - bright_flux_160
	else:
	    unres_flux_160 = -1

	#calc 'unres' mag
	unres_mag_275[i] = flux_to_mag(unres_flux_275)
	unres_mag_336[i] = flux_to_mag(unres_flux_336)
	unres_mag_475[i] = flux_to_mag(unres_flux_475)
	unres_mag_814[i] = flux_to_mag(unres_flux_814)
	unres_mag_110[i] = flux_to_mag(unres_flux_110)
	unres_mag_160[i] = flux_to_mag(unres_flux_160)

	#calc unres color
	unres_color_48[i] = unres_mag_475[i] - unres_mag_814[i]
	unres_color_34[i] = unres_mag_336[i] - unres_mag_475[i]
	unres_color_38[i] = unres_mag_336[i] - unres_mag_814[i]

	#find int color
	if int_mag_475[clust_place].size == 1:  
	    int_color_48[i] = int_mag_475[clust_place] - int_mag_814[clust_place]
            int_mag_475_keep[i] = int_mag_475[clust_place]
            int_mag_814_keep[i] = int_mag_814[clust_place]
	else:
	    int_color_48[i] = 99.99

	if int_mag_336[clust_place].size == 1:  
	    int_color_34[i] = int_mag_336[clust_place] - int_mag_475[clust_place]
            int_mag_336_keep[i] = int_mag_336[clust_place]
            int_mag_275_keep[i] = int_mag_275[clust_place]
	else:
	    int_color_34[i] = 99.99

	if int_mag_336[clust_place].size == 1:  
	    int_color_38[i] = int_mag_336[clust_place] - int_mag_814[clust_place]
            int_mag_110_keep[i] = int_mag_110[clust_place]
            int_mag_160_keep[i] = int_mag_160[clust_place]
	else:
	    int_color_38[i] = 99.99


    except IOError:
	print name +' phot does not exist'


#save results to file
dt = np.dtype([('ap_id', 'i'), ('int_mag_275', 'd'), ('int_mag_336', 'd'), ('int_mag_475', 'd'), ('int_mag_814', 'd'), ('int_mag_110', 'd'), ('int_mag_160', 'd'), ('int_color_34', 'd'), ('int_color_38', 'd'), ('int_color_48', 'd'), ('unres_mag_275', 'd'), ('unres_mag_336', 'd'), ('unres_mag_475', 'd'), ('unres_mag_814', 'd'), ('unres_mag_110', 'd'), ('unres_mag_160', 'd'), ('unres_color_34', 'd'), ('unres_color_38', 'd'), ('unres_color_48', 'd')])

a = np.zeros(len(ap_id), dt)

a['ap_id'] = int_id
a['int_mag_275'] = int_mag_275_keep
a['int_mag_336'] = int_mag_336_keep
a['int_mag_475'] = int_mag_475_keep
a['int_mag_814'] = int_mag_814_keep
a['int_mag_110'] = int_mag_110_keep
a['int_mag_160'] = int_mag_160_keep
a['int_color_34'] = int_color_34
a['int_color_38'] = int_color_38
a['int_color_48'] = int_color_48
a['unres_mag_275'] = unres_mag_275
a['unres_mag_336'] = unres_mag_336
a['unres_mag_475'] = unres_mag_475
a['unres_mag_814'] = unres_mag_814
a['unres_mag_110'] = unres_mag_110
a['unres_mag_160'] = unres_mag_160
a['unres_color_34'] = unres_color_34
a['unres_color_38'] = unres_color_38
a['unres_color_48'] = unres_color_48

np.savetxt('unres_results512.txt', a, '%15s')


