import random as random
import numpy as np

'''computes chi2 min for one cluster
based off mf.pro idl code by Morgan Fouesneau'''

def getModel(model_flux):
    '''get appropriate model, model units in flux
    model shape = (number of ages, number of filters)'''

    #load model
    #model_flux = np.loadtxt(model)

    #convert model flux into magnitudes
    model =  -2.5*np.log10(model_flux)
    
    return model


def getData(data_flux):
    '''get data from file, data units in flux
    data shape = (number of filters)'''

    #load data
    #data_flux = np.loadtxt(data)

    obs =  -2.5*np.log10(data_flux)

    return obs


def computeMass(obs, err, mod):
    '''compute mass scaled from flux'''

    ferr = 10**(-0.4*obs)*(1.-10**(-0.4*err)) 		#sigmas in flux
    fobs = 10**(-0.4*obs)			 	#input flux
    fmod = 10**(-0.4*mod)	 	         	#model flux
    mass = (np.sum(fobs*fmod/ferr**2. )) / np.sum((fmod/ferr)**2.)
    
    return mass


def compute_chi2(obs, err, mod, mass):
    '''compute chi^2'''
	
    dfree = np.float(len(obs)-1.) 			#degree of freedom
	
    ferr  = 10**(-0.4*obs)*(1.-10**(-0.4*err)) 		#sigmas in flux
    fobs  = 10**(-0.4*obs)			 	#input flux
    fmod  = mass*10**(-0.4*mod)		 		#model flux

    val = np.sum(((fobs - fmod)**2) / (ferr)**2 )/(dfree)
    
    return val


def getAvVect():
    '''return extinction values from Cardelli extinction law'''			
	
    R475 = 1.19119
    R814 = 0.60593
    R160 = 0.20443
    R336 = 1.65798
    Avmag = np.array([ R336, R475, R814, R160 ])

    return Avmag


def find_min(model, Av_arr, mass_arr, chi2_arr):
    '''find minimum chi2 and model and Av for minimum chi2'''

    chi2_min = np.min(chi2_arr)
    min_place = np.where(chi2_arr == chi2_min)
    Av_min = Av_arr[min_place[0]]
    model_min = model[min_place[1]]
    mass_min = mass_arr[min_place]

    return [chi2_min, Av_min, model_min, mass_min]


#m='model_filename.blah'
m=np.random.rand(10,4)
clust_model = getModel(m)

#d='data_filename.blah'
d=np.random.rand(4)
clust_obs = getData(d)

#create a 50 element array from 0 to 3.43, spaced by 0.07
nAv = 50
dAv = 0.07
Av_arr = np.arange(0, nAv*dAv, dAv)
Avmag = getAvVect()

err = 0.05  #???

nModels= clust_model.shape[0]

#define arrays to be filled in
mass_arr = np.zeros([nAv, nModels])
chi2_arr = np.zeros([nAv, nModels])

#loop all possible Av
for iAv in range(nAv):
    
    #loop over all age models
    for kMod in range(nModels):

        this_mod = clust_model[kMod]
        dereddened_obs = clust_obs - (Av_arr[iAv] * Avmag)
        mass_arr[iAv, kMod] = computeMass(dereddened_obs, err, this_mod)
        chi2_arr[iAv, kMod] = compute_chi2(dereddened_obs, err, this_mod, mass_arr[iAv, kMod])

#find min values for this cluster
chi2_min, Av_min, model_min, mass_min = find_min(clust_model,  Av_arr, mass_arr, chi2_arr)
print chi2_min, Av_min, model_min, mass_min


