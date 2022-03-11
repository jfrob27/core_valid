from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import numpy as np

def data_to_CDS(filespath, name):
	#Read fits
	HDU = fits.open(filespath+name+'.fits')
	image = HDU[0].data
	header = HDU[0].header
	
	mapdict = {'name': [name],
			   'image': [image],
			   'width': [image.shape[1]],
			   'height': [image.shape[0]]}
	#Read catalogue
	cores = pd.read_csv(filespath+name+'.csv',sep=';')
	ra = cores['WCS_ACOOR'].to_numpy()
	dec = cores['WCS_DCOOR'].to_numpy()
	wcs = WCS(header)
	reso = header['CDELT2'] * 3600.
	x, y = wcs.all_world2pix(ra,dec,1)
	w = cores['AFWHM01'].to_numpy()/reso
	h = cores['BFWHM01'].to_numpy()/reso
	angles = cores['THETA01'].to_numpy()-90
	angles = (angles * np.pi)/180.
	catsz = ra.size
	validation = np.repeat('Undefined',catsz)
	
	dict1 = dict(ra=ra,dec=dec,x=x,y=y,width=w,height=h,angles=angles,validation=validation)
	
	#Catalogue of new added cores
	dict2 = dict(x=[],y=[],width=[],height=[],angles=[])
	
	return mapdict, dict1, dict2