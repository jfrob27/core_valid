from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import numpy as np
from turbustat.statistics import PDF
from astropy.visualization.interval import PercentileInterval

def data_to_CDS(filespath, name):
	#Read fits
	HDU = fits.open(filespath+name+'.fits')
	image = HDU[0].data
	header = HDU[0].header
	
	interval = PercentileInterval(99.5)
	vmin, vmax = interval.get_limits(image)
	
	mapdict = {'name': [name],
			   'image': [image],
			   'width': [image.shape[1]],
			   'height': [image.shape[0]]}
	
	#Compute PDF
	pdf = PDF(image, bins=None)
	pdf.run(verbose=False, do_fit=False);
	
	pdfdict = {'PDFbin': pdf.bins,
			   'PDF': pdf.pdf}
	
	vdict = {'vmin': [vmin,vmin],
			 'vmax': [vmax,vmax],
			 'y': [0,pdf.pdf.max()]}
	
	#Read catalogue
	cores = pd.read_csv(filespath+name+'.csv')
	if 'WCS_ACOOR' in list(cores.columns):
		ra = cores['WCS_ACOOR'].to_numpy()
		dec = cores['WCS_DCOOR'].to_numpy()
		wcs = WCS(header)
		reso = header['CDELT2'] * 3600.
		x, y = wcs.all_world2pix(ra,dec,1)
		w = cores['AFWHM01'].to_numpy()/reso
		h = cores['BFWHM01'].to_numpy()/reso
		angles = cores['THETA01'].to_numpy()-90
		angles = (angles * np.pi)/180.
	else:
		ra = cores['ra'].to_numpy()
		dec = cores['dec'].to_numpy()
		wcs = WCS(header)
		reso = header['CDELT2'] * 3600.
		x, y = wcs.all_world2pix(ra,dec,1)
		w = cores['width'].to_numpy()/reso
		h = cores['height'].to_numpy()/reso
		angles = cores['angles'].to_numpy()-90
		angles = (angles * np.pi)/180.	
	catsz = ra.size
	validation = np.repeat('Undefined',catsz)

	dict1 = dict(ra=ra,dec=dec,x=x,y=y,width=w,height=h,angles=angles,validation=validation)
	
	#Catalogue of new added cores
	dict2 = dict(x=[],y=[],width=[],height=[],angles=[])
	
	return mapdict, dict1, dict2, pdfdict, vdict