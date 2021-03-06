#To run
#bokeh serve --show core_validation.py

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.interval import PercentileInterval
from skimage import measure

from bokeh.plotting import figure
from bokeh.models.mappers import LinearColorMapper, LogColorMapper
from bokeh.models import ColumnDataSource, Ellipse, Button, DataTable, TableColumn, Slider, \
	HoverTool, BoxEditTool, Div, Select, RadioButtonGroup, RangeSlider, FileInput, TextInput, \
	PrintfTickFormatter, CheckboxGroup, CDSView, IndexFilter
from bokeh.layouts import column, row
from bokeh.models.widgets import Panel, Tabs
from bokeh.transform import factor_cmap
from bokeh.io import curdoc

from read_write_data import data_to_CDS

import numpy as np
import pandas as pd
import os
from datetime import datetime

#Read maps
#########
valfold = './valid_cat/'
filespath = './data/'
time = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

names = []
files = os.listdir(filespath)
for file in files:
	if file[-5:] == ".fits":
		names.append(file[:-5])
names.sort()

username = TextInput(value="Dave_Bowman", title="Username:", width=300)
#if not os.path.isdir(valfold+username.value+'/'):
#    os.makedirs(valfold+username.value+'/')

listcat = []
cats = os.listdir(valfold+username.value+'/')
for file in cats:
	if file[-12:] == "ValidCat.csv":
		listcat.append(file)
listcat.sort(reverse=True)
		
mapdict, dict1, dict2, pdfdict, vdict = data_to_CDS(filespath, names[0]) 
sourcemap = ColumnDataSource(mapdict)
source = ColumnDataSource(dict1)
source2 = ColumnDataSource(dict2)
sourcepdf = ColumnDataSource(pdfdict)
sourcev = ColumnDataSource(vdict)

#Percentil contrast
percent = [95, 98, 99, 99.5, 99.99]

#Contour levels
isomin = 10
isomax = 90
isostep = 5
isonum = int((isomax - isomin)/isostep)
levels = np.linspace(isomin, isomax, isonum)
dlevels = ((np.array(levels) * np.max(sourcemap.data['image'][0]))/100.)
isodict = {'level' : [],
		   'x': [],
		   'y': []}
for dlevel in dlevels.tolist():
	contours = measure.find_contours(sourcemap.data['image'][0], dlevel)
	for contour in contours:
		isodict['level'].append(dlevel)
		isodict['x'].append(contour[:,1].tolist())
		isodict['y'].append(contour[:,0].tolist())
iso = ColumnDataSource(isodict)

#Interaction functions
######################

def selectmap(attr, old, new):
	global mapdict, dict1, dict2, pdfdict, vdict, time
	mapdict, dict1, dict2, pdfdict, vdict = data_to_CDS(filespath, new)
	sourcemap.data = mapdict
	source.data = dict1
	source2.data = dict2
	sourcepdf.data = pdfdict
	sourcev.data = vdict
	time = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
	
	#Contrast
	interval = PercentileInterval(percent[radio_button_group.active])
	vmin, vmax = interval.get_limits(sourcemap.data['image'][0])
	if vmin < 0.:
		vmin = 0.
	color_mapper.low = vmin
	color_mapper.high = vmax
	
	#Contours
	levels = np.linspace(isomin, isomax, isonum)
	dlevels = ((np.array(levels) * np.max(sourcemap.data['image'][0]))/100.)
	newisodict = {'level' : [],
			   'x': [],
			   'y': []}
	for dlevel in dlevels.tolist():
		contours = measure.find_contours(sourcemap.data['image'][0], dlevel)
		for contour in contours:
			newisodict['level'].append(dlevel)
			newisodict['x'].append(contour[:,1].tolist())
			newisodict['y'].append(contour[:,0].tolist())
	iso.data = newisodict
	range_contours.update(value = (isomin,isomax))
	step_slider.update(value = isostep)
	
	plot.plot_height = np.int32(pw *(sourcemap.data['height'][0]/sourcemap.data['width'][0]))
	plot.x_range.start = 0
	plot.x_range.end = sourcemap.data['width'][0]
	plot.y_range.start = 0
	plot.y_range.end = sourcemap.data['height'][0]
	plot.title.text = sourcemap.data['name'][0]
	
	distri.x_range.start = sourcev.data['vmin'][0]
	distri.x_range.end = sourcev.data['vmax'][0]
	distri.y_range.start = 0
	distri.y_range.end = sourcepdf.data['PDF'].max()
	
	range_slider.update(start = sourcev.data['vmin'][0],
					    end = sourcev.data['vmax'][0],
					    value = (vmin,vmax))
	

def select_core(attr, old, new):
	#Change status
	#sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
	#for key in ['ra','dec','x','y','width','height','angles','validation']:
	#	for index in range(len(source.data['ra'])):
	#		sdict[key].append(source.data[key][index])
	for index in new:
		#if sdict['validation'][index] == 'Undefined':
		#	sdict['validation'][index] = 'VALID'
		#	source.data = sdict
		#elif sdict['validation'][index] == 'VALID':
		#	sdict['validation'][index] = 'REJECT'
		#	source.data = sdict
		#else:
		#	sdict['validation'][index] = 'Undefined'
		#	source.data = sdict
		#Change select status button
		if source.data['validation'][index] != 'Modified':
			status_button.active = STATUS.index(source.data['validation'][index])
		else:
			status_button.active = None
		#Update parameters
		ellx.value = source.data['x'][index]
		ellx.start = source.data['x'][index] - 50.
		ellx.end = source.data['x'][index] + 50.
		elly.value = source.data['y'][index]
		elly.start = source.data['y'][index] - 50.
		elly.end = source.data['y'][index] + 50.
		ellw.value = source.data['width'][index]
		ellh.value = source.data['height'][index]
		ella.value = source.data['angles'][index]
	#Save catalog
	#validate = pd.DataFrame.from_dict(sdict)
	#if not os.path.isdir(valfold+username.value+'/'):
	#	os.makedirs(valfold+username.value+'/')
	#validate.to_csv(valfold+username.value+'/'+'{}_{}_{}_ValidCat.csv'
	#				.format(select.value,username.value,time),index=False)
	#if '{}_{}_{}_ValidCat.csv'.format(select.value,username.value,time) not in listcat:
	#	listcat.append('{}_{}_{}_ValidCat.csv'.format(select.value,username.value,time))
	#	select_saved.options = listcat
	#	select_saved.update(value = '{}_{}_{}_ValidCat.csv'
	#						.format(select.value,username.value,time))
		
def select_status(attr, old, new):
	sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
	for key in ['ra','dec','x','y','width','height','angles','validation']:
		for index in range(len(source.data['ra'])):
			sdict[key].append(source.data[key][index])
	for index in source.selected.indices:
		if sdict['validation'][index] != 'Modified':
			sdict['validation'][index] = STATUS[status_button.active]
			source.data = sdict
	#Save catalog
	validate = pd.DataFrame.from_dict(sdict)
	if not os.path.isdir(valfold+username.value+'/'):
		os.makedirs(valfold+username.value+'/')
	validate.to_csv(valfold+username.value+'/'+'{}_{}_{}_ValidCat.csv'
					.format(select.value,username.value,time),index=False)
	if '{}_{}_{}_ValidCat.csv'.format(select.value,username.value,time) not in listcat:
		listcat.append('{}_{}_{}_ValidCat.csv'.format(select.value,username.value,time))
		select_saved.options = listcat
		select_saved.update(value = '{}_{}_{}_ValidCat.csv'
							.format(select.value,username.value,time))

def select_newcore(attr, old, new):
	for index in new:
		ellx.value = source2.data['x'][index]
		ellx.start = source2.data['x'][index] - 50.
		ellx.end = source2.data['x'][index] + 50.
		elly.value = source2.data['y'][index]
		elly.start = source2.data['y'][index] - 50.
		elly.end = source2.data['y'][index] + 50.
		ellw.value = source2.data['width'][index]
		ellh.value = source2.data['height'][index]
		if source2.data['angles'][index] == None:
			ella.value = 0
		else:
			ella.value = source2.data['angles'][index]
			
def reset_core():
	sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
	for key in ['ra','dec','x','y','width','height','angles','validation']:
		for index in range(len(source.data['ra'])):
			sdict[key].append(source.data[key][index])
	for index in source.selected.indices:
		for key in ['ra','dec','x','y','width','height','angles','validation']:
			sdict[key][index] = dict1[key][index]
	source.data = sdict
	#Change select status button
	status_button.active = 0
			
def centering():
	subheigth = 300
	subwidth = np.int32(subheigth *(sourcemap.data['width'][0]/sourcemap.data['height'][0]))
	if len(source.selected.indices) != 0:
		index = source.selected.indices[0]
		plot.x_range.start = source.data['x'][index] - subwidth/2
		plot.x_range.end = source.data['x'][index] + subwidth/2
		plot.y_range.start = source.data['y'][index] - subheigth/2
		plot.y_range.end = source.data['y'][index] + subheigth/2
	elif len(source2.selected.indices) != 0:
		index = source2.selected.indices[0]
		plot.x_range.start = source2.data['x'][index] - subwidth/2
		plot.x_range.end = source2.data['x'][index] + subwidth/2
		plot.y_range.start = source2.data['y'][index] - subheigth/2
		plot.y_range.end = source2.data['y'][index] + subheigth/2
	
def reset_view():
	plot.x_range.start = 0
	plot.x_range.end = sourcemap.data['width'][0]
	plot.y_range.start = 0
	plot.y_range.end = sourcemap.data['height'][0]
	
def save_newcat(attr, old, new):
	NewCat = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[])
	for Nindex in range(len(source2.data['x'])):
		Nx = source2.data['x'][Nindex]
		Ny = source2.data['y'][Nindex]
		HDU = fits.open(filespath+select.value+'.fits')
		header = HDU[0].header
		wcs = WCS(header)
		Nra, Ndec = wcs.all_pix2world(Nx,Ny,1)
		NewCat['ra'].append(Nra)
		NewCat['dec'].append(Ndec)
		for key in ['x','y','width','height','angles']:
			NewCat[key].append(source2.data[key][Nindex])
		newcat = pd.DataFrame.from_dict(NewCat)
		newcat.to_csv(valfold+username.value+'/'+'{}_{}_{}_NewCat.csv'
					  .format(select.value,username.value,time),index=False)
		
def setlistcat(attr, old, new):
	global listcat, dict1, dict2
	source.data = dict1
	source2.data = dict2
	listcat = []
	if not os.path.isdir(valfold+username.value+'/'):
		os.makedirs(valfold+username.value+'/')
	cats = os.listdir(valfold+username.value+'/')
	#print("listcat=",cats)
	for file in cats:
		if file[-12:] == "ValidCat.csv":
			listcat.append(file)
	listcat.sort(reverse=True)
	select_saved.options = listcat
	select_saved.update(value = listcat[0])
	

def save_valid():
	sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
	for key in ['ra','dec','x','y','width','height','angles','validation']:
		for index in range(len(source.data['ra'])):
			sdict[key].append(source.data[key][index])
	validate = pd.DataFrame.from_dict(sdict)
	if not os.path.isdir(valfold+username.value+'/'):
            os.makedirs(valfold+username.value+'/')
	validate.to_csv(valfold+username.value+'/'+'{}_{}_{}_ValidCat.csv'
					.format(select.value,username.value,time),index=False)
	if '{}_{}_{}_ValidCat.csv'.format(select.value,username.value,time) not in listcat:
		listcat.append('{}_{}_{}_ValidCat.csv'.format(select.value,username.value,time))
		select_saved.options = listcat
		select_saved.update(value = '{}_{}_{}_ValidCat.csv'
							.format(select.value,username.value,time))
		
	NewCat = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[])
	for Nindex in range(len(source2.data['x'])):
		Nx = source2.data['x'][Nindex]
		Ny = source2.data['y'][Nindex]
		HDU = fits.open(filespath+select.value+'.fits')
		header = HDU[0].header
		wcs = WCS(header)
		Nra, Ndec = wcs.all_pix2world(Nx,Ny,1)
		NewCat['ra'].append(Nra)
		NewCat['dec'].append(Ndec)
		for key in ['x','y','width','height','angles']:
			NewCat[key].append(source2.data[key][Nindex])
		newcat = pd.DataFrame.from_dict(NewCat)
		newcat.to_csv(valfold+username.value+'/'+'{}_{}_{}_NewCat.csv'
					  .format(select.value,username.value,time),index=False)
		
def load_valid():
	global time
	#validfile = './valid_cat/{}_{}_ValidCat.csv'.format(select.value,time)
	#newcatfile = './valid_cat/{}_{}_NewCat.csv'.format(select.value,time)
	file_input = select_saved.value
	time = file_input[-32:-13]
	validfile = valfold+username.value+'/'+file_input
	validate = pd.read_csv(validfile)
	cdict = pd.DataFrame.to_dict(validate,orient='list')
	source.data = cdict

	newcatfile = valfold+username.value+'/'+file_input[:-12]+'NewCat.csv'
	if os.path.exists(newcatfile) == True:
		fields = ['x','y','width','height','angles']
		newcatfile = valfold+username.value+'/'+file_input[:-12]+'NewCat.csv'
		newcat = pd.read_csv(newcatfile, usecols=fields)
		expdict = pd.DataFrame.to_dict(newcat,orient='list')
		source2.data = expdict

'''def mapscale(attr, old, new):
	global color_mapper
	vmin = range_slider.value[0]
	vmax = range_slider.value[1]
	if scaling_button.active == 0:
		color_mapper = LinearColorMapper(palette="Inferno256",low=vmin,high=vmax)
	else:
		color_mapper = LogColorMapper(palette="Inferno256",low=vmin,high=vmax)
	bkg.glyph.color_mapper = color_mapper'''
		
def percentil(attr, old, new):
	#global color_mapper
	interval = PercentileInterval(percent[radio_button_group.active])
	vmin, vmax = interval.get_limits(sourcemap.data['image'][0])
	if vmin < 0.:
		vmin = 0.
	color_mapper.low = vmin
	color_mapper.high = vmax
	#bkg.glyph.color_mapper = color_mapper
	range_slider.update(value = (vmin, vmax))
	newvdict = {'vmin': [vmin,vmin],
				'vmax': [vmax,vmax],
			    'y': sourcev.data['y']}
	sourcev.data = newvdict
	
def adjustcontrast(attr, old, new):
	#global color_mapper
	newvdict = {'vmin': [range_slider.value[0],range_slider.value[0]],
				'vmax': [range_slider.value[1],range_slider.value[1]],
			    'y': sourcev.data['y']}
	sourcev.data = newvdict
	color_mapper.low = range_slider.value[0]
	color_mapper.high = range_slider.value[1]
	#bkg.glyph.color_mapper = color_mapper
	
def Xelliparam(attr, old, new):
	if not (len(source.selected.indices) == 0):
		sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
		for key in ['ra','dec','x','y','width','height','angles','validation']:
			for index in range(len(source.data['ra'])):
				sdict[key].append(source.data[key][index])
		for index in source.selected.indices:
			sdict['x'][index] = ellx.value
			HDU = fits.open(filespath+select.value+'.fits')
			header = HDU[0].header
			wcs = WCS(header)
			Nra, Ndec = wcs.all_pix2world(ellx.value,elly.value,1)
			sdict['ra'][index] = Nra
			sdict['dec'][index] = Ndec
			source.data = sdict
		if (abs(sdict['x'][index] - dict1['x'][index]) > 0.00001) or \
		   (abs(sdict['y'][index] - dict1['y'][index]) > 0.00001):
			sdict['validation'][index] = 'Modified'
			status_button.active = None
	if not (len(source2.selected.indices) == 0):
		NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
		for Mindex in range(len(source2.data['x'])):
			for key in ['x','y','width','height','angles']:
					NewCoreCat[key].append(source2.data[key][Mindex])
		for index in source2.selected.indices:
			NewCoreCat['x'][index] = ellx.value
			source2.data = NewCoreCat
			
def Yelliparam(attr, old, new):
	if not (len(source.selected.indices) == 0):
		sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
		for key in ['ra','dec','x','y','width','height','angles','validation']:
			for index in range(len(source.data['ra'])):
				sdict[key].append(source.data[key][index])
		for index in source.selected.indices:
			sdict['y'][index] = elly.value
			HDU = fits.open(filespath+select.value+'.fits')
			header = HDU[0].header
			wcs = WCS(header)
			Nra, Ndec = wcs.all_pix2world(ellx.value,elly.value,1)
			sdict['ra'][index] = Nra
			sdict['dec'][index] = Ndec
			source.data = sdict
		if (abs(sdict['x'][index] - dict1['x'][index]) > 0.00001) or \
		   (abs(sdict['y'][index] - dict1['y'][index]) > 0.00001):
			sdict['validation'][index] = 'Modified'
			status_button.active = None
	if not (len(source2.selected.indices) == 0):
		NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
		for Mindex in range(len(source2.data['x'])):
			for key in ['x','y','width','height','angles']:
					NewCoreCat[key].append(source2.data[key][Mindex])
		for index in source2.selected.indices:
			NewCoreCat['y'][index] = elly.value
			source2.data = NewCoreCat
	
def Welliparam(attr, old, new):
	if not (len(source.selected.indices) == 0):
		sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
		for key in ['ra','dec','x','y','width','height','angles','validation']:
			for index in range(len(source.data['ra'])):
				sdict[key].append(source.data[key][index])
		for index in source.selected.indices:
			sdict['width'][index] = ellw.value
			source.data = sdict
		if abs(sdict['width'][index] - dict1['width'][index]) > 0.001:
			sdict['validation'][index] = 'Modified'
			status_button.active = None
	if not (len(source2.selected.indices) == 0):
		NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
		for Mindex in range(len(source2.data['x'])):
			for key in ['x','y','width','height','angles']:
					NewCoreCat[key].append(source2.data[key][Mindex])
		for index in source2.selected.indices:
			NewCoreCat['width'][index] = ellw.value
			source2.data = NewCoreCat
		
def Helliparam(attr, old, new):
	if not (len(source.selected.indices) == 0):
		sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
		for key in ['ra','dec','x','y','width','height','angles','validation']:
			for index in range(len(source.data['ra'])):
				sdict[key].append(source.data[key][index])
		for index in source.selected.indices:
			sdict['height'][index] = ellh.value
			source.data = sdict
		if abs(sdict['height'][index] - dict1['height'][index]) > 0.001:
			sdict['validation'][index] = 'Modified'
			status_button.active = None
	if not (len(source2.selected.indices) == 0):
		NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
		for Mindex in range(len(source2.data['x'])):
			for key in ['x','y','width','height','angles']:
					NewCoreCat[key].append(source2.data[key][Mindex])
		for index in source2.selected.indices:
			NewCoreCat['height'][index] = ellh.value
			source2.data = NewCoreCat
		
def Aelliparam(attr, old, new):
	if not (len(source.selected.indices) == 0):
		sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
		for key in ['ra','dec','x','y','width','height','angles','validation']:
			for index in range(len(source.data['ra'])):
				sdict[key].append(source.data[key][index])
		for index in source.selected.indices:
			sdict['angles'][index] = ella.value
			source.data = sdict
		if abs(sdict['angles'][index] - dict1['angles'][index]) > 0.001:
			sdict['validation'][index] = 'Modified'
			status_button.active = None
	if not (len(source2.selected.indices) == 0):
		NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
		for Mindex in range(len(source2.data['x'])):
			for key in ['x','y','width','height','angles']:
					NewCoreCat[key].append(source2.data[key][Mindex])
		for index in source2.selected.indices:
			NewCoreCat['angles'][index] = ella.value
			source2.data = NewCoreCat
			
def adjustcontours(attr, old, new):
	newisodict = {'level' : [],
				  'x': [],
				  'y': []}
	isonum = int((new[1] - new[0])/step_slider.value)
	levels = np.linspace(new[0], new[1], isonum)
	dlevels = ((np.array(levels) * np.max(sourcemap.data['image'][0]))/100.)
	for dlevel in dlevels.tolist():
		contours = measure.find_contours(sourcemap.data['image'][0], dlevel)
		for contour in contours:
			newisodict['level'].append(dlevel)
			newisodict['x'].append(contour[:,1].tolist())
			newisodict['y'].append(contour[:,0].tolist())
	iso.data = newisodict
	
def stepcontours(attr, old, new):
	newisodict = {'level' : [],
				  'x': [],
				  'y': []}
	isonum = int((range_contours.value[1] - range_contours.value[0])/new)
	levels = np.linspace(range_contours.value[0], range_contours.value[1], isonum)
	dlevels = ((np.array(levels) * np.max(sourcemap.data['image'][0]))/100.)
	for dlevel in dlevels.tolist():
		contours = measure.find_contours(sourcemap.data['image'][0], dlevel)
		for contour in contours:
			newisodict['level'].append(dlevel)
			newisodict['x'].append(contour[:,1].tolist())
			newisodict['y'].append(contour[:,0].tolist())
	iso.data = newisodict
			
#Plot image and Table
#####################
	
ph = 500
pw = np.int32(ph *(sourcemap.data['width'][0]/sourcemap.data['height'][0]))

select = Select(title="Region:", value=names[0], options=names)
select.on_change("value", selectmap)
plot = figure(plot_width = pw, plot_height = ph, x_range=(0,sourcemap.data['width'][0]), y_range=(0,sourcemap.data['height'][0]), match_aspect = True,tools=["tap,wheel_zoom,pan"],title=names[0])
plot.axis.visible = False

interval = PercentileInterval(99.5)
vmin, vmax = interval.get_limits(sourcemap.data['image'][0])
if vmin < 0.:
	vmin = 0.
color_mapper = LinearColorMapper(palette="Inferno256",low=vmin,high=vmax)

bkg = plot.image(image='image',
           dh='height', dw='width', x=0, y=0, source=sourcemap, color_mapper=color_mapper)

FSTATUS = ['Undefined', 'VALID','REJECT','Modified']
COLOR = ['gray','red','blue','red']
renderer = plot.ellipse(x="x", y="y", width="width", height="height", angle="angles", line_color=factor_cmap('validation', COLOR, FSTATUS),fill_alpha=0,source=source,legend_field='validation')

#Trace first position and size for new cores
expert1 = plot.rect('x', 'y', 'width', 'height', source=source2, line_alpha=0, fill_alpha=0)

#Adjust ellipse in rectangle
expert2 = plot.ellipse(x="x", y="y", width="width", height="height", angle="angles", line_color='green',fill_alpha=0,source=source2,legend_label='New cores')

#Contours
#################
#contour_text = Div(text="""Levels:"""+np.str_(np.unique(np.array(iso.data['level']))), width=400)
plot.multi_line(xs="x", ys="y", color='#02FF34',source=iso, line_width=1.0, legend_label='Contours')
contours_title = Div(text="""<h2>Plot contours</h2>""", width=400)

range_contours = RangeSlider(start=0, end=100, value=(levels[0], levels[-1]), step=.1, title="Min/Max % isocontours", width = 400)
range_contours.on_change("value",adjustcontours)

step_slider = Slider(start=1, end=50, value= isostep, step=1, title="Contour's step", max_width = np.int32(pw/2))
step_slider.on_change("value",stepcontours)
		
plot.legend.location = "top_right"
plot.legend.click_policy= "hide"

#Tables & validation button
###########################

STATUS = ['Undefined', 'VALID','REJECT']
status_button = RadioButtonGroup(labels=STATUS, active=0, width=300)
status_button.on_change("active",select_status)

columns = [
        TableColumn(field="ra", title="ra"),
        TableColumn(field="dec", title="dec"),
		TableColumn(field="validation", title="status")
    ]

data_table = DataTable(source=source, columns=columns, width=300, height=ph)

Ncolumns = [
        TableColumn(field="width", title="width"),
        TableColumn(field="height", title="height"),
		TableColumn(field="angles", title="angle")
    ]

Ndata_table = DataTable(source=source2, columns=Ncolumns, width=300, height=ph)

#Plot PDF
#####################

contrast_title = Div(text="""<h2>Adjust contrast</h2>""", width=400)

distri = figure(plot_width = 400, plot_height = 250, 
				x_range=(sourcev.data['vmin'][0],sourcev.data['vmax'][0]), 
				y_range=(0,sourcepdf.data['PDF'].max()), match_aspect = True, 
				tools=["box_zoom,pan,reset"],title="Map's histogram")
distri.line(x='PDFbin', y='PDF',line_width=1.5,source=sourcepdf)
range_slider = RangeSlider(start=sourcev.data['vmin'][0], end=sourcepdf.data['PDFbin'].max(), value=(sourcev.data['vmin'][0],sourcev.data['vmax'][0]), step=.1, title="Adjust contrast", width = 400, format=PrintfTickFormatter(format="%.2e"))
range_slider.on_change("value",adjustcontrast)
distri.line(x='vmin',y='y',source=sourcev,color='orange',line_width=1.5)
distri.line(x='vmax',y='y',source=sourcev,color='orange',line_width=1.5)
distri.xaxis.formatter=PrintfTickFormatter(format="%e")
distri.yaxis.formatter=PrintfTickFormatter(format="%e")


#Side Column : Description, Interactions and Buttons
#########################

text = Div(text="""<h1><img src="https://ipag.osug.fr/~robitaij/figures/logo_crop.png" 
            width="149" height="123" align="left" margin=20px>
			LHYRICA Project</h1></br>
This application allows you to validate cores catalogs in order to train LHYRICA's neural network.
""",
width=400)

#username = TextInput(value="Dave_Bowman", title="Username:",max_width = np.int32(pw/2))
username.on_change("value",setlistcat)

'''LABELSCALES = ["Linear", "Log"]
scaling_button = RadioButtonGroup(labels=LABELSCALES, active=0)
scaling_button.on_change("active",mapscale)'''

LABELS = ["95%", "98%", "99%","99.5%","99.99%"]
radio_button_group = RadioButtonGroup(labels=LABELS, active=3)
radio_button_group.on_change("active",percentil)
#slider = Slider(start=50, end=100, value= 99.5, step=0.01, title="Percentile Interval (image contrast)", max_width = np.int32(pw/2))
#slider.on_change("value",percentil)

expert1.data_source.selected.on_change("indices", select_newcore)
ellx = Slider(start=2, end=100, value= 20, step=0.5, title="Ellipse x position", max_width = np.int32(pw/2))
elly = Slider(start=2, end=100, value= 20, step=0.5, title="Ellipse y position", max_width = np.int32(pw/2))
ellw = Slider(start=2, end=100, value= 20, step=0.1, title="Ellipse width", max_width = np.int32(pw/2))
ellh = Slider(start=2, end=100, value= 20, step=0.1, title="Ellipse height", max_width = np.int32(pw/2))
ella = Slider(start=-np.pi/2, end=np.pi/2, value= 0, step=0.1, title="Ellipse angle", max_width = np.int32(pw/2))
ellx.on_change("value",Xelliparam)
elly.on_change("value",Yelliparam)
ellw.on_change("value",Welliparam)
ellh.on_change("value",Helliparam)
ella.on_change("value",Aelliparam)

renderer.data_source.selected.on_change("indices", select_core)

RCbutton = Button(label="Reset core",button_type="default",max_width = np.int32(pw/2))
RCbutton.on_click(reset_core)

RVbutton = Button(label="Reset view",button_type="default",max_width = np.int32(pw/2))
RVbutton.on_click(reset_view)

Sbutton = Button(label="Save Cat",button_type="success",max_width = np.int32(pw/2))
Sbutton.on_click(save_valid)

selval=listcat[0] if len(listcat) != 0  else ''
select_saved = Select(title="Select saved catalog:", value=selval, options=listcat)
Lbutton = Button(label="Load Cat",button_type="primary",max_width = np.int32(pw/2))
Lbutton.on_click(load_valid)

boxedit = BoxEditTool(renderers=[expert1], empty_value = 0.)
expert2.data_source.on_change('data',save_newcat)

center = Button(label="center object",button_type="default",max_width = np.int32(pw/2))
center.on_click(centering)

tooltips = [('index','$index'),('ra', '@ra'), ('dec', '@dec')]
plot.add_tools(HoverTool(tooltips=tooltips),boxedit)

#Tutorial
tuto = Div(text="""<h2>Memo</h2>
Enter a username to create a catalogue repository and keep the same username for later sessions.</br>
 </br>
The catalog on the left panel is linked to the cores displayed on the map. This means that cores can be validated/rejected with a click on the catalog or on the map.</br>
 </br>
Ellipses and contours can be activated/hidden by clicking in the plot's legend.
 </br>
<img src="https://docs.bokeh.org/en/latest/_images/WheelZoom.png"> Zoom with mouse wheel or track pad</br>
<img src="https://docs.bokeh.org/en/latest/_images/Tap.png"> Valid/reject catalog's cores</br>
<img src="https://docs.bokeh.org/en/latest/_images/BoxEdit.png"> Shift+click to add new cores</br>
<img src="https://docs.bokeh.org/en/latest/_images/Hover.png"> Activate/deactivate core's tag</br>""", width=400)

#Dashboard Display

tab1 = Panel(child=data_table, title='catalog')
tab2 = Panel(child=Ndata_table, title='new cores')
tabs = Tabs(tabs=[tab1,tab2])
present = column(text,username,status_button,tabs)

#Orion configuration
#mapint = column(plot,row(column(ellx, elly, ellw, ellh, ella),column(RCbutton,center,RVbutton,select,select_saved,Lbutton,Sbutton)))
#contrast = column(contrast_title, radio_button_group, distri, range_slider, tuto)
#panel = row(present, mapint, contrast)

#Taurus configuration
contrast = column(radio_button_group, distri, range_slider, contours_title, range_contours, step_slider)
mapint = column(plot,row(column(ellx, elly, ellw, ellh, ella),column(RCbutton,center,RVbutton,select,select_saved,Lbutton,Sbutton),contrast,tuto))
panel = row(present, mapint)

curdoc().add_root(panel)