#To run
#bokeh serve --show core_validation.py

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.interval import PercentileInterval

from bokeh.plotting import figure
from bokeh.models.mappers import LinearColorMapper
from bokeh.models import ColumnDataSource, Ellipse, Button, DataTable, TableColumn, Slider, HoverTool, BoxEditTool, Div, Select, CheckboxGroup
from bokeh.layouts import column, row
from bokeh.models.widgets import Panel, Tabs
from bokeh.transform import factor_cmap
from bokeh.io import curdoc

from read_write_data import data_to_CDS

import numpy as np
import pandas as pd
import os

#Read maps
#########
valfold = './valid_cat/'
filespath = './data/'
files = os.listdir(filespath)

names = []
for file in files:
	if file[-5:] == ".fits":
		names.append(file[:-5])
		
mapdict, dict1, dict2 = data_to_CDS(filespath, names[0]) 
sourcemap = ColumnDataSource(mapdict)
source = ColumnDataSource(dict1)
source2 = ColumnDataSource(dict2)

#Interaction functions
######################

def selectmap(attr, old, new):
	mapdict, dict1, dict2 = data_to_CDS(filespath, new)
	sourcemap.data = mapdict
	source.data = dict1
	source2.data = dict2
	
	plot.plot_height = np.int32(pw *(sourcemap.data['height'][0]/sourcemap.data['width'][0]))
	plot.x_range.start = 0
	plot.x_range.end = sourcemap.data['width'][0]
	plot.y_range.start = 0
	plot.y_range.end = sourcemap.data['height'][0]
	plot.title.text = sourcemap.data['name'][0]

def select_core(attr, old, new):
	sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
	for key in ['ra','dec','x','y','width','height','angles','validation']:
		for index in range(len(source.data['ra'])):
			sdict[key].append(source.data[key][index])
	for index in new:
		if sdict['validation'][index] == 'Undefined':
			sdict['validation'][index] = 'VALID'
			source.data = sdict
		elif sdict['validation'][index] == 'VALID':
			sdict['validation'][index] = 'REJECT'
			source.data = sdict
		else:
			sdict['validation'][index] = 'Undefined'
			source.data = sdict
	validate = pd.DataFrame.from_dict(sdict)
	#valfold = './valid_cat/'
	if not os.path.isdir(valfold):
            os.makedirs(valfold)
	validate.to_csv(valfold+'{}_ValidCat.csv'.format(select.value),index=False)

def select_newcore(attr, old, new):
	for index in new:
		ellw.value = source2.data['width'][index]
		ellh.value = source2.data['height'][index]
		if source2.data['angles'][index] == None:
			ella.value = 0
		else:
			ella.value = source2.data['angles'][index]
			
def centering():
	subheigth = 300
	subwidth = np.int32(subheigth *(sourcemap.data['width'][0]/sourcemap.data['height'][0]))
	for index in source.selected.indices:
		plot.x_range.start = source.data['x'][index] - subwidth/2
		plot.x_range.end = source.data['x'][index] + subwidth/2
		plot.y_range.start = source.data['y'][index] - subheigth/2
		plot.y_range.end = source.data['y'][index] + subheigth/2
	
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
		newcat.to_csv(valfold+'{}_NewCat.csv'.format(select.value),index=False)

def save_valid():
	sdict = dict(ra=[],dec=[],x=[],y=[],width=[],height=[],angles=[],validation=[])
	for key in ['ra','dec','x','y','width','height','angles','validation']:
		for index in range(len(source.data['ra'])):
			sdict[key].append(source.data[key][index])
	#for index in range(catsz):
	validate = pd.DataFrame.from_dict(sdict)
	#valfold = './valid_cat/'
	if not os.path.isdir(valfold):
            os.makedirs(valfold)
	validate.to_csv(valfold+'{}_ValidCat.csv'.format(select.value),index=False)
		
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
		newcat.to_csv(valfold+'{}_NewCat.csv'.format(select.value),index=False)
		
def load_valid():
	#global cdict, expdict
	validfile = './valid_cat/{}_ValidCat.csv'.format(select.value)
	newcatfile = './valid_cat/{}_NewCat.csv'.format(select.value)
	if os.path.exists(validfile):
		validate = pd.read_csv(validfile)
		cdict = pd.DataFrame.to_dict(validate,orient='list')
		source.data = cdict
	
	if os.path.exists(newcatfile):
		fields = ['x','y','width','height','angles']
		newcat = pd.read_csv(newcatfile, usecols=fields)
		expdict = pd.DataFrame.to_dict(newcat,orient='list')
		source2.data = expdict

def percentil(attr, old, new):
	interval = PercentileInterval(new)
	vmin, vmax = interval.get_limits(sourcemap.data['image'][0])
	color_mapper.low = vmin
	color_mapper.high = vmax
	
def Welliparam(attr, old, new):
	NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
	for Mindex in range(len(source2.data['x'])):
		for key in ['x','y','width','height','angles']:
				NewCoreCat[key].append(source2.data[key][Mindex])
	for index in source2.selected.indices:
		NewCoreCat['width'][index] = ellw.value
		source2.data = NewCoreCat
		
def Helliparam(attr, old, new):
	NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
	for Mindex in range(len(source2.data['x'])):
		for key in ['x','y','width','height','angles']:
				NewCoreCat[key].append(source2.data[key][Mindex])
	for index in source2.selected.indices:
		NewCoreCat['height'][index] = ellh.value
		source2.data = NewCoreCat
		
def Aelliparam(attr, old, new):
	NewCoreCat = dict(x=[],y=[],width=[],height=[],angles=[])
	for Mindex in range(len(source2.data['x'])):
		for key in ['x','y','width','height','angles']:
				NewCoreCat[key].append(source2.data[key][Mindex])
	for index in source2.selected.indices:
		NewCoreCat['angles'][index] = ella.value
		source2.data = NewCoreCat

#Plot image and Table
#####################
	
ph = 700
pw = np.int32(ph *(sourcemap.data['width'][0]/sourcemap.data['height'][0]))

select = Select(title="Region:", value=names[0], options=names)
select.on_change("value", selectmap)
plot = figure(plot_width = pw, plot_height = ph, x_range=(0,sourcemap.data['width'][0]), y_range=(0,sourcemap.data['height'][0]), match_aspect = True,tools=["tap,wheel_zoom,pan"],title=names[0])
plot.axis.visible = False

interval = PercentileInterval(99.5)
vmin, vmax = interval.get_limits(sourcemap.data['image'][0])
color_mapper = LinearColorMapper(palette="Inferno256",low=vmin,high=vmax)

bkg = plot.image(image='image',
           dh='height', dw='width', x=0, y=0, source=sourcemap, color_mapper=color_mapper)

STATUS = ['Undefined', 'VALID','REJECT']
COLOR = ['gray','red','blue']
renderer = plot.ellipse(x="x", y="y", width="width", height="height", angle="angles", line_color=factor_cmap('validation', COLOR, STATUS),fill_alpha=0,source=source,legend_field='validation')

#Trace first position and size for new cores
expert1 = plot.rect('x', 'y', 'width', 'height', source=source2, line_alpha=0, fill_alpha=0)

#Adjust ellipse in rectangle
expert2 = plot.ellipse(x="x", y="y", width="width", height="height", angle="angles", line_color='green',fill_alpha=0,source=source2,legend_label='New cores')

plot.legend.location = "top_right"
plot.legend.click_policy= "hide"

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

#Side Column : Description, Interactions and Buttons
#########################

text = Div(text="""<h1><img src="https://ipag.osug.fr/~robitaij/figures/logo_crop.png" 
            width="149" height="123" align="left" margin=20px>
			LHYRICA Project</h1></br>
This application allows you to validate cores catalogues in order to train LHYRICA's neural network.
""",
width=400)

#LABELS = ["95%", "98%", "99%","99.5%","99.99%"]
#radio_button_group = RadioButtonGroup(labels=LABELS, active=3, max_width = np.int32(pw/2))
#radio_button_group.on_change("active",percentil)
slider = Slider(start=90, end=100, value= 99.5, step=0.01, title="Percentile Interval (image contrast)", max_width = np.int32(pw/2))
slider.on_change("value",percentil)

expert1.data_source.selected.on_change("indices", select_newcore)
ellw = Slider(start=2, end=100, value= 20, step=0.1, title="Ellipse width", max_width = np.int32(pw/2))
ellh = Slider(start=2, end=100, value= 20, step=0.1, title="Ellipse height", max_width = np.int32(pw/2))
ella = Slider(start=-np.pi/2, end=np.pi/2, value= 0, step=0.1, title="Ellipse angle", max_width = np.int32(pw/2))
ellw.on_change("value",Welliparam)
ellh.on_change("value",Helliparam)
ella.on_change("value",Aelliparam)

renderer.data_source.selected.on_change("indices", select_core)

Rbutton = Button(label="Reset view",button_type="default",max_width = np.int32(pw/2))
Rbutton.on_click(reset_view)

Sbutton = Button(label="Save Cat",button_type="success",max_width = np.int32(pw/2))
Sbutton.on_click(save_valid)

Lbutton = Button(label="Load Cat",button_type="primary",max_width = np.int32(pw/2))
Lbutton.on_click(load_valid)

boxedit = BoxEditTool(renderers=[expert1], empty_value = 0.)
expert2.data_source.on_change('data',save_newcat)

center = Button(label="center object",button_type="default",max_width = np.int32(pw/2))
center.on_click(centering)

tooltips = [('index','$index'),('ra', '@ra'), ('dec', '@dec'),('x', '@x'), ('y', '@y')]
plot.add_tools(HoverTool(tooltips=tooltips),boxedit)

#Dashboard Display

tab1 = Panel(child=data_table, title='catalogue')
tab2 = Panel(child=Ndata_table, title='new cores')
tabs = Tabs(tabs=[tab1,tab2])
present = column(text, tabs)
mapint = column(plot,row(column(slider, ellw, ellh, ella),column(center,Rbutton,select,Lbutton,Sbutton)))
panel = row(present, mapint)

curdoc().add_root(panel)