import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from scipy.fftpack import fft,fftfreq#forFourierTransform

"""
Parameters:
	saveImages = True - it will save images of every plot in the folder outImages_*
	plotTimeEveryAxis = True - it will put the time on every axis, otherwise only as a title on the main figure

"""

saveImages = False
#saveImages = True

plotTimeEveryAxis = True

useHlines = True


def relimAxis(ax, title, setLimits = False):

	"""
		it will relimit the specified axis by title in concordance with getXLimits and getYLimits functions from domain_plot_limits.py
	"""
	from domain_plot_limits import  getYLimits ,getXLimits
	ylim = getYLimits(title)
	xlim = getXLimits(title)
	if(xlim is None or ylim is None):
		ax.relim()
	if xlim is None and ylim is None:
		ax.autoscale_view(True,True,True)
	else:
		#I only have to set this limits once at the beginning
		if setLimits:	
			if ylim:
				from matplotlib.ticker import FormatStrFormatter
				ax.set_ylim(ylim["minY"],ylim["maxY"])
				ax.yaxis.set_major_formatter(FormatStrFormatter('%f'))
			if xlim:
				ax.set_xlim(xlim["minX"],xlim["maxX"])
		#afterwards the limits are already set and I only have to autoscale 
		if not xlim:
			ax.autoscale_view(True,True,False)
		else:
			ax.autoscale_view(True,False,True)


#ylim = {"pres":{ "maxY": 1.0005, "minY": 0.9995} , "vel" : { "maxY": 0.00035, "minY": -0.00035}, "rho":{ "maxY": 1.0004, "minY": 0.9996}, 'rhoCurve': { "maxY": 0.00025, "minY": -0.00025} , 'velFFT' : {"maxX" : 80, "minX" : -80} } 
#inhom1
#ylim = {"pres":{ "maxY": 1.0006, "minY": 0.9995} , "vel" : { "maxY": 0.0015, "minY": -0.0015}, "rho":{ "maxY": 1.0004, "minY": 0}, 'rhoCurve': { "maxY": 0.00035, "minY": -0.00035} , 'velFFT' : {"maxX" : 80, "minX" : -80} } 
#inhom2
#ylim = {"pres":{ "maxY": 1.0008, "minY": 0.9992} , "vel" : { "maxY": 0.0015, "minY": -0.0015}, "rho":{ "maxY": 1.3, "minY": 0}, 'rhoCurve': { "maxY": 0.0020, "minY": -0.0020}, 'velFFT' : {"maxX" : 80, "minX" : -80} } 
#ylim = None


def getRandomColor():
	from random import random
	red = random()
	blue = random()
	green = random()	
	return "#%02x%02x%02x" % (red*255, green*255, blue*255)

def getColorFromArray(array):
	return "#%02x%02x%02x" % (array[0][0]*255, array[0][1]*255, array[0][2]*255)



	


class VisualPlot:

	"""
	class responsible with ploting
	"""


	def __init__(self, z, titles, iniValues):
		if saveImages:
			from common import createFolder	
			self.dirname = createFolder("outImages")
		nplots = len(titles)
		fig, ax =  plt.subplots(nplots,1,True)
		#fig.set_size_inches(300,200)
		#fig.set_figwidth(300)
		#fig.set_figheight(200)
		if(nplots == 1):
			ax = [ax]
		self.figures = [fig]
		self.lines = {}
		self.axes = {}
		for i in range(0, nplots):
			title = titles[i]
			self.addAxis(z, ax[i], title, iniValues[i])
		#now I piut time in every figure
		if not plotTimeEveryAxis:
			self.plotTitle = ax[0].set_title("Time 0")
		wm = plt.get_current_fig_manager()
		wm.window.wm_geometry("1000x900+50+50")
		fig.subplots_adjust(right=0.8)
		plt.draw()
		plt.show(block=False)

	def afterInit(self):
		#import time
		#time.sleep(5)
		#save initial figures to files
		if saveImages:
			numFig = 0
			for fig in self.figures:
				os.mkdir(os.path.join(self.dirname, "Fig%d" % numFig))
				fig.savefig(os.path.join(self.dirname, "Fig%d" % numFig , "img000000.png"))
				numFig +=1


	def addAxis(self, z, ax, title, vals, xlabel="z", linestyle="-"):
		ax.set_xlabel(xlabel)
		ax.set_ylabel(title)
		ax.grid(True)
		self.axes[title] = ax
		shape = np.shape(vals)
		plotLegend = False
		#we can plot multiple graphs on the same axis : example numerical and analytical
		if(len(shape)==1):
			if useHlines:
				ax.hlines(np.min(vals), z[0], z[len(z)-1])
				ax.hlines(np.max(vals), z[0], z[len(z)-1])
			if(linestyle == '-'):
				l, = ax.plot(z, vals, lw=2, color='b')
			else:
				l, = ax.plot(z, vals,  markersize=3, linestyle="None", marker="o",  color='b')
			#l, = ax.plot(z, vals[i], lw=2, color='b', markersize=5, linestyle="-", marker="o")
			self.lines[title] = l
		elif(len(shape)==2):
			self.lines[title] = []
			if useHlines:
				ax.hlines(np.min(vals[0]), z[0], z[len(z)-1])
				ax.hlines(np.max(vals[0]), z[0], z[len(z)-1])
			for i in range(0, shape[0]):
				if(linestyle == '-'):
					if(len(z.shape) == 2):
						l, = ax.plot(z[i], vals[i], lw=2, color=getRandomColor(),  label="%d" % i)
					else:
						l, = ax.plot(z, vals[i], lw=2, color=getRandomColor(),  label="%d" % i)
				else:
					if(len(z.shape) == 2):
						l, = ax.plot(z[i], vals[i], color=getRandomColor(),  markersize=3, linestyle="None", marker="o",  label="%d" % i)
					else:
						l, = ax.plot(z, vals[i], color=getRandomColor(),  markersize=3, linestyle="None", marker="o",  label="%d" % i)
				plotLegend = True
				#l, = ax.plot(z, vals[i], lw=2, color=getRandomColor(), markersize=5, linestyle="-", marker="o", label="%d" % i)
				self.lines[title].append(l)
#		ylim = getYLimits(title)
#		xlim = getXLimits(title)
#		if not xlim and not ylim:
#			ax.relim()
#			ax.autoscale_view(True,True,True)
#		else:
#			if (ylim):
#				from matplotlib.ticker import FormatStrFormatter
#				ax.set_ylim(ylim["minY"],ylim["maxY"])
#				ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
#			if xlim:
#				ax.set_xlim(xlim["minX"],xlim["maxX"])
		if(plotLegend):
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
		relimAxis(ax, title, True)	
		#ax.relim()
		#ax.autoscale_view(True,True,True)

	def markPoint(self, axTitle, pointName, value):
		#print("markPoint on axTitle = %s, pointName = %s, value = %E" % (axTitle, pointName, value))
		from common import testKeyInDict
		if not hasattr(self, 'markPoints'):
			self.markPoints = {}
		if(testKeyInDict(pointName, self.markPoints)):
			self.markPoints[pointName].remove()
			#keep color
			color = getColorFromArray(self.markPoints[pointName].get_color())
			del  self.markPoints[pointName]
		else:
			#generate random color
			color = getRandomColor()
		l = self.lines[axTitle]
		if hasattr(l, '__len__'):
			l = l[0]
		yvals = l.get_ydata()
		minValue =  np.min(yvals)
		maxValue =  np.max(yvals)
		delta = 0.2
		#I have to make the following test because
		#sometimes (in the case of riemann problem and initial velocity  0 )
		#because zC is only marked once at the beginning  when velocity is 0 for all z
		#TODO I choose 1 but it might be too small
		if(maxValue == minValue):
			maxValue = minValue + 1
		self.markPoints[pointName] = self.axes[axTitle].vlines(value, minValue - delta, maxValue + delta, color=color, label=pointName)
		self.axes[axTitle].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
			

	def addGraph(self, z, title, vals, xlabel="z", linestyle="-"):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		self.addAxis(z, ax, title, vals, xlabel, linestyle)
		self.figures.append(fig)
		fig.subplots_adjust(right=0.8)
		plt.draw()
		plt.show(block=False)

	def plotAxisTwin(self, z, title, vals, newtitle):
		ax2 = self.axes[title].twinx()
		ax2.set_ylabel(newtitle)
		ax2.plot(z, vals, color=getRandomColor())
		plt.draw()
		plt.show(block=False)

	def plotAxis(self, z, title, vals, label=None):
		ax = self.axes[title]
		ax.plot(z, vals, color=getRandomColor(),label=label)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
		plt.draw()
		plt.show(block=False)

	def updateValues(self, title, newValues, newTime, newZ = None, sortValues = False):
		#TODO newZ as tuple
		if not newZ is None:
			sortValues = True
		#print("updateValues %s" % title)
		shape = np.shape(newValues)
		#we can plot multiple graphs on the same axis : example numerical and analytical: see addAxis before!!

		def getXYSort(x,y,sortValues):
			if sortValues:
				argSort = np.argsort(x)
				resx = x[argSort]	
				resy = y[argSort]
			else:
				resx = x
				resy = y
			return (resx,resy)
		

		if(len(shape)==1):
			if not newZ is None:
				#I know newz will be in form [(zararray,sortBool)] no test
				newZP, newValuesP = getXYSort(newZ[0], newValues, newZ[1])
				self.lines[title].set_xdata(newZP)
			else:
				newValuesP = newValues


			self.lines[title].set_ydata(newValuesP)

			#print(" ".join(map(str, newValues)))
		elif(len(shape)==2):
			#print(" ".join(map(str, newValues[0])))
			nlines = shape[0]
			if(hasattr(self.lines[title], "__len__") and len(self.lines[title])==nlines):
				for i in range(0, nlines):
					if not newZ is None:
					#I know newz will be in form [(zararray,sortBool)] with length equal to y length no test
						newZP, newValuesP = getXYSort(newZ[i][0], newValues[i], newZ[i][1])
						self.lines[title][i].set_xdata(newZP)
					else:
						newValuesP = newValues[i]
					self.lines[title][i].set_ydata(newValuesP)
		if plotTimeEveryAxis:	
			self.axes[title].set_title("Time %4.4f" % newTime)
		relimAxis(self.axes[title], title)	
		
		
	def afterUpdateValues(self, newTime):
		if not plotTimeEveryAxis:
			self.plotTitle.set_text("Time %4.4f" % newTime)
		numFig = 0
		for fig in self.figures:
			#this will print over old one:
			# I should save the var
			#http://stackoverflow.com/questions/10559144/matplotlib-suptitle-prints-over-old-title
			#fig.suptitle("Time %4.4f" % newTime)
			fig.canvas.draw()
			if saveImages:
				#make name compatible with ffmpeg
				#ffmpeg -r 1 -i img%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
				#the above does not work
				#ffmpeg -r 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4
				#convert HANGS!!
				#convert -antialias -delay 1x2 *.png mymovie.mp4
				imgname = "%4.4f" % newTime
				if(len(imgname) == 6):
					imgname = "0"+imgname
				imgname = imgname.replace(".", "")	
				fig.savefig(os.path.join(self.dirname, "Fig%d" % numFig , "img%s.png"%imgname))
			numFig +=1
		#import time
		#time.sleep(5)

	def finish(self):
		import time
		time.sleep(10)
		#pass







