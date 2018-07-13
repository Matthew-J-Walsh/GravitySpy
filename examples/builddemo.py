from pycbc.waveform import get_td_waveform
from pycbc.noise.reproduceable import normal
from gwpy.timeseries import TimeSeries
from gravityspy.plot.plot import plot_qtransform
from gwpy.segments import Segment
from numpy import random
import numpy, pycbc
import os.path
import datetime
from PIL import Image, ImageDraw
from shutil import copyfile

SAMPLE_RATE = 16384
baseTime = 1161907217
plotTimeRanges = [0.5, 1.0, 2.0, 4.0]
plotNormalizedERange = [0, 25.5]
detectorName = 'H1'
searchQRange = [4, 64]
searchFrequencyRange = [10, 2048]
majorfolder = "Simulated_binaries_for_demonstration/"
timeseriesfolder = majorfolder + "Timescale/"
spectrogramsfolder = majorfolder + "Spectrogram/"
visualfolder = majorfolder + "VisualRepresentation/"
trainingfolder = majorfolder + "TrainingSet/"

flow = 30.0
delta_f = 1.0/16
flen = int(8192 / delta_f) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)
delta_t = 1.0/16384
tsamples = int(32 / delta_t)
ts = pycbc.noise.noise_from_psd(tsamples, delta_t, psd, seed=1980)
noise = TimeSeries.from_pycbc(ts)

#noise = 1e-21*TimeSeries(random.normal(scale=.1, size=16384*100), sample_rate=16384)

#draws a visual representation of specified demo
def drawdemo(m1, m2, d, s1x, s1y, s1z, s2x, s2y, s2z, ecc, inc, visualreploc):
	distancebetween = pow((m1+m2)*pow((1.0/31557600)*(1.0/4),2.0),1.0/3)
	img = Image.new('RGB', (500, 500), color = 'white')
	d = ImageDraw.Draw(img)
	d.pieslice([140,240,160,260],0,360,fill='black',outline='black')
	d.pieslice([340,240,360,260],0,360,fill='black',outline='black')
	d.text([150,40], 'distance = '+str(distancebetween)+' AUs',fill='black')
	return img

#takes image locations and places them into a trainingset form in the trainingset folder
def converttotraining(appro, m1, m2, d=1000, s1x=0, s1y=0, s1z=0, s2x=0, s2y=0, s2z=0, ecc=0, inc=0):
	filename = detectorName+'_'+str(m1)+str(m2)+str(d)+"s"+str(s1x)+str(s1y)+str(s1z)+str(s2x)+str(s2y)+str(s2z)+str(ecc)+str(inc)
	traininglocs = [trainingfolder + filename + "_spectrogram_0.5.png",
			trainingfolder + filename + "_spectrogram_1.0.png",
			trainingfolder + filename + "_spectrogram_2.0.png",
			trainingfolder + filename + "_spectrogram_4.0.png",]

	[timeserieslocp, timeserieslocn, superfigloc, subfiglocs, visualreploc] = gendemo(appro, m1, m2, d, s1x, s1y, s1z, s2x, s2y, s2z, ecc, inc, 0)

	for i in range(0,4):
		copyfile(subfiglocs[i], traininglocs[i])
	

#generates demo drawings of certian inputs if they dont already exist and then returns their locations
def gendemo(appro, m1, m2, d=1000, s1x=0, s1y=0, s1z=0, s2x=0, s2y=0, s2z=0, ecc=0, inc=0, overwrite=0):
	demoid = detectorName+'-'+str(m1)+'-'+str(m2)+"-"+str(d)+"-s-"+str(s1x)+"-"+str(s1y)+"-"+str(s1z)+"-"+str(s2x)+"-"+str(s2y)+"-"+str(s2z)+"-"+str(ecc)+"-"+str(inc)
	timeserieslocp = timeseriesfolder + "ts-pure-" + demoid + ".png"
	timeserieslocn = timeseriesfolder + "ts-noisy-" + demoid + ".png"
	superfigloc = spectrogramsfolder + "sps-" + demoid + ".png"
	subfiglocs = [spectrogramsfolder + "sp-" + demoid + "05.png", spectrogramsfolder + "sp-" + demoid + "10.png", spectrogramsfolder + "sp-" + demoid + "20.png", spectrogramsfolder + "sp-" + demoid + "40.png"]
	visualreploc = visualfolder + "vis-" + demoid + ".png"
	logfileloc = majorfolder + "log-" + demoid + ".txt"
	
	#it already exists
	if os.path.isfile(logfileloc) and overwrite==0:
		return [timeserieslocp, timeserieslocn, superfigloc, subfiglocs, visualreploc]

	#make the data
	hp, hc = get_td_waveform(approximant=appro, delta_t=1./SAMPLE_RATE, f_lower=30, mass1=m1, mass2=m2, distance=d, spin1x=s1x, spin1y=s1y, spin1z=s1z, spin2x=s2x, spin2y=s2y, spin2z=s2z, eccentricity=ecc, inclination=inc)
	hp.start_time -= hp.start_time
	signal = TimeSeries.from_pycbc(hp)
	signal.t0 = 50
	plot = signal.plot()
	plot.save(timeserieslocp)
	data = noise.inject(signal)
	data.shift(baseTime)
	plot = data.plot()
	plot.save(timeserieslocn)

	centerTime = baseTime + 50 + hp.duration #- 3277.0 / 16384.0 #this is ~.2 seconds
	startTime = centerTime

	TIMESERIES = data
	specsgrams = []

	for iTimeWindow in plotTimeRanges:
		#print("w")
		durForPlot = iTimeWindow/2
		try:
			outseg = Segment(centerTime - durForPlot, centerTime + durForPlot)
			qScan = TIMESERIES.q_transform(qrange=tuple(searchQRange), frange=tuple(searchFrequencyRange), gps=centerTime, search=0.5, tres=0.002, fres=0.5, outseg=outseg, whiten=True)
			qValue = qScan.q
			qScan = qScan.crop(centerTime-iTimeWindow/2, centerTime+iTimeWindow/2)
		except:
			#print("failed")
			outseg = Segment(centerTime - 2*durForPlot, centerTime + 2*durForPlot)
			qScan = TIMESERIES.q_transform(qrange=tuple(searchQRange), frange=tuple(searchFrequencyRange),gps=centerTime, search=0.5, tres=0.002,fres=0.5, outseg=outseg, whiten=True)
			qValue = qScan.q
			qScan = qScan.crop(centerTime-iTimeWindow/2, centerTime+iTimeWindow/2)
		specsgrams.append(qScan)

	indFigAll, superFig = plot_qtransform(specsgrams, plotNormalizedERange, plotTimeRanges, detectorName, startTime)

	superFig.save(superfigloc)
	indFigAll[0].save(subfiglocs[0])
	indFigAll[1].save(subfiglocs[1])
	indFigAll[2].save(subfiglocs[2])
	indFigAll[3].save(subfiglocs[3])

	#draw the aproximation
	img = drawdemo(m1, m2, d, s1x, s1y, s1z, s2x, s2y, s2z, ecc, inc, visualreploc)
	img.save(visualreploc)

	#make the log file
	f = open(logfileloc, 'w')
	f.write(demoid+"\n")
	f.write(str(datetime.datetime.now()))
	f.close()

	return [timeserieslocp, timeserieslocn, superfigloc, subfiglocs, visualreploc]










