from pycbc.waveform import get_td_waveform
from pycbc.noise.reproduceable import normal
from gwpy.timeseries import TimeSeries
from gravityspy.plot.plot import plot_qtransform
from gwpy.segments import Segment
from numpy import random
import gc

import numpy, pycbc.psd

SAMPLE_RATE = 16384


lmrange = [21, 25]
hmrange = [130, 130]
baseTime = 1161907217

plotTimeRanges = [0.5, 1.0, 2.0, 4.0]
plotNormalizedERange = [0, 25.5]
detectorName = 'L1'
searchQRange = [4, 64]
searchFrequencyRange = [10, 2048]

noise = 1e-21*TimeSeries(random.normal(scale=.1, size=16384*100), sample_rate=16384)

for i in range(lmrange[0], lmrange[1], 1):
	print("lm "+str(i))
	hp, hc = get_td_waveform(approximant="TaylorT4", delta_t=1./SAMPLE_RATE, f_lower=30, mass1=i,mass2=3, distance=1000)
	hp.start_time -= hp.start_time
	signal = TimeSeries.from_pycbc(hp)
	signal.t0 = 50
	data = noise.inject(signal)
	data.shift(baseTime)
	plot = data.plot()
	plot.save("fakes/ts- lowmass: 1 - "+str(i)+".png")

	print("t-plot saved. duration = "+str(hp.duration))
	centerTime = baseTime + 50 + hp.duration
	startTime = centerTime
	print(str(data.t0)+" "+str(startTime))

	TIMESERIES = data
	specsgrams = []

	for iTimeWindow in plotTimeRanges:
		print("w")
		durForPlot = iTimeWindow/2
		try:
        		outseg = Segment(centerTime - durForPlot, centerTime + durForPlot)
        		qScan = TIMESERIES.q_transform(qrange=tuple(searchQRange), frange=tuple(searchFrequencyRange), gps=centerTime, search=0.5, tres=0.002, fres=0.5, outseg=outseg, whiten=True)
        		qValue = qScan.q
        		qScan = qScan.crop(centerTime-iTimeWindow/2, centerTime+iTimeWindow/2)
    		except:
			print("failed")
        		outseg = Segment(centerTime - 2*durForPlot, centerTime + 2*durForPlot)
        		qScan = TIMESERIES.q_transform(qrange=tuple(searchQRange), frange=tuple(searchFrequencyRange),gps=centerTime, search=0.5, tres=0.002,fres=0.5, outseg=outseg, whiten=True)
        		qValue = qScan.q
			qScan = qScan.crop(centerTime-iTimeWindow/2, centerTime+iTimeWindow/2)
		specsgrams.append(qScan)

	print("compiling figure")
	indFigAll, superFig = plot_qtransform(specsgrams, plotNormalizedERange, plotTimeRanges, detectorName, startTime)

	superFig.save("fakes/spec- lowmass 1 -"+str(i)+".png")
	indFigAll[0].save("fakes/specs1- lowmass 1 -"+str(i)+".png")
	indFigAll[1].save("fakes/specs2- lowmass 1 -"+str(i)+".png")
	indFigAll[2].save("fakes/specs3- lowmass 1 -"+str(i)+".png")
	indFigAll[3].save("fakes/specs4- lowmass 1 -"+str(i)+".png")

	gc.collect()



print("lowmass compiled")

for i in range(hmrange[0], hmrange[1], 5):
	print("hm -"+str(i))
	hp, hc = get_td_waveform(approximant="TaylorT4", delta_t=1./SAMPLE_RATE, f_lower=30, mass1=i,mass2=35, distance=1000)
	hp.start_time -= hp.start_time
	signal = TimeSeries.from_pycbc(hp)
	signal.t0 = 50
	data = noise.inject(signal)
	data.shift(baseTime)
	plot = data.plot()
	plot.save("fakes/ts- highmass: 35 - "+str(i)+".png")

	print("t-plot saved")
	centerTime = baseTime + 50 + hp.duration
	startTime = centerTime

	TIMESERIES = data
	specsgrams = []

	for iTimeWindow in plotTimeRanges:
		durForPlot = iTimeWindow/2
	        outseg = Segment(centerTime - durForPlot, centerTime + durForPlot)
		qScan = TIMESERIES.q_transform(qrange=tuple(searchQRange), frange=tuple(searchFrequencyRange), gps=centerTime, search=0.5, tres=0.002, fres=0.5, outseg=outseg, whiten=True)
		qValue = qScan.q
		qScan = qScan.crop(centerTime-iTimeWindow/2, centerTime+iTimeWindow/2)
		specsgrams.append(qScan)

	indFigAll, superFig = plot_qtransform(specsgrams, plotNormalizedERange, plotTimeRanges, detectorName, startTime)

	superFig.save("spec- highmass 35 -"+str(i)+".png")
	indFigAll[0].save("specs1- highmass 35 -"+str(i)+".png")
	indFigAll[1].save("specs2- highmass 35 -"+str(i)+".png")
	indFigAll[2].save("specs3- highmass 35 -"+str(i)+".png")
	indFigAll[3].save("specs4- highmass 35 -"+str(i)+".png")


print("highmass compiled")

#signal = TimeSeries.from_pycbc(hp)
##signal.taper()

#signal.t0 = 10
#data = noise.inject(signal)
#data.shift(baseTime)

#plot = data.plot()
#plot.save("ztimeseries.png")

#code from test_plot.py
#TIMESERIES = data

#i = 1
"""for sig in lowmass:
	TIMESERIES = sig
	specsgrams = []

	for iTimeWindow in plotTimeRanges:
		durForPlot = iTimeWindow/2
	        outseg = Segment(centerTime - durForPlot, centerTime + durForPlot)
		qScan = TIMESERIES.q_transform(qrange=tuple(searchQRange), frange=tuple(searchFrequencyRange), gps=centerTime, search=0.5, tres=0.002, fres=0.5, outseg=outseg, whiten=True)
		qValue = qScan.q
		qScan = qScan.crop(centerTime-iTimeWindow/2, centerTime+iTimeWindow/2)
		specsgrams.append(qScan)

	indFigAll, superFig = plot_qtransform(specsgrams, plotNormalizedERange, plotTimeRanges, detectorName, startTime)

	superFig.save("spec- lowmass 1 -"+str(i)+".png")
	indFigAll[0].save("specs1- lowmass 1 -"+str(i)+".png")
	indFigAll[1].save("specs2- lowmass 1 -"+str(i)+".png")
	indFigAll[2].save("specs3- lowmass 1 -"+str(i)+".png")
	indFigAll[3].save("specs4- lowmass 1 -"+str(i)+".png")
	i += 1

print("lowmass made")

i = 35
for sig in highmass:
	TIMESERIES = sig
	specsgrams = []

	for iTimeWindow in plotTimeRanges:
		durForPlot = iTimeWindow/2
	        outseg = Segment(centerTime - durForPlot, centerTime + durForPlot)
		qScan = TIMESERIES.q_transform(qrange=tuple(searchQRange), frange=tuple(searchFrequencyRange), gps=centerTime, search=0.5, tres=0.002, fres=0.5, outseg=outseg, whiten=True)
		qValue = qScan.q
		qScan = qScan.crop(centerTime-iTimeWindow/2, centerTime+iTimeWindow/2)
		specsgrams.append(qScan)

	indFigAll, superFig = plot_qtransform(specsgrams, plotNormalizedERange, plotTimeRanges, detectorName, startTime)

	superFig.save("spec- highmass 35 -"+str(i)+".png")
	indFigAll[0].save("specs1- highmass 35 -"+str(i)+".png")
	indFigAll[1].save("specs2- highmass 35 -"+str(i)+".png")
	indFigAll[2].save("specs3- highmass 35 -"+str(i)+".png")
	indFigAll[3].save("specs4- highmass 35 -"+str(i)+".png")
	i += 5

"""






