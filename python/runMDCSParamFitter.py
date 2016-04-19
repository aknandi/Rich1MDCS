import GaudiPython
import Gaudi
import datetime
from GaudiPython import gbl
from datetime import *
from Gaudi.Configuration import *
from Configurables import LHCbApp, MDMFitParameters, EventClockSvc,DDDBConf, CondDB, CondDBAccessSvc
from Configurables import Brunel
import os,sys
from optparse import OptionParser, OptionGroup

def getUNIXTime(dtime):
	# Note dtime must be a date in CET (CERN) time.
	import time
	t = time.mktime(dtime.timetuple())
	zone = time.tzname[0]
	if zone not in ['GMT','CET'] : raise Exception('Unknown time zone '+zone)
	offset = 0
	if zone == 'GMT' : offset = -3600
	return int( (t+offset) * 1e9 )

def umsSvc():
	import GaudiPython
	GaudiPython.AppMgr().createSvc('UpdateManagerSvc')
	return GaudiPython.AppMgr().service('UpdateManagerSvc','IUpdateManagerSvc')

def iDetDataSvc():
	import GaudiPython
	return GaudiPython.AppMgr().service('DetectorDataSvc','IDetDataSvc')

if __name__ == '__main__' :
	#play with the magnet
	#from Configurables import UpdateManagerSvc
	#UpdateManagerSvc().ConditionsOverride += ["Conditions/Online/LHCb/Magnet/Measured := double Current = 0.1 ; int Polarity = 1;"]
	parser = OptionParser()
	parser.set_usage("python runMDCSParamFitter.py -R RUNNUMBER [extra options]");
	parser.add_option("-R", "--RunNumber", dest="runNumber", help="process RUNNUMBER", metavar="RUNNUMBER", type="string")
	parser.add_option("-M", "--MagPolarity", dest="polarity", help="Scan POLARITY", metavar="POLARITY", type="string")
	parser.add_option("-P", "--Rich1Panel", dest="panel", help="RICH1 PANEL", metavar="PANEL", type="string")
	parser.add_option("-I", "--InputFile", dest="inputfile", help="Input file (Cleaned ntuple)", metavar="FILE", type="string")
	#Additional -S systematics option (added 20/2/14)
	#Takes an integer from shell script loop, and adds it to the end of the file name
	#Saves each fit to a separate nTuple which we can later hadd
	parser.add_option("-S", "--Systematics", dest="systematics", help="Systematic Fit Number (provide with a loop in sheel script)", metavar="SYST", type="string")
	#New random seed for systematics jobs
	parser.add_option("-T", "--RandomSeed", dest="randomSeed", help="Random seed for systematics", metavar="RANDOMSEED", type="string")
	#New radius change factor for systematics jobs
	parser.add_option("-F", "--RadiusFactor", dest="radiusFactor", help="Radius change factor for systematics", metavar="RADIUSCHANGE", type="string")
	group =OptionGroup(parser, "Other Options","These options may not be always necessary")
	group.add_option("--DDDBtag", dest="DDDBtag", help="Set custom DDDB tag",type="string",metavar="TAG")
	group.add_option("--CondDBtag", dest="CondDBtag", help="Set custom CondDB tag",type="string",metavar="TAG")
	group.add_option("-O", "--OutputXML", dest="outputxml", help="set output xml file name",metavar="OUTPUTNTUPLE")
	group.add_option("-N", "--TupleName", dest="tuplename", help="set name of the output ntuple",metavar="FILENAME")

	parser.add_option_group(group)

	(options, args) = parser.parse_args()

	dict_opts=options.__dict__
	
	if dict_opts["runNumber"]==None:
		print "Please provide a run number!"
		print "run: \"python runMDCSPeakFinder.py -h\" for help"
		sys.exit()
	if dict_opts["polarity"]==None:
		print "Please specify magnet polarity!!"
		print "run: \"python runMDCSPeakFinder.py -h\" for help"
		sys.exit()
	if dict_opts["panel"]==None:
		print "Please specify a panel!!"
		print "run: \"python runMDCSPeakFinder.py -h\" for help"
		sys.exit()
	blacklist=[]
	for opts, val in dict_opts.iteritems():
	        if dict_opts[opts]==None:
			blacklist.append(opts)
	for bkeys in blacklist:
		del dict_opts[bkeys]	
		
#	print dict_opts
	outuple="MDCSOutput_Run"+dict_opts["runNumber"]+"_Magnet"+dict_opts["polarity"]+"_"+dict_opts["panel"]+"_RadiusChange_"+dict_opts["radiusFactor"]+"percent_"+dict_opts["systematics"]+".root"
	if dict_opts.has_key("tuplename"):
		outuple=dict_opts["tuplename"]
	outxml="MDCSParameters_Run"+dict_opts["runNumber"]+"_Magnet"+dict_opts["polarity"]+"_"+dict_opts["panel"]+"_RadiusChange_"+dict_opts["radiusFactor"]+"percent_"+dict_opts["systematics"]+".xml"
	if dict_opts.has_key("outputxml"):
		outxml=dict_opts["outputxml"]

	
	
	MDCSFitter = MDMFitParameters("MDCSFitParameters_run"+dict_opts["runNumber"]+"_Magnet"+dict_opts["polarity"]+"_"+dict_opts["panel"])
	MDCSFitter.SetRunNumber=int(dict_opts["runNumber"])
	MDCSFitter.SetCleanedNTuple =dict_opts["inputfile"]
	MDCSFitter.OutputTuple =outuple
	MDCSFitter.Output_MDMS_DB_Xml=outxml
	###Set the random seed - new MARCH 2014 (means both U and D box get the same seed for a given fit)
	MDCSFitter.SetRandomSeed = int(dict_opts["randomSeed"])
	###Set the radius change factor for the systematics fits (in percent, later divided by 100 in the code)
	MDCSFitter.SetRadiusFactor = float(dict_opts["radiusFactor"])
	
	mg=0
	if dict_opts["polarity"]=="Up":
		mg=1
	if dict_opts["polarity"]=="Down":
		mg=-1
	
	pn=1
	if dict_opts["panel"]=="U":
		pn=0
	print int(mg)		
	MDCSFitter.BField=int(mg)
	MDCSFitter.Panel=int(pn)
	
#	LHCbApp().DataType="2011"
	LHCbApp().Simulation=False
#	LHCbApp().DDDBtag="HEAD"
#	LHCbApp().CondDBtag="HEAD"

	if dict_opts.has_key("DDDBtag"):
		LHCbApp().DDDBtag= dict_opts["DDDBtag"]
	if dict_opts.has_key("CondDBtag"):
		LHCbApp().CondDBtag= dict_opts["CondDBtag"]
		#CondDB().addLayer(
			#CondDBAccessSvc("HPDAlign",
				#	ConnectionString="sqlite_file:2012-RootFiles-RunAligned-Sobel-AveragePol0-HPDAlign-13092012 2.db/LHCBCOND",
					#DefaultTAG="HEAD") )
		
	ApplicationMgr().TopAlg += [MDCSFitter]
	ApplicationMgr().ExtSvc += ['DataOnDemandSvc']
	ApplicationMgr().EvtMax = 1
#	ApplicationMgr().OutputLevel=DEBUG
	EventSelector().PrintFreq = 10

	appMgr = GaudiPython.AppMgr()
	appMgr.HistogramPersistency="NONE"
#Date of magnet up scan
        #mydate_UP=datetime(2012,4,27,12,0,0)
#Date of magnet down scan
        #mydate_UP=datetime(2012,3,27,12,0,0)
	#Old magnet DOWN scan for reference
	#mydate_UP=datetime(2010,8,30,12,0,0)
	#magnet DOWN scan November 2012
	#mydate_UP=datetime(2012,11,29,12,0,0)
	#July date 2012 for testing sensitivity to choice of DB tags
	#mydate_UP=datetime(2012,7,20,12,0,0)

	#Most recent scans taken a day apart in January
	mydate_UP=datetime(2013,1,30,12,0,0)
	
	
	unixTime_UP=getUNIXTime(mydate_UP)
	iDetDataSvc().setEventTime( gbl.Gaudi.Time(unixTime_UP))
	umsSvc().newEvent()
	
	appMgr.initialize()

	

#	appMgr.run(1)
#	appMgr.stop()
	appMgr.finalize()
