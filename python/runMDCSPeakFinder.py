from optparse import OptionParser, OptionGroup
import os, sys, getpass


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
	
	parser = OptionParser()
	parser.set_usage("python runMDCSPeakFinder.py -R RUNNUMBER [extra options]");
	parser.add_option("-R", "--RunNumber", dest="runNumber", help="process RUNNUMBER", metavar="RUNNUMBER", type="int")
	
	group =OptionGroup(parser, "Other Options","These options may not be always necessary")
	group.add_option("-D", "--DEBUG", action="store_true", dest="DEBUG", help="print debug information")
	#group.add_option("-G","--runWithGanga", action="store_true", dest="ganga", help="Submit the job on the batch system through Ganga")
	group.add_option("-H", "--StoreHistos", action="store_true", dest="StoreHistos", help="store root files with hitmap histograms information")
	group.add_option("--DDDBtag", dest="DDDBtag", help="Set custom DDDB tag",type="string",metavar="TAG")
	group.add_option("--CondDBtag", dest="CondDBtag", help="Set custom CondDB tag",type="string",metavar="TAG")
	group.add_option("--addCondDBLayer", dest="addCondDBLayer", help="Add custom CondDB layer",type="string",metavar="LAYER")
	group.add_option("--addCondDBLayer_tag", dest="addCondDBLayer_tag", help="Set custom CondDB layer tag, if different from \"HEAD\"",type="string",metavar="TAG")
	group.add_option("-O", "--OutputDirectory", dest="outputdir", default="", help="set output directory directory (you'll need some space..) [default: %default]",metavar="DIRECTORYPATH")
	group.add_option("-N", "--TupleName", dest="TupleName", help="set name of the output ntuple",metavar="FILENAME")
	group.add_option("--HistoName", dest="HistoName", help="set name of the output histogram file",metavar="FILENAME")
	group.add_option("-T", "--TotSteps", dest="totsteps", help="total number of steps in the scan [default: %default]", metavar="NSTEPS", type="int",default=2200)
	group.add_option("--EvtPerStep", dest="nEvtsPerStep", help="number of events in each step [default: %default]", metavar="NEVENTS", type="int",default=6000)
	group.add_option("--IgnoreHeartBeat", dest="IgnoreHeartBeat", help="Set IgnoreHeartBeat=True of events in CondDB", action="store_true")
	
	group_one =OptionGroup(parser, "1 Job Options","You must provide these options in case you do not want to process the entire scan but just one subjob")
	group_one.add_option("-S", "--StartStep", dest="start", help="first calibration step to process [default: %default]", metavar="STEP", type="int")
	group_one.add_option("-E", "--EndStep", dest="end", help="last calibration step to process [default: %default]", metavar="STEP", type="int")
	
	parser.add_option_group(group)
	parser.add_option_group(group_one)
	
	(options, args) = parser.parse_args()
	
	dict_opts=options.__dict__
	
	if dict_opts["runNumber"]==None:
		print "Please provide a run number!"
		print "run: \"python runMDCSPeakFinder.py -h\" for help"
		sys.exit()

	
	#clean dictionary from non-specified options
	blacklist=[]
	for opts, val in dict_opts.iteritems():
		if dict_opts[opts]==None:
			blacklist.append(opts)
	for bkeys in blacklist:
		del dict_opts[bkeys]
	
	if dict_opts.has_key("end") and dict_opts.has_key("start"):
		print "Running only one subjob"
	else:
		print "Running all steps"
		dict_opts["start"]=0
		dict_opts["end"]=dict_opts["totsteps"]
		if not dict_opts.has_key('TupleName'):
			print "setting tuple name"
			dict_opts["TupleName"]="NTuples_Run%i_AllSteps.root"%(dict_opts["runNumber"])
		if not dict_opts.has_key('HistoName'):
			print "setting histogram name"
			dict_opts["HistoName"]="Histos_Run%i_AllSteps.root"%(dict_opts["runNumber"])
	#print dict_opts
	from rawDataToNtuple import *
	rawDataToNtuple(dict_opts)


		#Set system time
	appMgr = GaudiPython.AppMgr()
	appMgr.HistogramPersistency="NONE"
	appMgr.initialize()
#Date of magnet up scan
        #mydate_UP=datetime(2012,4,27,12,0,0)
#Date of magnet down scan
        #mydate_UP=datetime(2012,3,27,12,0,0)
	#Old magnet DOWN scan for reference
	#mydate_UP=datetime(2010,8,30,12,0,0)
	#magnet DOWN scan November 2012
	#mydate_UP=datetime(2012,11,29,12,0,0)
	#July date 2012 for testing sensitivity to choice of DB -tags
	#mydate_UP=datetime(2013,1,30,12,0,0)
	#Date for 2010 reprocessing (in March 2014)
	#mydate_UP= datetime.datetime(2010,4,10,12,0,0)

	#mydate_UP= datetime.datetime(2016,3,16,14,46,0)
	
	#unixTime_UP=getUNIXTime(mydate_UP)
	#iDetDataSvc().setEventTime( gbl.Gaudi.Time(unixTime_UP))
	#umsSvc().newEvent()
	
	#appMgr.initialize()

	
#	appMgr.execute()
#	appMgr.run(1)
#	appMgr.stop()
	appMgr.finalize()
