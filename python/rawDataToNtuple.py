#!/usr/bin/env python
import os
import sys
import glob
from array import array
import socket
import getpass
import datetime

#if(os.environ['ENVIRONMENT']=='BATCH'):
#	sys.argv.append('-b')

import GaudiPython
import Gaudi
from Gaudi.Configuration import *


from ROOT import TSystem,TStyle,TROOT,TColor,FileStat_t,gStyle,gROOT,gSystem

#gSystem = TSystem.gSystem
#gStyle = TROOT.gStyle
#gROOT = TROOT.gROOT
gROOT.SetStyle("Plain")

nColors=51
s = [0.00, 0.25, 0.50, 0.75, 1.00]
r = [0.99, 0.00, 0.87, 1.00, 0.70]
g = [0.99, 0.81, 1.00, 0.20, 0.00]
b = [0.99, 1.00, 0.12, 0.00, 0.00]
ss = array('d',s)
rr = array('d',r)
gg = array('d',g)
bb = array('d',b)
TColor.CreateGradientColorTable(len(ss), ss, rr, gg, bb, nColors)
gStyle.SetNumberContours(nColors)

def exists(file):
	x=FileStat_t()
	return (not gSystem.GetPathInfo("%s"%(file),x))

def getData(runNumber,fstep,lstep,tsteps,islocal=False):
	print "Retrieving data..."
	if islocal:
		print "Searching in local cluster"
	else:
		print "Searching on CASTOR"
	
	preambles={}
	now = datetime.datetime.now()
	
	# First look in the users local eos area RichMDCS/{year}/
	if islocal:
		preamble_string="root://eoslhcb.cern.ch//eos/lhcb/user/%s/%s/RichMDCS/"%(getpass.getuser()[0],getpass.getuser())
		post_string = ""
	# Look for the files on castor if the data cannot be found locally
	else:
		preamble_string="root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/"
		post_string="/RAW/FULL/RICH1/TEST"

	# Construct full file path for a range of years
	for yr in range(2015,now.year+1):
		preambles[str(yr)]=preamble_string+str(yr)+post_string
			
	LIST=[]
	the_year=2009
	# Look for the raw files in each of the years
	for yr,preamble in preambles.iteritems():
		file1 = "%s/%d/%06d_%010d.raw"%(preamble,runNumber,runNumber,1)
                print file1
		if exists(file1):
			the_year=int(yr)
			for nfl in range(1,20):
				file = "%s/%d/%06d_%010d.raw"%(preamble,runNumber,runNumber,nfl)
				fileURL = "mdf:%s/%d/%06d_%010d.raw"%(preamble,runNumber,runNumber,nfl)
				if exists(file):
					LIST.append(fileURL)
			break

	if(len(LIST)==0):
		print "No data "+str(runNumber)+" files found at any of "
		print preambles
		#sys.exit()
	
	# BIT OF JIGGING TO TRY TO GET JOB TO START AT AN "INTELLIGENT" PLACE IN THE LIST OF FILES
	nfiles=len(LIST)
	fs=float(fstep)
	ls=float(lstep)
	tot=float(tsteps)
	
	frac_s=fs/tot
	frac_e=ls/tot
	
	firstfile=int(nfiles*frac_s)-10
	lastfile=int(nfiles*frac_e)+1
	
	if lastfile>(nfiles-1):
		lastfile=nfiles-1
	
	if firstfile<0:
		firstfile=0
	
	print "N files: "+str(nfiles)
	print "First file: "+str(firstfile)
	print "Last file: "+str(lastfile)
	
	DATA=[]
	for i in range(firstfile,lastfile+1):
		DATA.append("DATAFILE=\'%s\' SVC=\'LHCb::MDFSelector\'"%LIST[i])
		print DATA[i]
	return {"DATA":DATA,"year":the_year}

def rawDataToNtuple(options):
#	print options
	required_options=["runNumber","start","end","outputdir","nEvtsPerStep","totsteps"]
	
	for check_opts in required_options:
		if not options.has_key(check_opts):
			print "Please specify minimal options!"
			print "Option \'"+check_opts+"\' is missing!"
			sys.exit()
	
	start=options["start"]
	end=options["end"]
	runNumber=options["runNumber"]
	outputdir=options["outputdir"]
	totsteps=options["totsteps"]
	nEvtsPerStep=options["nEvtsPerStep"]
	
	from Configurables import DDDBConf, CondDB, CondDBAccessSvc, NTupleSvc, EventClockSvc, Brunel, LHCbApp
#	if options.has_key("IgnoreHeartBeat"):
#		CondDB().IgnoreHeartBeat = options["IgnoreHeartBeat"]
		
	if options.has_key("addCondDBLayer"):
		altag="HEAD"
		if options.has_key("addCondDBLayer_tag"):
			altag=options["addCondDBLayer_tag"]
		CondDB().addLayer(CondDBAccessSvc("myCond",ConnectionString = "sqlite_file:"+options["addCondDBLayer"]+"/LHCBCOND",DefaultTAG = altag))

	# Need this line so as to not get db errors- should be fixed properly at some point
	CondDB().IgnoreHeartBeat = True
	CondDB().EnableRunStampCheck = False	

#	customDBs = glob.glob('/group/rich/ActiveDBSlices/*.db')
#	for db in customDBs:
#		CondDB().addLayer( CondDBAccessSvc(os.path.basename(db), ConnectionString="sqlite_file:"+db+"/LHCBCOND", DefaultTAG="HEAD") )
		
#	importOptions('$STDOPTS/DecodeRawEvent.py')
	#importOptions("$STDOPTS/RootHist.opts")
	#importOptions("$STDOPTS/RawDataIO.opts")
	#DEBUG by DisplayingHitMaps=False
	from Configurables import MDMRich1Algorithm
	mdmAlg=MDMRich1Algorithm("Rich1MDCS")
	mdmAlg.NumberOfEventsPerStep=nEvtsPerStep
	mdmAlg.StoreHistos=False
	mdmAlg.DEBUG = False
	
	if options.has_key("StoreHistos"):
		mdmAlg.StoreHistos = options["StoreHistos"]
	
	if options.has_key("DEBUG"):
		mdmAlg.DEBUG = options["DEBUG"]
	
	print "start step: "+str(start)
	print "stop step: "+str(end)
	print "processing "+str(nEvtsPerStep*(end-start))+" events"
	
	tuplestring="NTuple_Run%i_Steps%04d-%04d.root"%(runNumber,start,end)
		
	if options.has_key("TupleName"):
		tuplestring = options["TupleName"]
		
	histoname = "Histos_Run%i_Steps%04d-%04d.root"%(runNumber,start,end)
	
	if options.has_key("HistoName"):
		histoname = options["HistoName"]
	
	if outputdir!="":
		tuplestring = "%s/%s"%(outputdir,tuplestring)
		histoname = "%s/%s"%(outputdir,histoname)
	
	tuplename = "RICHTUPLE1 DATAFILE=\'%s\' TYP=\'ROOT\' OPT=\'NEW\'"%(tuplestring)

	# Currently put in manually. Edit here to use correct db and conddb tags 			
	LHCbApp().DDDBtag="dddb-20150724"	
	LHCbApp().CondDBtag="cond-20160123"
	
	if options.has_key("DDDBtag"):
		LHCbApp().DDDBtag= options["DDDBtag"]
	if options.has_key("CondDBtag"):
		LHCbApp().CondDBtag= options["CondDBtag"]
	
	#customDBs = glob.glob('/group/rich/ActiveDBSlices/*.db')
	#for db in customDBs:
	#	CondDB().addLayer( CondDBAccessSvc(os.path.basename(db), ConnectionString="sqlite_file:"+db+"/LHCBCOND", DefaultTAG="HEAD") )	
	
	ApplicationMgr().TopAlg += [mdmAlg]
	ApplicationMgr().ExtSvc += ['DataOnDemandSvc']
	ApplicationMgr().EvtMax = end*nEvtsPerStep

	# Timing information for application
	from Configurables import AuditorSvc, SequencerTimerTool
	ApplicationMgr().ExtSvc += ['AuditorSvc']
	ApplicationMgr().AuditAlgorithms = True
	AuditorSvc().Auditors += ['TimingAuditor']
	SequencerTimerTool().OutputLevel =  4

	LHCbApp().TimeStamp = True

	HistogramPersistencySvc().OutputFile = histoname
	NTupleSvc().Output= [tuplename]
	EventSelector().PrintFreq = 10
	EventSelector().PrintFreq=nEvtsPerStep
	EventSelector().FirstEvent=start*nEvtsPerStep
	print "First event: "+str(start*nEvtsPerStep)
	print "Last event: "+str(end*nEvtsPerStep)
	
	#get data, will look in local cluster first, then on castor
	isLocal=True
	if options.has_key("isLocal"):
		isLocal=options["isLocal"]
	DATA_and_Year=(getData(runNumber,start,end,totsteps,isLocal))

	DATA=DATA_and_Year["DATA"]
	if not len(DATA)>0:
		print "Data not found in local, switching to CASTOR"
		DATA_and_Year=getData(runNumber,start,end,totsteps,not isLocal)
		DATA=DATA_and_Year["DATA"]
	if not len(DATA)>0:
		print "DATA not found anywhere!"
		sys.exit()
	LHCbApp.DataType = str(DATA_and_Year["year"])
	
	EventSelector().Input=DATA
	EventClockSvc().EventTimeDecoder = "OdinTimeDecoder"
	
	appMgr = GaudiPython.AppMgr()
	appMgr.HistogramPersistency="ROOT"
	#appMgr.OutputLevel=DEBUG
	evtSvc= appMgr.evtSvc()
	
	esel = appMgr.evtsel()
	esel.PrintFreq=nEvtsPerStep
	appMgr.initialize()

	appMgr.run(nEvtsPerStep*(end-start))
	appMgr.stop()
	appMgr.finalize()
	

if __name__ == '__main__' :

	narg = sys.argv.__len__()
	if (narg<6):
		print "expecting: rawDataToNtuple.py <run number> <first step> <last step> <events per step> <tot steps>"
		sys.exit()
	runNumber=int(sys.argv[1])
	start=0
	end=999999
	nSteps=100000000
	nEvtsPerStep=6000
	start=int(sys.argv[2])
	end=int(sys.argv[3])
	nEvtsPerStep=int(sys.argv[4])
	totsteps=int(sys.argv[5])
	
	options={}
	options["runNumber"]=runNumber
	options["start"]=start
	options["end"]=end
	options["nEvtsPerStep"]=nEvtsPerStep
	options["totsteps"]=totsteps
	options["outputdir"]="/tmp/"+getpass.getuser()+"/"
	#options["CondDBtag"]="head-20110901"
	#options["DDDBtag"]="head-20110823"
	rawDataToNtuple(options)
	
