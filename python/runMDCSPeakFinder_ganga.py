from optparse import OptionParser, OptionGroup
import os, sys, getpass
	
parser = OptionParser()
parser.set_usage("python runMDCSPeakFinder_ganga.py -R RUNNUMBER -V VERSION [extra options]");
parser.add_option("-R", "--RunNumber", dest="--RunNumber", help="process RUNNUMBER", metavar="RUNNUMBER", type="int")
parser.add_option("-V", "--BrunelVersion", dest="--BrunelVersion", help="Brunel VERSION", metavar="VERSION", type="string")

group =OptionGroup(parser, "Other Options","These options may not be always necessary")
group.add_option("-D", "--DEBUG", action="store_true", dest="--DEBUG", help="print debug information")
group.add_option("-C","--dataToCastor", action="store_true", dest="outputdata", help="Store data to Castor")
group.add_option("-H", "--StoreHistos", action="store_true", dest="--StoreHistos", help="store root files with hitmap histograms information")
group.add_option("-B", "--Backend", dest="bkend",default="LSF",help="Choose ganga backend (available options are LSF at CERN and Interactive at both CERN and Online Cluster) [default: %default]",choices=["LSF","Interactive"],metavar="BACKEND")
group.add_option("-Q","--Queue",dest="queue",type="string", help="Set LSF queue [default = 8nh]")
group.add_option("--DDDBtag", dest="--DDDBtag", help="Set custom DDDB tag",type="string",metavar="TAG")
group.add_option("--CondDBtag", dest="--CondDBtag", help="Set custom CondDB tag",type="string",metavar="TAG")
group.add_option("--addCondDBLayer", dest="--addCondDBLayer", help="Add custom CondDB layer",type="string",metavar="LAYER")
group.add_option("--addCondDBLayer_tag", dest="--addCondDBLayer_tag", help="Set custom CondDB layer tag, if different from \"HEAD\"",type="string",metavar="TAG")
#group.add_option("-O", "--OutputDirectory", dest="--OutputDirectory", help="set output directory directory (you'll need some space..) [default: %default]",metavar="DIRECTORYPATH")
group.add_option("-N", "--TupleName", dest="--TupleName", help="set name of the output ntuple",metavar="FILENAME")
group.add_option("--HistoName", dest="--HistoName", help="set name of the output histogram file",metavar="FILENAME")
group.add_option("-T", "--TotSteps", dest="--TotSteps", help="total number of steps in the scan [default: %default]", metavar="NSTEPS", type="int",default=2200)
group.add_option("--EvtPerStep", dest="--EvtPerStep", help="number of events in each step [default: %default]", metavar="NEVENTS", type="int",default=6000)
group.add_option("--IgnoreHeartBeat", action="store_true", dest="--IgnoreHeartBeat", help="Set IgnoreHeartBeat=True in CondDB")

group_one =OptionGroup(parser, "1 Job Options","You must provide these options in case you do not want to process the entire scan but just one subjob")
group_one.add_option("-S", "--StartStep", dest="--StartStep", help="first calibration step to process", metavar="STEP", type="int")
group_one.add_option("-E", "--EndStep", dest="--EndStep", help="last calibration step to process", metavar="STEP", type="int")

parser.add_option_group(group)
parser.add_option_group(group_one)


flag_list=["--DEBUG","outputdata","--StoreHistos","--IgnoreHeartBeat"]
	
(options, args) = parser.parse_args()

dict_opts=options.__dict__

#print "options"
#print dict_opts

if dict_opts["--RunNumber"]==None:
	print "Please provide a run number!"
	print "run: \"python runMDCSPeakFinder.py -h\" for help"
	sys.exit()

if dict_opts["--BrunelVersion"]==None:
	print "Please provide a Brunel version!"
	print "run: \"python runMDCSPeakFinder.py -h\" for help"
	sys.exit()
#clean dictionary from non-specified options
blacklist=[]
for opts, val in dict_opts.iteritems():
	if dict_opts[opts]==None:
		blacklist.append(opts)
for bkeys in blacklist:
	del dict_opts[bkeys]

if dict_opts.has_key("--EndStep") and dict_opts.has_key("--StartStep"):
	print "Running only one subjob"
	
else:
	print "Running all steps"
	if not dict_opts.has_key('--TupleName'):
		dict_opts["--TupleName"]="NTuples_Run%i_AllSteps.root"%(dict_opts["--RunNumber"])
	if not dict_opts.has_key('--HistoName'):
		dict_opts["--HistoName"]="Histos_Run%i_AllSteps.root"%(dict_opts["--RunNumber"])

if not dict_opts.has_key('--TupleName'):
	dict_opts["--TupleName"]="MDCS_ganga_tuple_Run%i.root"%(dict_opts["--RunNumber"])
if not dict_opts.has_key('--HistoName'):
	dict_opts["--HistoName"]="MDCS_ganga_histo_Run%i.root"%(dict_opts["--RunNumber"])

#print dict_opts
#build ganga job
#app=GaudiPython(project="Brunel",version=dict_opts["--BrunelVersion"])
app=GaudiPython()
app.user_release_area=os.environ['User_release_area']
app.project="Brunel"
app.version=dict_opts["--BrunelVersion"]
#app.script = [os.path.dirname(os.path.abspath(sys.argv[0]))+"/runMDCSPeakFinder.py"]
app.script = ["/afs/cern.ch/user/a/anandi/cmtuser/Brunel_v47r5/Rich/Rich1MDCS/python/runMDCSPeakFinder.py"]

for opts, val in dict_opts.iteritems():
	if opts!="--BrunelVersion":
		if opts not in flag_list:
			app.args.append(str(opts)+"="+str(val))
		else:
			app.args.append(str(opts))

j = Job(application=app)
j.name="RICH1MDCS_PeakFinder"
if dict_opts.has_key("outputdata"):
	if dict_opts["outputdata"]==True:
		j.outputdata=[str(dict_opts["--TupleName"]),str(dict_opts["--HistoName"])]

the_queue="8nh"
if dict_opts.has_key("queue"):
	the_queue=dict_opts["queue"]

j.backend = LSF( queue = the_queue)
if dict_opts["bkend"]=="Interactive":
	j.backend = Interactive()
if dict_opts["bkend"]=="LSF":
	j.backend = LSF(queue = the_queue)
#j.submit()
