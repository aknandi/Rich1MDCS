#Donal Hill 9 Jan 2016
Usage of the RICH1MDCS software


The MDCS analysis is divided into 3 steps:
1) Peak finding
2) Pattern cleaning
3) Calibration parameters fitting


1. Peak finding

A Brunel algorithm (MDMRich1Algorithm) runs on raw data (the famous MDCS scan).
Peaks candidates are formed from hits in the HPD.

How to run: GuadiPython script under python/runMDCSParamFitter.py

The script can be run on the online cluster(recommended) or lxplus. 
It will search for data in the local disk first (online cluster only), if the data is not found
it will search on CASTOR

The script supports several options, just run
"python runMDCSParamFitter.py --help"
to get all of them.
Some are mandatory, others are optional.
A "standard" jobs (read "what you should try first") will look like:

"python runMDCSPeakFinder.py --RunNumber=94224"

This will run on all steps and give an ntuple to be processed in the cleaning. The run number here is the run number of the MDMS scan, which you can find from the LHCb logbook.

If necessary only the needed steps can be processed:
python runMDCSPeakFinder.py --RunNumber=94224 --StartStep=2100 --StopStep=2200

but I have never needed to do this.

Database tags (and even CondDB layer) can be passed explicitly:
"python runMDCSPeakFinder.py --RunNumber=94224 --DDDBtag=head-20110823 --CondDBtag=head-20110901 --addCondDBLayer=/pathtothelayer/layer.db"
This may be needed if you run on very old or when some hardware has been changed recently.
As a debugging tool you can set the DEBUG verbose mode with the "-D" option and the storage of hitmap histograms with "-H" (I suggest to store hitmaps only
if you run for<=100 steps as you may run out of memory)

Another scripts accepts almost the same options but submits the jobs to LSF (lxplus only) through ganga, run with
ganga runMDCSPeakFinder_ganga.py
This maybe useful for old jobs, running on castor is pretty slow so you may just send the jobs and 
wait for ganga to tell you that it is finished.

This option was old advice from Andrea, and I never used it.


2. Pattern cleaning

This is a simple ROOT/rooFit C++ script. It is compiled automatically when you "make" the package.
The executable is under scripts/PatternCleaning.exe
The default output filename of the peak finder is something like "NTuples_Run92442_AllSteps.root"
Then you can clean the pattern with
"PatternCleaning.exe SOMEPATH/NTuples_Run92442_AllSteps.root"
After ~3 minutes you will get 2 ntuples, 1 for each HPD box (U or D):
NTuple_Run92442_afterProcess_U.root
NTuple_Run92442_afterProcess_D.root

These ntuples contain various control plots and an ntuple with the cleaned peak that will
be further processed in 3.
After a scan is taken, you will be asked to check whether it is good or not by the RICH guys. In order to quickly check that
everything is ok, check the canvas "can_All_Peaks" in the root file. All the peaks found on each HPD should appear here. There may be missing hpds.
This may either be caused by wrong DB tags during step 1 or simply because the are disabled for
whatever reason. Nothing to worry about in the latter case, otherwise go back to step 1 logs and see
if there is something wrong.

You can also look at the pattern after the cleaning "can_Final"
If everything looks normal proceed to step 3, otherwise you may want to check that the fits to the distributions
make sense:
"rooFit_xdiff"
"rooFit_prop"
"rooFit_npes"
"rooFit_max"
Also the likelihood distribution may be useful:
"DLL_tot"
You may also change the DLL cut by changing the value of the global variable "double DLL_CUT=-2"
at the beginning of scripts/Functions4Cleaning.C (and re-compile of course)
I either use 0 or -2 depending on the scan.


3. Calibration parameters fitting
Again you need to run a GaudiPython script: python/runMDCSParamFitter.py
The algorithm called is "MDMFitParameters"
Upper and lower box are processed separately (option --Rich1Panel):
"python runMDCSParamFitter.py --RunNumber=92442 --MagPolarity=Down --Rich1Panel=D --InputFile=~/NTuple_Run92442_afterProcess_D.root"
"python runMDCSParamFitter.py --RunNumber=92442 --MagPolarity=Down --Rich1Panel=U --InputFile=~/NTuple_Run92442_afterProcess_D.root"
Just run "python runMDCSParamFitter.py --help" to list all options.

For my work I added some more options to this, which you need to add to the command line for it to run properly. You can just use the values below here, and don't worry about what they are set to. I added them in to experiment with the code, but they don't do anything important to the output:

--RandomSeed=0 --RadiusFactor=1 --Systematics=0
 

You will get a root file with control plots (way more than what you need) and an 4 xml files, two with the new and two with the previus parameters (I would ignore the previous one, as I don't think it gets filled reliably!)
There is one xml file for each HPD box.

In order to produce a full set of parameters for the CondDB you need to process two scans of opposite polarity (Up or Down) and
copy and paste the content of the xml files to produce a "HPD.xml" file that you may then give to some rich guy (usually Chris Jones).
You can also check them by yourself with a Brunel job in which you re-run the rich reconstruction with the newly calculated parameters.

 I think there's an example HPD.xml in doc/, where you will see for distinct blocks of parameters. Mag Up U box, Mag Up D box, Mag Down U box and Mag Down D box.

I have code that makes HPD.xml files into a DB slice under the "CreateDBLayer" folders, but I think it's maybe better if you ask Chris and Antonis for the latest code to do this. You can send them the HPD.xml perhaps and get them to make the DB slice. I think this is definitely safer in the beginning, as they can also assess the MDMS parameters within 

Sorry for the delay and I'll try my best to answer any issues via email if I can!

Donal (9/1/16)

