#include "Functions4Cleaning.C" 
#include "mystyle.h"
#include "TPRegexp.h"
#include "TObjArray.h"

void CleanRun(TString box_string, TString inputfile, TString outputfile){
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TFile *infile=NULL;
	infile=new TFile(Form("%s",inputfile.Data()));
	
	TTree *tree=(TTree*)infile->Get("RICH/ToolSvc.MDMHitManager/Map");
	
	TFile *outfile=new TFile(Form("%s",outputfile.Data()),"RECREATE");
	int thebox=-1;
	if(box_string=="U") thebox=0;
	else thebox=1;
	
	gStyle->SetOptTitle(1);
	
	//fill hpd map
	TList *HPDMap=new TList();
	FillHPDMap(HPDMap,tree,thebox);
	
	//get initial peakmap
	AllPeaksPlot(HPDMap, thebox);
	OnePeakPerStep(HPDMap, thebox);
	
	//flag IFB-affected hpds
	double cutoff[2];
	double sigmas=2;
	TH2I* h2=GetHPDMapIFB(HPDMap, cutoff, box_string);
	
	FlagIFBHPDs(HPDMap, cutoff, sigmas);
	
	cout << "Ion feedback map ready"<<endl;
// 	outfile->cd("..");
// 	//create distributions
	FillDistVar(HPDMap);
	FillXDiff(HPDMap);
	CreateDistributions(HPDMap);
	CreateDistributions_trees(HPDMap);
// 	//create ll
	BuildLL(HPDMap,outfile);
	CreateDLLDistributions(HPDMap,thebox);
	StripMultiPeakSteps(HPDMap,thebox);
	
	FinalTree(HPDMap,thebox);
	
	//other plots
	NPeaksPerHPD(HPDMap);
	LightSpotFeatures(HPDMap);
	PlotHPDsPatterns(HPDMap);
	
	outfile->Write();
	outfile->Close();
	
	delete tree;
	infile->Close();
// 	exit(1);
}


int main(int argc, const char* argv[]){
	mystyle();
// 	cout << argc << endl;
	if(argc<2){
		cout << "Please give at least an input file!"<<endl;
		exit(1);
	}
	TString s_infile=argv[1];
	TString string_up, string_down;
	if(argc==2){
		TPRegexp guess_run("Run[0-9]*");
		TObjArray* subStrL=guess_run.MatchS(s_infile.Data());
		cout << "Guessing the scan number..."<<endl;
		TString s_run;
		if((subStrL->GetEntries())>0){
			s_run = ((TObjString *)subStrL->At(0))->GetString();
			cout << "Processing "<<s_run << endl;
		}
		string_up=Form("NTuple_%s_afterProcess_%s.root",s_run.Data(),"U");
		string_down=Form("NTuple_%s_afterProcess_%s.root",s_run.Data(),"D");
		
		cout << "The output files will be: "<<endl;
		cout << string_up << endl;
		cout << string_down << endl;
		cout << endl;
		cout << "Processing.."<<endl;
		CleanRun("U",s_infile,string_up);
		CleanRun("D",s_infile,string_down);
// 		cout << nmatches << endl;
// 		cout << sub_s.Data() << endl;
	}
// 	else{
// 		
// 		FILE * pFile;
// 		pFile = fopen (Form("/scratch/rich/CASTOR/MDMS/Raw_Processed/NTuples_Run%i_AllSteps.root",s_infile.Data()),"r");
// 		if (pFile==NULL)
// 		{
// 			cout << "Error opening file!" << endl;
// 			cout << "You probably need to run on raw data first."<<endl;
// 			exit (1);
// 		}
// 		else
// 		{
// 			
// 			
// 	// 		cout << s_run<<endl;
// 	// 		cout << s_run.Atoi()<<endl;
// 			CleanRun(s_run.Atoi(),"D"); //lower box
// 			CleanRun(s_run.Atoi(),"U"); //upper box
// 		}
// 	}
	
}
