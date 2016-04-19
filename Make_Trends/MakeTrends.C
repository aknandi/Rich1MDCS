#include <algorithm>
//#include "/home/andrea/mystyle.h"
#include <math.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TROOT.h"
#include "TString.h"
#include "TVector3.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TF1.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <utility>
#include "TTimeStamp.h"
#include "TDatime.h"
#include "TMap.h"
#include <map>

using namespace std;

double r_min[4]={-0.1,5,-0.2,-0.05};
double r_max[4]={0.1,5.8,0.2,0.05};
double a_min[4]={-0.5,-0.3,-0.02,-0.0005};
double a_max[4]={0.5,0.3,0.02,0.0005};

map<int,TTimeStamp> GetMap(){
	map<int,TTimeStamp> tmap;
	/*tmap[65428] = TTimeStamp(2010,2,8,0,0,0);
	tmap[53625] = TTimeStamp(2009,7,25,0,0,0);
	tmap[53652] = TTimeStamp(2009,7,27,0,0,0);
	tmap[69453] = TTimeStamp(2010,4,2,0,0,0);
	tmap[69454] = TTimeStamp(2010,4,2,12,0,0);
	tmap[70254] = TTimeStamp(2010,4,16,0,0,0);
	tmap[70748] = TTimeStamp(2010,4,26,0,0,0);
	tmap[72065] = TTimeStamp(2010,5,20,0,0,0);
	tmap[78335] = TTimeStamp(2010,8,30,0,0,0);
	tmap[78347] = TTimeStamp(2010,8,30,12,0,0);
	tmap[78512] = TTimeStamp(2010,9,2,0,0,0);
	tmap[78656] = TTimeStamp(2010,9,4,0,0,0);
	tmap[62152] = TTimeStamp(2009,10,18,0,0,0);
	tmap[65402] = TTimeStamp(2010,2,8,0,0,0);
	tmap[63489] = TTimeStamp(2009,12,6,0,0,0);
	tmap[92442] = TTimeStamp(2011,05,27,12,0,0);*/
	tmap[113865] = TTimeStamp(2012,04,27,12,0,0);
	tmap[136273] = TTimeStamp(2013,01,30,12,0,0);
	tmap[136319] = TTimeStamp(2013,01,31,12,0,0);
	
	return tmap;
}

double getrundate(int run){
	map<int,TTimeStamp> tmap=GetMap();
	return tmap[run].AsDouble();
}
TString getrundate_string(int run){
	map<int,TTimeStamp> tmap=GetMap();
	TDatime dt(tmap[run].AsString("s"));
	TString ts=Form("%i/%i/%i",dt.GetDay(),dt.GetMonth(),dt.GetYear());
	return ts;
}

void GetRuns(TString polarity, TString Box, vector<TString> *files, vector<int> *run_number){
	
	DIR *dp;
	TString dir=".";
	struct dirent *dirp;
	if((dp  = opendir(dir.Data())) == NULL) {
		cout << "Error(" << errno << ") opening " << dir << endl;
		//              return errno;
	}
	
	while ((dirp = readdir(dp)) != NULL) {
		TString tmps=string(dirp->d_name);
		TString regex="";
		regex="\\bMDMSOutput_Run([0-9]*)_Magnet"+polarity+"_"+Box+".root\\b";
// 		cout << tmps << endl;
// 		cout << regex << endl;
		TPRegexp myregex(regex.Data());
		if(tmps.Contains(myregex)){
// 			cout << "yeah"<< endl;
			TObjArray *subStrL = myregex.MatchS(tmps);
			const int nrSubStr = subStrL->GetLast()+1;
// 			cout << nrSubStr << endl;
			if (nrSubStr > 1) {
				const TString scut1 = ((TObjString *)subStrL->At(1))->GetString();
				if(scut1.Atoi()!=78335 || 1){
					(*files).push_back(tmps);
	// 				cout << scut1.Atof() << endl;
					(*run_number).push_back(scut1.Atoi());
				}
			}
		}
	}
	
	cout << files->size() << " of polarity "<<polarity<< " for " << Box <<" Box found:" <<endl;
	sort(run_number->begin(),run_number->end());
	sort(files->begin(),files->end());
	int i=0;
	for(i=0;i<files->size();i++) cout << (*files)[i] << endl;
	for(i=0;i<files->size();i++) cout << (*run_number)[i] << endl;
	
	closedir(dp);
}

double getDatapoint(TString file, int hpd, int col, TString param=""){
	TFile *infile=new TFile(file.Data());
	TTree *tree=(TTree*)infile->Get("parameters");
	
	double par;
	int tmpcol,tmphpd;
	tree->SetBranchAddress(param.Data(),&par);
	tree->SetBranchAddress("hpd",&tmphpd);
	tree->SetBranchAddress("col",&tmpcol);
	unsigned int i=0;
// 	cout << "new"<<endl;
	double var=-1;
// 	cout << "hc " <<hpd << "\t"<<col<<endl;
	for(i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);
// 		cout << tmpcol << "\t"<<tmphpd<<endl;

		if(tmpcol==col && tmphpd==hpd){
// 			cout << "Datapoint found"<<endl;
			var=par;
			break;
		}
	}

	infile->Close();
// 	delete tree;
// 	delete infile;
// 	cout << var << endl;
	return var;
}

TF1* getFunction(TString file, int hpd, int col, TString type, int runnumber){
	TFile *infile=new TFile(file.Data());
	TTree *tree=(TTree*)infile->Get("parameters");
	
	TF1* f1=new TF1();
	if(type=="r"){
		f1=new TF1(Form("%s_func_hpd%i_col%i_%i",type.Data(),hpd,col,runnumber),"[0]+x*[1]+x*x*[2]+x*x*x*[3]",0,40);
		f1->SetTitle(Form("%s_func_hpd%i_col%i;R_{C} (mm);#Delta R (mm)",type.Data(),hpd,col));
		f1->SetMaximum(5);
		f1->SetMinimum(-5);
	}
	if(type=="a"){
		f1=new TF1(Form("%s_func_hpd%i_col%i_%i",type.Data(),hpd,col,runnumber),"x*[0]+x*x*[1]+x*x*x*[2]+x*x*x*x*[3]",0,40);
		f1->SetTitle(Form("%s_func_hpd%i_col%i;R_{C} (mm);#Delta #phi * R_{C} (mm)",type.Data(),hpd,col));
		f1->SetMaximum(5);
		f1->SetMinimum(-5);
	}
	
	f1->SetParameter(0,0);
	f1->SetParameter(1,0);
	f1->SetParameter(2,0);
	f1->SetParameter(3,0);
	
	double par_r[4];
	double par_a[4];
	int tmpcol,tmphpd;
	int k=0;
	for(k=0;k<4;k++) tree->SetBranchAddress(Form("r%i",k),&par_r[k]);
	for(k=0;k<4;k++) tree->SetBranchAddress(Form("a%i",k),&par_a[k]);
	tree->SetBranchAddress("hpd",&tmphpd);
	tree->SetBranchAddress("col",&tmpcol);
	unsigned int i=0;
// 	cout << "new"<<endl;
	double var=-1;
// 	cout << "hc " <<hpd << "\t"<<col<<endl;
	for(i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);
// 		cout << tmpcol << "\t"<<tmphpd<<endl;

		if(tmpcol==col && tmphpd==hpd){
			if(type=="r"){
				f1->SetParameter(0,par_r[0]);
				f1->SetParameter(1,(par_r[1]-1./0.18)*0.18);
				f1->SetParameter(2,par_r[2]*(pow(0.18,2)));
				f1->SetParameter(3,par_r[3]);
			}
			if(type=="a"){
				
				double a1=par_a[1]/par_r[1];
				double a2=(par_a[2]-a1*par_r[2])/(pow(par_r[1],2));
				double a3=(par_a[3]-a1*par_r[3]-2*a2*par_r[1]*par_r[2])/(pow(par_r[1],3));
				
				f1->SetParameter(0,par_a[0]);
				f1->SetParameter(1,a1);
				f1->SetParameter(2,a2);
				f1->SetParameter(3,a3);
			}
			break;
		}
	}

	infile->Close();
// 	delete tree;
// 	delete infile;
// 	cout << var << endl;
	return f1;
	
}

void MakeTrends(TString polarity="Off", TString Box="U"){
	vector<TString> files; 
	vector<int> run_number;
	GetRuns(polarity,Box,&files,&run_number);
	
	TFile *outfile=new TFile(Form("MDMSTrends_Magnet%s_%s.root",polarity.Data(),Box.Data()),"RECREATE");
	
	TGraph gr[4][14][7];
	TGraph ga[4][14][7];
	TF1 r_fun[10][14][7];
	TF1 a_fun[10][14][7];
	
// 	cout << "crashing"<< endl;
	
	TLegend *l=new TLegend(0.6,0.6,0.9,0.9);
	l->SetFillStyle(0);
	l->SetBorderSize(0);
	
	int i=0,k=0,h=0,c=0;
	for(k=0;k<4;k++){
		for(h=0;h<14;h++){
			for(c=0;c<7;c++){
				gr[k][h][c]=TGraph(files.size());
				gr[k][h][c].SetName(Form("Trend_radial%i_hpd%i_col%i",k,h,c));
				gr[k][h][c].SetTitle(Form("Trend_radial%i_col%i",k,c));
				gr[k][h][c].SetMarkerStyle(h+20);
				if((h+1)%10!=0){
					gr[k][h][c].SetMarkerColor(h+1);
					gr[k][h][c].SetLineColor(h+1);
				}
				else{
					gr[k][h][c].SetMarkerColor(h+10);
					gr[k][h][c].SetLineColor(h+10);
				}
				gr[k][h][c].SetLineWidth(2);
				ga[k][h][c]=TGraph(files.size());
				ga[k][h][c].SetName(Form("Trend_axial%i_hpd%i_col%i",k,h,c));
				ga[k][h][c].SetTitle(Form("Trend_axial%i_col%i",k,c));
				ga[k][h][c].SetMarkerStyle(h+20);
				ga[k][h][c].SetMarkerColor(h+1);
				ga[k][h][c].SetLineColor(h+1);
				ga[k][h][c].SetLineWidth(2);
				if(k==0 && c==0){
					l->AddEntry(&gr[k][h][c],Form("HPD %i",h),"lp");
				}
			}
		}
	}
// 		TGraph *g=TGraphparname_a=Form("a%i",k);
// 		TString parname_r=Form("r%i",k);
	
	
	
	for(i=0;i<files.size();i++){
		for(k=0;k<4;k++){
			for(h=0;h<14;h++){
				for(c=0;c<7;c++){
					TString parname_a=Form("a%i",k);
					TString parname_r=Form("r%i",k);
					gr[k][h][c].SetPoint(i,double(getrundate(run_number[i])),getDatapoint(files[i], h, c, parname_r));
					ga[k][h][c].SetPoint(i,double(getrundate(run_number[i])),getDatapoint(files[i], h, c, parname_a));
					gr[k][h][c].GetYaxis()->SetTitle(Form("#rho_{%i}",k));
					ga[k][h][c].GetYaxis()->SetTitle(Form("#theta_{%i}",k));
					
				}
			}
		}
	}
	
	for(k=0;k<4;k++){
		for(h=0;h<14;h++){
			for(c=0;c<7;c++){
				outfile->Add(&gr[k][h][c]);
				outfile->Add(&ga[k][h][c]);
			}
		}
	}
	
	TLegend *lfun=new TLegend(0.6,0.1,0.9,0.4);
	lfun->SetFillStyle(0);
	lfun->SetBorderSize(0);
	
	for(i=0;i<files.size();i++){
		for(h=0;h<14;h++){
			for(c=0;c<7;c++){
				r_fun[i][h][c]= *getFunction(files[i], h, c, "r",run_number[i]);
				r_fun[i][h][c].SetLineColor(i+1);
				if(i+1==1) r_fun[i][h][c].SetLineColor(i+1+9);
				if(i+1==5) r_fun[i][h][c].SetLineColor(i+1+9);
				r_fun[i][h][c].SetLineStyle(i+1);
				r_fun[i][h][c].SetLineWidth(4);
				a_fun[i][h][c]= *getFunction(files[i], h, c, "a",run_number[i]);
				a_fun[i][h][c].SetLineColor(i+1);
				if(i+1==1) a_fun[i][h][c].SetLineColor(i+1+18);
				if(i+1==5) a_fun[i][h][c].SetLineColor(i+1+18);
				a_fun[i][h][c].SetLineStyle(i+1);
				a_fun[i][h][c].SetLineWidth(4);
				outfile->Add(&r_fun[i][h][c]);
				outfile->Add(&a_fun[i][h][c]);
			}
		}
		if(run_number[i]!=92442) lfun->AddEntry(&r_fun[i][0][0],Form("%s - 18KV",getrundate_string(run_number[i]).Data()),"l");
		else lfun->AddEntry(&r_fun[i][0][0],Form("%s - 15KV",getrundate_string(run_number[i]).Data()),"l");
	}
	
	for(k=0;k<4;k++){
		for(c=0;c<7;c++){
			TCanvas *tmpcan_r=new TCanvas(Form("Canvas_r%i_col%i",k,c),Form("Canvas_r%i_col%i",k,c));
			for(h=0;h<14;h++){
				if(h==0) {
					gr[k][h][c].SetMinimum(r_min[k]);
					gr[k][h][c].SetMaximum(r_max[k]);
					gr[k][h][c].Draw("APL");
					gr[k][h][c].GetXaxis()->SetTimeDisplay(1);
					gr[k][h][c].GetXaxis()->SetNdivisions(-505);
					gr[k][h][c].GetXaxis()->SetTimeFormat("%m/%y");
					gr[k][h][c].GetXaxis()->SetTimeOffset(0,"gmt");
					gr[k][h][c].GetYaxis()->SetTitle(Form("#rho_{%i}",k));
				}
				else gr[k][h][c].Draw("PL");
			}
			l->Draw();
			outfile->Add(tmpcan_r);
			TCanvas *tmpcan_a=new TCanvas(Form("Canvas_a%i_col%i",k,c),Form("Canvas_a%i_col%i",k,c));
			for(h=0;h<14;h++){
				if(h==0) {
					ga[k][h][c].SetMinimum(a_min[k]);
					ga[k][h][c].SetMaximum(a_max[k]);
					ga[k][h][c].Draw("APL");
					ga[k][h][c].GetXaxis()->SetTimeDisplay(1);
					ga[k][h][c].GetXaxis()->SetNdivisions(-505);
					ga[k][h][c].GetXaxis()->SetTimeFormat("%m/%y");
					ga[k][h][c].GetXaxis()->SetTimeOffset(0,"gmt");
					ga[k][h][c].GetYaxis()->SetTitle(Form("#theta_{%i}",k));
				}
				else ga[k][h][c].Draw("PL");
			}
			l->Draw();
			outfile->Add(tmpcan_a);
		}
	}
	for(c=0;c<7;c++){
		for(h=0;h<14;h++){
			TCanvas *tmpcanf_a=new TCanvas(Form("Canvas_fa_col%i_hpd%i",c,h),Form("Canvas_fa_col%i_hpd%i",c,h));
			for(i=0;i<files.size();i++){
				if(i==0)a_fun[i][h][c].Draw();
				else a_fun[i][h][c].Draw("SAME");
			}
			lfun->Draw();
			outfile->Add(tmpcanf_a);
			TCanvas *tmpcanf_r=new TCanvas(Form("Canvas_fr_col%i_hpd%i",c,h),Form("Canvas_fr_col%i_hpd%i",c,h));
			for(i=0;i<files.size();i++){
				if(i==0)r_fun[i][h][c].Draw();
				else r_fun[i][h][c].Draw("SAME");
			}
			outfile->Add(tmpcanf_r);
			lfun->Draw();
		}
	}
	
	outfile->Write();
	outfile->Close();
	delete outfile;
	
}

int main(){
  //mystyle();
	MakeTrends("Off","U");
	MakeTrends("Off","D");
	MakeTrends("Down","U");
	MakeTrends("Down","D");
	MakeTrends("Up","U");
	MakeTrends("Up","D");
	return 0;
}
