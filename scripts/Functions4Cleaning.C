#include "CleanHPD.h"
#include "CleanHPD.cpp" 
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TArc.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TString.h"
#include "TNtuple.h"

#define nCOLs 7
#define nHPDs 14
#define nSTEPs 2200

using namespace RooFit;
using namespace std;

// int nCOLs=7;

double DLL_CUT=-2;

void FillHPDMap(TList *HPDMap,TTree *tree, int box){
	int i=0,k=0;
	for(i=0; i<nCOLs; i++){
		for(k=0; k<nHPDs; k++){
			HPDMap->Add(new CleanHPD(box,i,k,tree));
		}
	}
// 	return HPDMap;
}

void produceplot(TString title, TGraph* data, int thebox){
	int NHPD=14;
	int NCOL=7;
	float VSpace=28;
	float HTilt=16;
	float MagFactor=2;
	TArc *arcs=new TArc[NHPD*NCOL];
	
	data->SetMarkerStyle(6);
	data->SetMarkerSize(1);
	data->SetMarkerColor(4);
	data->SetLineColor(0);
	data->SetTitle("");
	data->SetMaximum(VSpace*7);
	data->SetMinimum(0);
	data->GetXaxis()->SetLimits(0,NHPD*32+HTilt);
	data->GetYaxis()->SetLimits(0,VSpace*NCOL+20);
	float width=(NHPD*32+HTilt)*MagFactor;
	float height=(VSpace*7)*MagFactor+20;
	float canmagfactor=1.3;
	
	TCanvas *can=new TCanvas("can_"+title,title,width*canmagfactor,height*canmagfactor);
	can->SetWindowSize(width*canmagfactor,height*canmagfactor);
	can->SetTopMargin(0.05);
	can->SetBottomMargin(0.05);
	can->SetLeftMargin(0.02);
	can->SetRightMargin(0.02);
	data->GetXaxis()->SetAxisColor(0);
	data->GetYaxis()->SetAxisColor(0);
	data->GetXaxis()->SetLabelColor(0);
	data->GetYaxis()->SetLabelColor(0);
	data->Draw("AP");
	
	//draw arcs
	int i,k,j;
	for(i=0;i<NHPD;i++){
		for(k=0;k<NCOL;k++){
			arcs[i+k].SetFillColor(19);
			//                      arcs[i+k].SetFillStyle(300);
			arcs[i+k].SetLineWidth(1);
			arcs[i+k].SetLineColor(13);
			// x centre, y centre, radius with default phimin=0 and phimax=360, i.e. a circle
			if(thebox==0) arcs[i+k].DrawArc(16+32*(13-i)+HTilt*(k%2),16+VSpace*(6-k),15.7);
			if(thebox==1) arcs[i+k].DrawArc(16+32*(i)+HTilt*((1+k)%2),16+VSpace*(k),15.7);
		}
	}
	data->Draw("P SAME");
	gPad->SetTicks(1);
	gPad->SetLineColor(0);
	gVirtualX->SetFillColor(0);
	gVirtualX->SetFillStyle(0);
	gVirtualX->SetLineColor(0);
	gVirtualX->SetLineWidth(0);
	can->Write();
// 	cout << thebox << endl;
	title+=Form("_box%i",thebox);
	data->Write(title+"_graph");
// 	gFile->WriteTObject(can);
// 	can->SaveAs(title+"_hpdplane.png");
}

void FitMax(TTree *t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=10;
	double high=140;
	RooRealVar *m=new RooRealVar("max","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_max","Mean of Gaussian",60.0,40.0,70.0);
	RooRealVar *sigma=new RooRealVar("sigma_max","Width of Gaussian",10.0,3.0,40.0);
	RooGaussian *gauss=new RooGaussian("gauss_max","gauss",*m,*mean,*sigma);
	//pol2bkg====
	RooRealVar *pol1=new RooRealVar("p1_max","p1_max",-2.0,-10.0,-0.5);
	RooRealVar *pol2=new RooRealVar("p2_max","p2_max",0.0,-2.0,-1.0);
	RooRealVar *pol3=new RooRealVar("p3_max","p3_max",0.0,-100.0,100.0);
// 	RooRealVar *pol4=new RooRealVar("p4_max","p4_max",0,-100,100);
// 	RooRealVar *pol5=new RooRealVar("p5_max","p5_max",0,-100,100);
// 	RooRealVar *pol6=new RooRealVar("p6_max","p6_max",0,-100,100);
// 	RooRealVar *pol7=new RooRealVar("p7_max","p7_max",0,-100,100);
// 	RooRealVar *pol8=new RooRealVar("p8_max","p8_max",0,-100,100);
// 	RooPolynomial *pol=new RooPolynomial("pol_max","pol_max",*m,*pol1);
// 	
	
	RooRealVar *exp=new RooRealVar("p1_max","p1_max",-0.04,-0.055,-0.015);
	RooExponential *pol=new RooExponential("pol_max","pol_max",*m,*exp);
	
	RooRealVar *nsig=new RooRealVar("nS_max","signal fraction",1000,500,100000);
	RooRealVar *nbkg=new RooRealVar("nB_max","background fraction",200,0,100000);
	RooRealVar *soverb=new RooRealVar("SoverB_max","background fraction",0.7,0.6,1.0);
	RooAddPdf *model=new RooAddPdf("model_max","model_max",*gauss,*pol,*soverb);
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_max","data",t1,RooArgList(*m));
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_max","rooFit_max");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
// 	model->paramOn(mframe,Label("Fit output: Gauss + Pol"));
// 	model->plotOn(mframe,VisualizeError(*r,1),FillColor(kRed));
	model->plotOn(mframe,Name("model_max"),LineColor(kBlue));
	model->plotOn(mframe,Name("Signal Gaussian_max"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Polynomial background_max"),Components(*pol),LineColor(kGreen));
	data->plotOn(mframe);
	mframe->SetName("mframe_max");
	mframe->Draw();
	mframe->SetXTitle("M");
	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("Signal Gaussian_max","Signal", "L");
	leg->AddEntry("Polynomial background_max","Background", "L");
	leg->Draw();
	can->Write();
	model->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	pol1->Write();
	pol2->Write();
	pol3->Write();
// 	exp->Write();
	pol->Write();
// 	can->SaveAs("testmax.eps","eps");
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_max_sig");
	TF1* f2 = pol->asTF( RooArgList(*m) );
	f2->Write("tf1_max_bkg");
// 	double xmax = f->GetMaximumX(); 
}


void FitRot(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=-0.4;
	double high=0.4;
	RooRealVar *m=new RooRealVar("rot","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_rot","Mean of Gaussian",0,-0.001,0.001);
	RooRealVar *sigma=new RooRealVar("sigma_rot","Width of Gaussian",0.001,0,0.01);
	RooGaussian *gauss=new RooGaussian("gauss_rot","gauss",*m,*mean,*sigma);
	//Gauss bkg====
	RooRealVar *mean_bkg=new RooRealVar("mean_rot_bkg","Mean of bkg Gaussian",0,-0.001,0.001);
	RooRealVar *sigma_bkg=new RooRealVar("sigma_rot_bkg","Width of bkg Gaussian",0.3,0.01,0.5);
	RooGaussian *gauss_bkg=new RooGaussian("gauss_rot_bkg","gauss bkg",*m,*mean,*sigma);
	
	RooRealVar *nsig=new RooRealVar("nS_rot","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_rot","background fraction",500,0,100000);
	RooAddPdf *model=new RooAddPdf("model_rot","model_rot",RooArgList(*gauss,*gauss_bkg),RooArgList(*nsig,*nbkg));
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_rot","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_rot","rooFit_rot");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
	model->paramOn(mframe);
	model->plotOn(mframe,VisualizeError(*r,1),LineColor(kCyan));
	model->plotOn(mframe,Name("model_rot"),LineColor(kBlue));
	model->plotOn(mframe,Name("Signal Gaussian_rot"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Gaussian background_rot"),Components(*gauss_bkg),LineColor(kGreen));
	mframe->SetName("mframe_rot");
	mframe->Draw();
	mframe->SetXTitle("Peak rot");
	can->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	model->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	mean_bkg->Write();
	sigma_bkg->Write();
	gauss_bkg->Write();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_rot_sig");
	TF1* f2 = gauss_bkg->asTF( RooArgList(*m) );
	f2->Write("tf1_rot_bkg");
// 	return model;
	//      frac.Print();
}

void FitXDiff(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=0;
	double high=25;
	RooRealVar *m=new RooRealVar("xdiff","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_xdiff","Mean of Gaussian",0.0,0.0,0.001);
	RooRealVar *sigma=new RooRealVar("sigma_xdiff","Width of Gaussian",2.0,0.5,2.0);
	RooGaussian *gauss=new RooGaussian("gauss_xdiff","gauss",*m,*mean,*sigma);
	//Gauss bkg====
	RooRealVar *mean_bkg=new RooRealVar("mean_xdiff_bkg","Mean of bkg Gaussian",0.0,-10.0,10.0);
	RooRealVar *sigma_bkg=new RooRealVar("sigma_xdiff_bkg","Width of bkg Gaussian",-10.0,-1000.0,0.0);
	RooExponential *exp_bkg=new RooExponential("exp_xdiff_bkg","exp bkg",*m,*sigma_bkg);
	
// 	RooChebychev *exp_bkg=new RooChebychev("exp_xdiff_bkg","exp bkg",*m,RooArgList(*mean_bkg));
	
	RooRealVar *nsig=new RooRealVar("nS_xdiff","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_xdiff","background fraction",500,0,100000);
	RooRealVar *soverb=new RooRealVar("SoverB_xdiff","background fraction",0.7,0.4,1.0);
	RooAddPdf *model=new RooAddPdf("model_xdiff","model_xdiff",*gauss,*exp_bkg,*soverb);
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_xdiff","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_xdiff","rooFit_xdiff");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
// 	model->paramOn(mframe);
// 	model->plotOn(mframe,VisualizeError(*r,1),LineColor(kCyan));
	model->plotOn(mframe,Name("model_xdiff"),LineColor(kBlue));
	model->plotOn(mframe,Name("Signal Gaussian_xdiff"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Exponential background_xdiff"),Components(*exp_bkg),LineColor(kGreen));
	mframe->SetName("mframe_xdiff");
	mframe->Draw();
	mframe->SetXTitle("x_{diff}");
	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("Signal Gaussian_xdiff","Signal", "L");
	leg->AddEntry("Exponential background_xdiff","Background", "L");
	leg->Draw();
	can->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	model->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
// 	can->SaveAs("testxdiff.eps","eps");
// 	mean_bkg->Write();
// 	sigma_bkg->Write();
	exp_bkg->Write();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_xdiff_sig");
	TF1* f2 = exp_bkg->asTF( RooArgList(*m) );
	f2->Write("tf1_xdiff_bkg");
// 	return model;
	//      frac.Print();
}


void FitProp(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=0.1;
	double high=1;
	RooRealVar *m=new RooRealVar("prop","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_prop","Mean of Gaussian",0.45,0.35,0.7);
	RooRealVar *sigma=new RooRealVar("sigma_prop","Width of Gaussian",0.1,0.01,0.2);
	RooGaussian *gauss=new RooGaussian("gauss_prop","gauss",*m,*mean,*sigma);
	
// 	RooRealVar *mean2=new RooRealVar("mean2_prop","Mean of 2nd Gaussian",0.45,0.3,0.7);
// 	RooRealVar *sigma2=new RooRealVar("sigma2_prop","Width of 2nd Gaussian",0.1,0.01,0.2);
// 	RooGaussian *gauss2=new RooGaussian("gauss2_prop","gauss2",*m,*mean2,*sigma2);
	//pol2bkg====
// 	RooRealVar pol1("p1","p1",0,-100,100);
// 	RooRealVar pol2("p2","p2",0,-100,100);
// 	RooRealVar pol3("p3","p3",0,-100,100);
// 	RooPolynomial pol("pol","pol",m,RooArgList(pol1,pol2,pol3));
	RooRealVar *mean_2=new RooRealVar("mean_prop_2","Mean of Gaussian",0.2,0.2,0.35);
	RooRealVar *sigma_2=new RooRealVar("sigma_prop_2","Width of Gaussian",0.2,0.01,0.2);
	RooGaussian *gauss_2=new RooGaussian("gauss_prop_2","gauss_2",*m,*mean_2,*sigma_2);
	
	RooRealVar *mean_landau=new RooRealVar("mean_landau_prop","Mean of Landau",0.2,0.01,0.6);
	RooRealVar *sigma_landau=new RooRealVar("sigma_landau_prop","Width of Landau",0.2,0.01,0.5);
	RooLandau *landau=new RooLandau("landau_prop","landau",*m,*mean_landau,*sigma_landau);
	RooRealVar *ng=new RooRealVar("ng_1_prop","signal 1 fraction",200,0,50000);
	RooRealVar *nl=new RooRealVar("nl_2_prop","signal 2 fraction",500,0,50000);
	RooRealVar *soverbtmp=new RooRealVar("SoverB_proptmp","background fraction",0.7,0.4,1.0);
	RooAddPdf *bkg=new RooAddPdf("prop_sig","prop_sig",*landau,*gauss_2,*soverbtmp);
	
	
	RooRealVar *nsig=new RooRealVar("nS_prop","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_prop","background fraction",500,0,100000);
	RooRealVar *soverb=new RooRealVar("SoverB_prop","background fraction",0.7,0.4,1.0);
	RooAddPdf *model=new RooAddPdf("model_prop","model_prop",*gauss,*bkg,*soverb);
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_prop","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_prop","rooFit_prop");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
// 	model->plotOn(mframe,VisualizeError(*r,1),LineColor(kRed));
	model->plotOn(mframe,Name("model_prop"),LineColor(kBlue));
// 	model->paramOn(mframe);
	model->plotOn(mframe,Name("Signal Gaussian"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Background Landau+Gaussian"),Components(*bkg),LineColor(kGreen));
	model->plotOn(mframe,Name("Background Gaussian"),Components(*gauss_2),LineColor(kGreen),LineWidth(1));
	model->plotOn(mframe,Name("Landau background"),Components(*landau),LineColor(kGreen),LineWidth(1));
	data->plotOn(mframe);
	mframe->SetName("mframe_prop");
	mframe->Draw();
	mframe->SetXTitle("proportion");
	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("Signal Gaussian","Signal", "L");
	leg->AddEntry("Background Landau+Gaussian","Background", "L");
	leg->Draw();
	can->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	model->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	gauss_2->Write();
	mean_landau->Write();
	sigma_landau->Write();
	landau->Write();
	bkg->Write();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_prop_sig");
	TF1* f2 = bkg->asTF( RooArgList(*m) );
	f2->Write("tf1_prop_bkg");
// 	can->SaveAs("testprop.eps","eps");
}


void FitSigma(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=0;
	double high=5;
	RooRealVar *m=new RooRealVar("sigma","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_sigma","Mean of Gaussian",1.2,0.5,1.7);
	RooRealVar *sigma=new RooRealVar("sigma_sigma","Width of Gaussian",0.5,0.01,1);
	RooGaussian *gauss=new RooGaussian("gauss_sigma","gauss",*m,*mean,*sigma);
	//pol2bkg====
// 	RooRealVar pol1("p1","p1",0,-100,100);
// 	RooRealVar pol2("p2","p2",0,-100,100);
// 	RooRealVar pol3("p3","p3",0,-100,100);
// 	RooPolynomial pol("pol","pol",m,RooArgList(pol1,pol2,pol3));
	RooRealVar *mean_landau=new RooRealVar("mean_landau_sigma","Mean of Landau",1.5,0.01,3.0);
	RooRealVar *sigma_landau=new RooRealVar("sigma_landau_sigma","Width of Landau",0.5,0.01,0.1);
	RooLandau *landau=new RooLandau("landau_sigma","landau",*m,*mean_landau,*sigma_landau);
	
	
	RooRealVar *nsig=new RooRealVar("nS_sigma","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_sigma","background fraction",500,0,100000);
	RooAddPdf *model=new RooAddPdf("model_sigma","model_sigma",RooArgList(*gauss,*landau),RooArgList(*nsig,*nbkg));
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_sigma","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_sigma","rooFit_sigma");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
	model->plotOn(mframe,VisualizeError(*r,1),Name("model_sigma"),LineColor(kCyan));
	model->plotOn(mframe,Name("model_sigma"),LineColor(kBlue));
	model->paramOn(mframe);
	model->plotOn(mframe,Name("Signal Gaussian"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Landau background"),Components(*landau),LineColor(kGreen));
	mframe->SetName("mframe_sigma");
	mframe->Draw();
	mframe->SetXTitle("Peak sigma");
	can->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	model->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	mean_landau->Write();
	sigma_landau->Write();
	landau->Write();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_sigma_sig");
	TF1* f2 = landau->asTF( RooArgList(*m) );
	f2->Write("tf1_sigma_bkg");
}

void FitChi2(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=0;
	double high=10;
	RooRealVar *m=new RooRealVar("chi2","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_chi2","Mean of Gaussian",3,1,05);
	RooRealVar *sigma=new RooRealVar("sigma_chi2","Width of Gaussian",1,0.1,2);
	RooGaussian *gauss=new RooGaussian("gauss_chi2","gauss",*m,*mean,*sigma);
	//pol2bkg====
// 	RooRealVar pol1("p1","p1",0,-100,100);
// 	RooRealVar pol2("p2","p2",0,-100,100);
// 	RooRealVar pol3("p3","p3",0,-100,100);
// 	RooPolynomial pol("pol","pol",m,RooArgList(pol1,pol2,pol3));
	RooRealVar *mean_landau=new RooRealVar("mean_landau_chi2","Mean of Landau",1.5,0.5,4);
	RooRealVar *sigma_landau=new RooRealVar("sigma_landau_chi2","Width of Landau",0.1,0.01,3);
	RooLandau *landau=new RooLandau("landau_chi2","landau",*m,*mean_landau,*sigma_landau);
	
	
	RooRealVar *nsig=new RooRealVar("nS_chi2","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_chi2","background fraction",500,0,100000);
	RooAddPdf *model=new RooAddPdf("model_chi2","model_chi2",RooArgList(*gauss,*landau),RooArgList(*nsig,*nbkg));
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_chi2","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_chi2","rooFit_chi2");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
	model->plotOn(mframe,VisualizeError(*r,1),LineColor(kRed));
	model->plotOn(mframe,Name("model_chi2"),LineColor(kBlue));
	model->paramOn(mframe);
	model->plotOn(mframe,Name("Signal Gaussian"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Landau background"),Components(*landau),LineColor(kGreen));
	mframe->SetName("mframe_chi2");
	mframe->Draw();
	mframe->SetXTitle("Peak chi2");
	can->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	m->Write();
	model->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	mean_landau->Write();
	sigma_landau->Write();
	landau->Write();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_chi2_sig");
	TF1* f2 = landau->asTF( RooArgList(*m) );
	f2->Write("tf1_chi2_bkg");
}

void FitNPEs(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=200;
	double high=5000;
	RooRealVar *m=new RooRealVar("npes","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_npes","Mean of Gaussian",2500.0,500.0,6000.0);
	RooRealVar *sigma=new RooRealVar("sigma_npes","Width of Gaussian",1000.0,600.0,1500.0);
	RooGaussian *gauss=new RooGaussian("gauss_npes","gauss",*m,*mean,*sigma);
	//pol2bkg====
// 	RooRealVar pol1("p1","p1",-10,-1000,10);
// 	RooRealVar pol2("p2","p2",-1,-1000,10);
// 	RooRealVar pol3("p3","p3",0,-10,10);
// 	RooPolynomial pol("pol","pol",m,RooArgList(pol1,pol2,pol3));
	
	RooRealVar *mean_landau=new RooRealVar("mean_landau_npes","Mean of Landau",500.0,0.0,5000.0);
	RooRealVar *sigma_landau=new RooRealVar("sigma_landau_npes","Width of Landau",1000.0,100.0,5000.0);
	RooLandau *landau=new RooLandau("landau_npes","landau",*m,*mean_landau,*sigma_landau);
// 	RooPolynomial *landau=new RooPolynomial("landau_npes","landau",*m,RooArgList(*mean_landau,*sigma_landau));
	
	RooRealVar *nsig=new RooRealVar("nS_npes","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_npes","background fraction",500,0,100000);
	RooRealVar *soverb=new RooRealVar("SoverB_npes","background fraction",0.7,0.4,1.0);
	RooAddPdf *model=new RooAddPdf("model_npes","model_npes",*gauss,*landau,*soverb);
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_npes","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_npes","rooFit_npes");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
// 	model->plotOn(mframe,VisualizeError(*r,1),LineColor(kRed));
	model->plotOn(mframe,Name("model_npes"),LineColor(kBlue));
// 	model->paramOn(mframe);
	model->plotOn(mframe,Name("Signal Gaussian"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Landau background"),Components(*landau),LineColor(kGreen));
	data->plotOn(mframe);
	mframe->SetName("mframe_npes");
	mframe->Draw();
	mframe->SetXTitle("N_{hits}");
	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("Signal Gaussian","Signal", "L");
	leg->AddEntry("Landau background","Background", "L");
	leg->Draw();
	can->Write();
	model->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	mean_landau->Write();
	sigma_landau->Write();
	landau->Write();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_npes_sig");
	TF1* f2 = landau->asTF( RooArgList(*m) );
	f2->Write("tf1_npes_bkg");
// 	can->SaveAs("testnpes.eps","eps");
// 	return model;
	//      frac.Print();
}

double FitInitDist(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=1;
	double high=2.5;
	RooRealVar *m=new RooRealVar("dist","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_dist","Mean of Gaussian",1.8,1.5,5);
	RooRealVar *sigma=new RooRealVar("sigma_dist","Width of Gaussian",0.3,0.001,0.5);
	RooGaussian *gauss=new RooGaussian("gauss_dist","gauss",*m,*mean,*sigma);
	//pol2bkg====
	RooRealVar *pol1=new RooRealVar("p1_dist","p1",0,-2000,2000);
// 	RooRealVar pol2("p2_dist","p2",0,-1000,1000);
// 	RooRealVar pol3("p3_dist","p3",200,-100,100);
// 	RooRealVar pol4("p4_dist","p4",200,-100,100);
// 	RooRealVar pol5("p5_dist","p5",200,-100,100);
	RooPolynomial *pol=new RooPolynomial("pol_dist","pol",*m,RooArgList(*pol1));
	
	RooRealVar *nsig=new RooRealVar("nS_dist","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_dist","background fraction",500,0,100000);
	RooAddPdf *model=new RooAddPdf("model_dist","model_dist",RooArgList(*gauss,*pol),RooArgList(*nsig,*nbkg));
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_dist","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_init_dist","rooFit_init_dist");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
	model->plotOn(mframe,VisualizeError(*r,1),LineColor(kCyan));
	model->plotOn(mframe,Name("model_dist"),LineColor(kBlue));
	model->paramOn(mframe);
	model->plotOn(mframe,Name("Signal Gaussian"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Polinomial background"),Components(*pol),LineColor(kGreen));
	mframe->SetName("mframe_init_dist");
	mframe->Draw();
	mframe->SetXTitle("Peak Dist Init");
	can->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	model->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	pol1->Write();
	pol->Write();
	cout << "mean dist"<<endl;
	mean->getVal();
	cout << "sigma dist"<<endl;
	sigma->getVal();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_dist_sig");
	TF1* f2 = pol->asTF( RooArgList(*m) );
	f2->Write("tf1_dist_bkg");
	return mean->getVal();
}

void FitDistVar(TTree* t1){
	//RooFit=============================================================== 
	
	//Configure
	double low=4;
	double high=18;
	RooRealVar *m=new RooRealVar("distvar","m",low,high);
	//GaussianSignal======
	RooRealVar *mean=new RooRealVar("mean_distvar","Mean of Gaussian",11,10,12);
	RooRealVar *sigma=new RooRealVar("sigma_distvar","Width of Gaussian",2,0.5,3);
	RooGaussian *gauss=new RooGaussian("gauss_distvar","gauss",*m,*mean,*sigma);
	
// 	RooRealVar *exp_weight=new RooRealVar("exp_weight","g1 fraction",500,0,100000);
// 	RooRealVar *expo_c_bkg=new RooRealVar("c_expo","exponential constant",-0.2,-0.05,-0.3);
// 	RooExponential *exp_bkg=new RooExponential("exp_distvar_bkg","expo bkg",*m,*expo_c_bkg);
	
	RooRealVar *g1_weight=new RooRealVar("g1_weight","g1 fraction",100,0,100000);
	RooRealVar *mean_bkg1=new RooRealVar("mean_distvar_bkg_1","Mean of Gaussian",6.0,5,7);
	RooRealVar *sigma_bkg1=new RooRealVar("sigma_distvar_bkg_1","Width of Gaussian",1,0.5,2);
	RooGaussian *gauss_bkg1=new RooGaussian("gauss_distvar_bkg_1","gauss",*m,*mean_bkg1,*sigma_bkg1);
	
	RooRealVar *g2_weight=new RooRealVar("g2_weight","g2 fraction",500,0,100000);
	RooRealVar *mean_bkg2=new RooRealVar("mean_distvar_bkg_2","Mean of Gaussian",8.5,7.2,9);
	RooRealVar *sigma_bkg2=new RooRealVar("sigma_distvar_bkg_2","Width of Gaussian",0.5,0.1,1.5);
	RooGaussian *gauss_bkg2=new RooGaussian("gauss_distvar_bkg_2","gauss",*m,*mean_bkg2,*sigma_bkg2);
	
	RooRealVar *g3_weight=new RooRealVar("g3_weight","g3 fraction",500,0,100000);
	RooRealVar *mean_bkg3=new RooRealVar("mean_distvar_bkg_3","Mean of Gaussian",11,10,12);
	RooRealVar *sigma_bkg3=new RooRealVar("sigma_distvar_bkg_3","Width of Gaussian",0.7,0.1,1.5);
	RooGaussian *gauss_bkg3=new RooGaussian("gauss_distvar_bkg_3","gauss",*m,*mean_bkg3,*sigma_bkg3);
	
	RooRealVar *g4_weight=new RooRealVar("g4_weight","g4 fraction",500,0,100000);
	RooRealVar *mean_bkg4=new RooRealVar("mean_distvar_bkg_4","Mean of Gaussian",13.5,12,15);
	RooRealVar *sigma_bkg4=new RooRealVar("sigma_distvar_bkg_4","Width of Gaussian",0.8,0.1,2);
	RooGaussian *gauss_bkg4=new RooGaussian("gauss_distvar_bkg_4","gauss",*m,*mean_bkg4,*sigma_bkg4);
	
	RooAddPdf *fourgauss_peaks=new RooAddPdf("fourgauss","fourgauss",RooArgList(*gauss_bkg1,*gauss_bkg2,*gauss_bkg3,*gauss_bkg4),RooArgList(*g1_weight,*g2_weight,*g3_weight,*g4_weight));
	
	RooRealVar *nsig=new RooRealVar("nS_distvar","signal fraction",500,0,100000);
	RooRealVar *nbkg=new RooRealVar("nB_distvar","background fraction",500,0,100000);
	RooAddPdf *model=new RooAddPdf("model_distvar","model_distvar",RooArgList(*gauss,*fourgauss_peaks),RooArgList(*nsig,*nbkg));
	
	//Converth1toRooFitdataclass:
	RooDataSet *data=new RooDataSet("data_distvar","data",t1,*m);
	RooFitResult *r=model->fitTo(*data,Save());
	
	TCanvas *can=new TCanvas("rooFit_distvar","rooFit_distvar");
	RooPlot* mframe=m->frame();
	data->plotOn(mframe);
	model->plotOn(mframe,VisualizeError(*r,1),FillColor(kCyan));
	model->plotOn(mframe,Name("model_distvar"),LineColor(kBlue));
	model->paramOn(mframe);
	model->plotOn(mframe,Name("Signal Gaussian"),Components(*gauss),LineColor(kRed));
	model->plotOn(mframe,Name("Polinomial background"),Components(*fourgauss_peaks),LineColor(kGreen));
	mframe->SetName("mframe_distvar");
	mframe->Draw();
	mframe->SetXTitle("Peak DistVar");
	can->Write();
// 	r->Write();
	nsig->Write();
	nbkg->Write();
	data->Write();
	model->Write();
	m->Write();
	mean->Write();
	sigma->Write();
	gauss->Write();
	fourgauss_peaks->Write();
// 	pol->Write();
// 	cout << "mean distvar"<<endl;
// 	mean->getVal();
// 	cout << "sigma distvar"<<endl;
// 	sigma->getVal();
	TF1* f = gauss->asTF( RooArgList(*m) );
	f->Write("tf1_distvar_sig");
	TF1* f2 = fourgauss_peaks->asTF( RooArgList(*m) );
	f2->Write("tf1_distvar_bkg");
}

TH2I* GetHPDMapIFB(TList *hpds, double *cutoff, TString box_string="U"){
	int thebox=-1;
	double npeaks_mean=0;
	double npeaks_rms=0;
	int count_peaks=0;
	double histo_overflow=2000;
	
	if(box_string=="U") thebox=0;
	else thebox=1;
	
	TH2I *h2=new TH2I(Form("badhpds2D_%s",box_string.Data()),Form("Number of reconstructed peaks %s Box;HPD;Column",box_string.Data()),nHPDs,0,nHPDs,nCOLs,0,nCOLs);
	TH1I *h1=new TH1I(Form("badhpds1D_%s",box_string.Data()),Form("Number of reconstructed peaks %s Box;Number of peaks;HPDs", box_string.Data()), histo_overflow,0, histo_overflow);
	
	
	
	
	int i=0;
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		if(hpd->GetNPeaks()>0){
			h1->Fill(hpd->GetNPeaks());
			npeaks_mean+=(double)hpd->GetNPeaks();
			count_peaks++;
			h2->Fill(hpd->GetHPD(),hpd->GetColumn(),hpd->GetNPeaks());
		}
	}
	
	npeaks_mean/=count_peaks;
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd2=(CleanHPD*)hpds->At(i);
		if(hpd2->GetNPeaks()>0){
			npeaks_rms+=((double)(hpd2->GetNPeaks())-npeaks_mean)*((double)(hpd2->GetNPeaks())-npeaks_mean);
		}
	}
	
	h2->SetMinimum(0);
// 	h2->SetMaximum(histo_overflow);
	npeaks_rms/=count_peaks;
	npeaks_rms=sqrt(npeaks_rms);
	
// 	h1->Write();
	
	cutoff[0]=npeaks_mean;
	cutoff[1]=npeaks_rms;
	cutoff[2]=count_peaks;
	cout << "mean npeaks: "<<cutoff[0]<<endl;
	cout << "rms npeaks: "<<cutoff[1]<<endl;
	cout << "good npeaks: "<<cutoff[2]<<endl;
	return h2;
}

void FlagIFBHPDs(TList *hpds, double *cutoff, double sigmas){
	int i=0;
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		if(hpd->GetNPeaks()>(cutoff[0]+sigmas*cutoff[1])){
			hpd->SetPassLevel("IFB");
			cout << "HPD " << hpd->GetHPD() << " on column " << hpd->GetColumn() << " has been flagged as IFB affected within "<<sigmas<<" sigmas."<<endl;
			TGraph *g=hpd->DrawPeakMap();
			g->SetName(Form("IFB_HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()));
			g->Write();
			TCanvas *tmpcan=new TCanvas(Form("can_IFB_Histo_HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()),Form("can_IFB_Histo_HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()),800,800);
			TH2D *h2=hpd->DrawPeakHisto();
			h2->SetName(Form("IFB_Histo_HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()));
// 			h2->SetMaximum(200);
			h2->SetMinimum(0);
			h2->Draw("colz");
			h2->Write();
			tmpcan->Write();
			tmpcan->Close();
			delete tmpcan;
		}
		else{
			TGraph *g=hpd->DrawPeakMap();
			g->SetName(Form("HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()));
			g->Write();
			TCanvas *tmpcan=new TCanvas(Form("can_Histo_HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()),Form("can_IFB_Histo_HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()),800,800);
			TH2D *h2=hpd->DrawPeakHisto();
			h2->SetName(Form("Histo_HPD%i_COL%i_BOX%i",hpd->GetHPD(),hpd->GetColumn(),hpd->GetBox()));
// 			h2->SetMaximum(200);
			h2->SetMinimum(0);
			h2->Draw("colz");
			h2->Write();
			tmpcan->Write();
			tmpcan->Close();
			delete tmpcan;
		}
	}
}

void SetDistVar(TList * List){
	int i=0,k=0;
	int count_peaks=0;
	
	for(i=0;i<List->GetSize();i++){
		CleanPeak *peak=(CleanPeak*)List->At(i);
		count_peaks=0;
		peak->SetDistVar(0);
		double tmpdist=0;
		for(k=0;k<List->GetSize();k++){
			CleanPeak *peak2=(CleanPeak*)List->At(k);
			if(abs(peak->GetCalibStep()-peak2->GetCalibStep())==1 || 1){
				tmpdist+=TMath::Sqrt( (peak->GetMuX()-peak2->GetMuX())*(peak->GetMuX()-peak2->GetMuX()) + (peak->GetMuY()-peak2->GetMuY())*(peak->GetMuY()-peak2->GetMuY()) );
				count_peaks++;
			}
		}
		if(count_peaks!=0) peak->SetDistVar(tmpdist/count_peaks);
// 		cout << "get dist var "<< peak->GetHPD() << "\t" << peak->GetColumn() <<endl;
// 		cout << "get dist var "<< peak->GetDistVar() <<endl;
	}
}

void FillDistVar(TList *hpds){
	int i=0,k=0;
	int tmpstep=-1;
	cout << "Filling dist var..." << endl;
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		TList *peaklist=(TList*)hpd->GetHitmap();
		tmpstep=-1;
// 		peakcounter=0;
// 		currentpeak=0;
		TList *MiniList=new TList();
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			MiniList->Add(peak);
		}
		SetDistVar(MiniList);
// 		MiniList->Clear("nodelete");
// 		peaklist->Clear("nodelete");
	}
	cout << "DONE" << endl;
}

double GetXmean_atStep(TList *hpds,int step){
	int i=0,k=0;
	double meanX=0;
	int count=0;
	bool mmh=false;
	for(i=0;i<hpds->GetSize();i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		if(hpd->GetPassLevel()!="IFB"){
			
			TList *peaklist=(TList*)hpd->GetHitmap();
			
			for(k=0;k<peaklist->GetSize();k++){
				CleanPeak *peak=(CleanPeak*)peaklist->At(k);
				if(peak->GetCalibStep()==step && peak->IsSingle() ){
					mmh=true;
					count++;
					meanX+=peak->GetMuX();  //yes, that's right, X on plots is actually Y in panel coordinates
				}
			}
// 			peaklist->Clear("nodelete");
		}
	}
	if(count>0) return meanX/(double)count;
	
// 	if(mmh) cout << meanX << "\t" << count << endl;
	
	return 9999;
}

void SetXDiff(CleanHPD* hpd,double meanX,int step){
	int i=0;
	TList *peaklist=(TList*)hpd->GetHitmap();
	for(i=0;i<peaklist->GetSize();i++){
		CleanPeak *peak=(CleanPeak*)peaklist->At(i);
		if(peak->GetCalibStep()==step) peak->SetXDiff(fabs(peak->GetMuX()-meanX)); //yes, that's right, X on plots is actually Y in panel coordinates
	}
// 	peaklist->Clear("nodelete");
}

void FillXDiff(TList *hpds){
	int i=0,k=0,j=0;
	int tmpstep=-1;
	cout << "Filling X diff..." << endl;
	for(k=0;k<nCOLs;k++){
		TList *MiniList=new TList();
		for(i=0; i<nCOLs*nHPDs; i++){
			CleanHPD *hpd=(CleanHPD*)hpds->At(i);
			if(hpd->GetColumn()==k){
				MiniList->Add(hpd);
			}
		}
		cout << "Column "<<k<<endl;
		for(i=0;i<nSTEPs;i++){
			double Xmean=GetXmean_atStep(MiniList,i);
// 			cout << "MeanX: "<<Xmean<<endl;
			for(j=0;j<hpds->GetSize();j++){
				CleanHPD* hpd=(CleanHPD*)hpds->At(j);
				if(hpd->GetColumn()==k){
					SetXDiff(hpd,Xmean,i);
				}
			}
		}
// 		MiniList->Clear("nodelete");
		
	}
// 	cout << "DONE" << endl;
}

void CreateNVDistribution(TList *hpds, double meandistance=0){
	int i=0,k=0,l=0;
	
	double max_newdistvar=20, min_newdistvar=0;
	//good peaks
	TH1D *h_newdistvar_good=new TH1D("h_newdistvar_good","Good Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	h_newdistvar_good->SetLineColor(4);
	
	//badpeaks
	TH1D *h_newdistvar_bad=new TH1D("h_newdistvar_bad","Bad Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	h_newdistvar_bad->SetLineColor(2);
	
	//totpeaks
	TH1D *h_newdistvar_tot=new TH1D("h_newdistvar_tot","Tot Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		TList *peaklist=(TList*)hpd->GetHitmap();
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			double peakradius=TMath::Sqrt((peak->GetMuX()-hpd->GetMeanX())*(peak->GetMuX()-hpd->GetMeanX())+(peak->GetMuY()-hpd->GetMeanY())*(peak->GetMuY()-hpd->GetMeanY()));
// 			cout << peak->GetMuX() << "\t" << peak->GetMeanX() << endl;
			//get rid of fake peaks
			bool flaggoodpeak=true;
			if(peak->GetMax()<10 || peak->GetSigma()<0.5 || peak->GetChi2()<0.5){
				peak->SetPassLevel("DIS_STEP1");
				flaggoodpeak=false;
			}
			
			//fill histos
			
			if(peak->IsSingle() && hpd->GetPassLevel()!="IFB" && flaggoodpeak){
				h_newdistvar_good->Fill(peak->GetDistVar()/meandistance);
			}
			else {
				if(flaggoodpeak){
					h_newdistvar_bad->Fill(peak->GetDistVar()/meandistance);
				}
			}
		}
// 		peaklist->Clear("nodelete");
	}
	h_newdistvar_tot->Add(h_newdistvar_good,h_newdistvar_bad);
	
	TCanvas *can_newdistvar=new TCanvas("can_newdistvar","can_newdistvar");
	h_newdistvar_tot->Draw();
	h_newdistvar_good->Draw("sames");
	h_newdistvar_bad->Draw("sames");
	can_newdistvar->Write();
	
}

void CreateDistributions_SingleHPDs(CleanHPD *hpd){
	int i=0,k=0,l=0;
	
	int ncol=hpd->GetColumn();
	int nhpd=hpd->GetHPD();
	
	double max_max=140, min_max=0;
	double max_npes=5000, min_npes=0;
	double max_sigma=5, min_sigma=0;
	double max_prop=1, min_prop=0;
	double max_rot=0.2, min_rot=-0.2;
	double max_rms=2, min_rms=0;
	double max_chi2=10, min_chi2=0;
	double max_distvar=20, min_distvar=0;
	double max_newdistvar=20, min_newdistvar=0;
	double max_theta=3.2, min_theta=-3.2;
	double max_rmsratio=5, min_rmsratio=0;
	double max_radius=20, min_radius=0;
	double max_rmsdistortion=5, min_rmsdistortion=-5;
	double max_xdiff=20, min_xdiff=0;
	//good peaks
	TH1D *h_max_good=new TH1D(Form("h_max_good_hpd%i_col%i",nhpd,ncol),"Good Peaks Max;max;",100,min_max,max_max);
	h_max_good->SetLineColor(4);
	h_max_good->SetLineWidth(2);
	TH1D *h_nPE_good=new TH1D(Form("h_nPE_good_hpd%i_col%i",nhpd,ncol),"Good Peaks nPE;nPE;",100,min_npes,max_npes);
	h_nPE_good->SetLineColor(4);
	h_nPE_good->SetLineWidth(2);
	TH1D *h_sigma_good=new TH1D(Form("h_sigma_good_hpd%i_col%i",nhpd,ncol),"Good Peaks #sigma ;#sigma (LHCb pixels);",100,min_sigma,max_sigma);
	h_sigma_good->SetLineColor(4);
	h_sigma_good->SetLineWidth(2);
	TH1D *h_xdiff_good=new TH1D(Form("h_xdiff_good_hpd%i_col%i",nhpd,ncol),"Good Peaks #xdiff ;#xdiff (LHCb pixels);",100,min_xdiff,max_xdiff);
	h_xdiff_good->SetLineColor(4);
	h_xdiff_good->SetLineWidth(2);
	TH1D *h_prop_good=new TH1D(Form("h_prop_good_hpd%i_col%i",nhpd,ncol),"Good Peaks Proportion; proportion;",100,min_prop,max_prop);
	h_prop_good->SetLineColor(4);
	h_prop_good->SetLineWidth(2);
	TH1D *h_rot_good=new TH1D(Form("h_rot_good_hpd%i_col%i",nhpd,ncol),"Good Peaks Rotation; rotation;",800,min_rot,max_rot);
	h_rot_good->SetLineColor(4);
	h_rot_good->SetLineWidth(2);
	TH1D *h_rms_good_X=new TH1D(Form("h_rms_good_X_hpd%i_col%i",nhpd,ncol),"Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	h_rms_good_X->SetLineColor(4);
	h_rms_good_X->SetLineWidth(2);
	TH1D *h_rms_good_Y_=new TH1D(Form("h_rms_good_Y__hpd%i_col%i",nhpd,ncol),"Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	h_rms_good_Y_->SetLineColor(4);
	h_rms_good_Y_->SetLineWidth(2);
	TH1D *h_chi2_good=new TH1D(Form("h_chi2_good_hpd%i_col%i",nhpd,ncol),"Good Peaks #chi^{2}; #chi^{2};",100,min_chi2,max_chi2);
	h_chi2_good->SetLineColor(4);
	h_chi2_good->SetLineWidth(2);
	TH1D *h_distvar_good=new TH1D(Form("h_distvar_good_hpd%i_col%i",nhpd,ncol),"Good Peaks dist var; dist var;",400,min_distvar,max_distvar);
	h_distvar_good->SetLineColor(4);
	h_distvar_good->SetLineWidth(2);
	TH1D *h_theta_good=new TH1D(Form("h_theta_good_hpd%i_col%i",nhpd,ncol),"Good Peaks #theta; #theta (rad);",100,min_theta,max_theta);
	h_theta_good->SetLineColor(4);
	h_theta_good->SetLineWidth(2);
	TH1D *h_rmsratio_good=new TH1D(Form("h_rmsratio_good_hpd%i_col%i",nhpd,ncol),"Good Peaks rmsratio; rmsy/rmsx;",200,min_rmsratio,max_rmsratio);
	h_rmsratio_good->SetLineColor(4);
	h_rmsratio_good->SetLineWidth(2);
	TH1D *h_radius_good=new TH1D(Form("h_radius_good_hpd%i_col%i",nhpd,ncol),"Good Peaks radius; radius (LHCb Pixels);",200,min_radius,max_radius);
	h_radius_good->SetLineColor(4);
	h_radius_good->SetLineWidth(2);
	TH1D *h_rmsdistortion_good=new TH1D(Form("h_rmsdistortion_good_hpd%i_col%i",nhpd,ncol),"Good Peaks rms distortion; rms distortion;",200,min_rmsdistortion,max_rmsdistortion);
	h_rmsdistortion_good->SetLineColor(4);
	h_rmsdistortion_good->SetLineWidth(2);
	TH1D *h_newdistvar_good=new TH1D(Form("h_newdistvar_good_hpd%i_col%i",nhpd,ncol),"Good Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	h_newdistvar_good->SetLineColor(4);
	h_newdistvar_good->SetLineWidth(2);
	TH2D *h_maxprop_good=new TH2D(Form("h_maxprop_good_hpd%i_col%i",nhpd,ncol),"Good peaks maxprop",25,min_max,max_max,25,min_prop,max_prop);
	h_maxprop_good->SetMarkerColor(4);
	h_maxprop_good->SetLineColor(4);
// 	h_maxprop_good->SetLineWidth(2);
	TH2D *h_maxnpes_good=new TH2D(Form("h_maxnpes_good_hpd%i_col%i",nhpd,ncol),"Good peaks maxnpes",25,min_max,max_max,25,min_npes,max_npes);
	h_maxnpes_good->SetMarkerColor(4);
	h_maxnpes_good->SetLineColor(4);
	TH2D *h_npesprop_good=new TH2D(Form("h_npesprop_good_hpd%i_col%i",nhpd,ncol),"Good peaks npesprop",25,min_npes,max_npes,25,min_prop,max_prop);
	h_npesprop_good->SetMarkerColor(4);
	h_npesprop_good->SetLineColor(4);
	TH2D *h_radiustheta_good=new TH2D(Form("h_radiustheta_good_hpd%i_col%i",nhpd,ncol),"Good peaks radiustheta",25,min_radius,max_radius,25,min_theta,max_theta);
	h_radiustheta_good->SetMarkerColor(4);
	h_radiustheta_good->SetLineColor(4);
	
	//badpeaks
	TH1D *h_max_bad=new TH1D(Form("h_max_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks Max;max;",100,min_max,max_max);
	h_max_bad->SetLineColor(2);
	h_max_bad->SetLineWidth(2);
	TH1D *h_nPE_bad=new TH1D(Form("h_nPE_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks nPE;nPE;",100,min_npes,max_npes);
	h_nPE_bad->SetLineColor(2);
	h_nPE_bad->SetLineWidth(2);
	TH1D *h_sigma_bad=new TH1D(Form("h_sigma_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks #sigma ;#sigma (LHCb pixels);",100,min_sigma,max_sigma);
	h_sigma_bad->SetLineColor(2);
	h_sigma_bad->SetLineWidth(2);
	TH1D *h_xdiff_bad=new TH1D(Form("h_xdiff_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks #xdiff ;#xdiff (LHCb pixels);",100,min_xdiff,max_xdiff);
	h_xdiff_bad->SetLineColor(2);
	h_xdiff_bad->SetLineWidth(2);
	TH1D *h_prop_bad=new TH1D(Form("h_prop_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks Proportion; proportion;",100,min_prop,max_prop);
	h_prop_bad->SetLineColor(2);
	h_prop_bad->SetLineWidth(2);
	TH1D *h_rot_bad=new TH1D(Form("h_rot_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks Rotation; rotation;",800,min_rot,max_rot);
	h_rot_bad->SetLineColor(2);
	h_rot_bad->SetLineWidth(2);
	TH1D *h_rms_bad_X=new TH1D(Form("h_rms_bad_X_hpd%i_col%i",nhpd,ncol),"Good Peaks RMS; RMS x (LHCb pixels);",100,min_rms,max_rms);
	h_rms_bad_X->SetLineColor(2);
	h_rms_bad_X->SetLineWidth(2);
	TH1D *h_rms_bad_Y=new TH1D(Form("h_rms_bad_Y_hpd%i_col%i",nhpd,ncol),"Good Peaks RMS; RMS y (LHCb pixels);",100,min_rms,max_rms);
	h_rms_bad_Y->SetLineColor(2);
	h_rms_bad_Y->SetLineWidth(2);
	TH1D *h_chi2_bad=new TH1D(Form("h_chi2_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks #chi^{2}; #chi^{2};",100,min_chi2,max_chi2);
	h_chi2_bad->SetLineColor(2);
	h_chi2_bad->SetLineWidth(2);
	TH1D *h_distvar_bad=new TH1D(Form("h_distvar_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks dist var; dist var;",400,min_distvar,max_distvar);
	h_distvar_bad->SetLineColor(2);
	h_distvar_bad->SetLineWidth(2);
	TH1D *h_theta_bad=new TH1D(Form("h_theta_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks #theta; #theta (rad);",100,min_theta,max_theta);
	h_theta_bad->SetLineColor(2);
	h_theta_bad->SetLineWidth(2);
	TH1D *h_rmsratio_bad=new TH1D(Form("h_rmsratio_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks rmsratio; rmsy/rmsx;",200,min_rmsratio,max_rmsratio);
	h_rmsratio_bad->SetLineColor(2);
	h_rmsratio_bad->SetLineWidth(2);
	TH1D *h_radius_bad=new TH1D(Form("h_radius_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks radius; radius (LHCb Pixels);",200,min_radius,max_radius);
	h_radius_bad->SetLineColor(2);
	h_radius_bad->SetLineWidth(2);
	TH1D *h_rmsdistortion_bad=new TH1D(Form("h_rmsdistortion_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks rms distortion; rms distortion;",200,min_rmsdistortion,max_rmsdistortion);
	h_rmsdistortion_bad->SetLineColor(2);
	h_rmsdistortion_bad->SetLineWidth(2);
	TH1D *h_newdistvar_bad=new TH1D(Form("h_newdistvar_bad_hpd%i_col%i",nhpd,ncol),"Bad Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	h_newdistvar_bad->SetLineColor(2);
	h_newdistvar_bad->SetLineWidth(2);
	TH2D *h_maxprop_bad=new TH2D(Form("h_maxprop_bad_hpd%i_col%i",nhpd,ncol),"Bad peaks maxprop",25,min_max,max_max,25,min_prop,max_prop);
	h_maxprop_bad->SetMarkerColor(2);
	h_maxprop_bad->SetLineColor(2);
	TH2D *h_maxnpes_bad=new TH2D(Form("h_maxnpes_bad_hpd%i_col%i",nhpd,ncol),"Bad peaks maxnpes",25,min_max,max_max,25,min_npes,max_npes);
	h_maxnpes_bad->SetMarkerColor(2);
	h_maxnpes_bad->SetLineColor(2);
	TH2D *h_npesprop_bad=new TH2D(Form("h_npesprop_bad_hpd%i_col%i",nhpd,ncol),"Bad peaks npesprop",25,min_npes,max_npes,25,min_prop,max_prop);
	h_npesprop_bad->SetMarkerColor(2);
	h_npesprop_bad->SetLineColor(2);
	TH2D *h_radiustheta_bad=new TH2D(Form("h_radiustheta_bad_hpd%i_col%i",nhpd,ncol),"Bad peaks radiustheta",25,min_radius,max_radius,25,min_theta,max_theta);
	h_radiustheta_bad->SetMarkerColor(2);
	h_radiustheta_bad->SetLineColor(2);
	
	//totpeaks
	TH1D *h_max_tot=new TH1D(Form("h_max_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks Max;M;",100,min_max,max_max);
	TH1D *h_nPE_tot=new TH1D(Form("h_nPE_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks nPE;N_{hits};",100,min_npes,max_npes);
	TH1D *h_sigma_tot=new TH1D(Form("h_sigma_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks #sigma ;#sigma (LHCb pixels);",100,min_sigma,max_sigma);
	TH1D *h_xdiff_tot=new TH1D(Form("h_xdiff_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks #xdiff ;x_{diff} (LHCb pixels);",100,min_xdiff,max_xdiff);
	TH1D *h_prop_tot=new TH1D(Form("h_prop_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks Proportion; proportion;",100,min_prop,max_prop);
	TH1D *h_rot_tot=new TH1D(Form("h_rot_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks Rotation; rotation;",800,min_rot,max_rot);
	TH1D *h_chi2_tot=new TH1D(Form("h_chi2_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks #chi^{2}; #chi^{2};",100,min_chi2,max_chi2);
	TH1D *h_distvar_tot=new TH1D(Form("h_distvar_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks dist var; dist var;",400,min_distvar,max_distvar);
	TH1D *h_theta_tot=new TH1D(Form("h_theta_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks #theta; #theta (rad);",100,min_theta,max_theta);
	TH1D *h_rmsratio_tot=new TH1D(Form("h_rmsratio_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks rmsratio; rmsy/rmsx;",200,min_rmsratio,max_rmsratio);
	TH1D *h_radius_tot=new TH1D(Form("h_radius_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks radius; radius (LHCb Pixels);",200,min_radius,max_radius);
	TH1D *h_rmsdistortion_tot=new TH1D(Form("h_rmsdistortion_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks rms distortion; rms distortion;",200,min_rmsdistortion,max_rmsdistortion);
	TH1D *h_newdistvar_tot=new TH1D(Form("h_newdistvar_tot_hpd%i_col%i",nhpd,ncol),"Tot Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	TH2D *h_maxprop_tot=new TH2D(Form("h_maxprop_tot_hpd%i_col%i",nhpd,ncol),"Tot peaks maxprop",25,min_max,max_max,25,min_prop,max_prop);
	TH2D *h_maxnpes_tot=new TH2D(Form("h_maxnpes_tot_hpd%i_col%i",nhpd,ncol),"Tot peaks maxnpes",25,min_max,max_max,25,min_npes,max_npes);
	TH2D *h_npesprop_tot=new TH2D(Form("h_npesprop_tot_hpd%i_col%i",nhpd,ncol),"Tot peaks npesprop",25,min_npes,max_npes,25,min_prop,max_prop);
	TH2D *h_radiustheta_tot=new TH2D(Form("h_radiustheta_tot_hpd%i_col%i",nhpd,ncol),"Tot peaks radiustheta",25,min_radius,max_radius,25,min_theta,max_theta);
	TH1D *h_rms_tot_X=new TH1D(Form("h_rms_tot_X_hpd%i_col%i",nhpd,ncol),"Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	TH1D *h_rms_tot_Y=new TH1D(Form("h_rms_tot_Y_hpd%i_col%i",nhpd,ncol),"Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	
	
	h_max_tot->SetLineWidth(2);
	h_nPE_tot->SetLineWidth(2);
	h_sigma_tot->SetLineWidth(2);
	h_xdiff_tot->SetLineWidth(2);
	h_prop_tot->SetLineWidth(2);
	h_rot_tot->SetLineWidth(2);
	h_chi2_tot->SetLineWidth(2);
	h_distvar_tot->SetLineWidth(2);
	h_theta_tot->SetLineWidth(2);
	h_rmsratio_tot->SetLineWidth(2);
	h_radius_tot->SetLineWidth(2);
	h_rmsdistortion_tot->SetLineWidth(2);
	h_newdistvar_tot->SetLineWidth(2);

	
		
	TList *peaklist=(TList*)hpd->GetHitmap();
	for(k=0; k<hpd->GetNPeaks(); k++){
		CleanPeak *peak=(CleanPeak*)peaklist->At(k);
		double peakradius=TMath::Sqrt((peak->GetMuX()-hpd->GetMeanX())*(peak->GetMuX()-hpd->GetMeanX())+(peak->GetMuY()-hpd->GetMeanY())*(peak->GetMuY()-hpd->GetMeanY()));
// 			cout << peak->GetMuX() << "\t" << peak->GetMeanX() << endl;
		//get rid of fake peaks
		bool flaggoodpeak=true;
		if(peak->GetMax()<10 || peak->GetSigma()<0.5 || peak->GetChi2()<0.5){
			peak->SetPassLevel("DIS_STEP1");
			flaggoodpeak=false;
		}
		
		//fill histos
		
		if(peak->IsSingle() && hpd->GetPassLevel()!="IFB" && flaggoodpeak){
			h_max_good->Fill(peak->GetMax());
			h_nPE_good->Fill(peak->GetNPEs());
			h_sigma_good->Fill(peak->GetSigma());
			h_xdiff_good->Fill(peak->GetXDiff());
			h_prop_good->Fill(peak->GetProp());
			h_rot_good->Fill(peak->GetRotation());
			h_chi2_good->Fill(peak->GetChi2());
			h_distvar_good->Fill(peak->GetDistVar());
			h_theta_good->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
			h_rmsratio_good->Fill(peak->GetSigmaY()/peak->GetSigmaX());
			h_radius_good->Fill(peakradius);
			h_rmsdistortion_good->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius)/(peak->GetSigmaY()/peak->GetSigmaX()));
			h_maxprop_good->Fill(peak->GetMax(),peak->GetProp());
			h_maxnpes_good->Fill(peak->GetMax(),peak->GetNPEs());
			h_npesprop_good->Fill(peak->GetNPEs(),peak->GetProp());
			h_radiustheta_good->Fill(peakradius,TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
		}
		else {
			if(flaggoodpeak){
				h_max_bad->Fill(peak->GetMax());
				h_nPE_bad->Fill(peak->GetNPEs());
				h_sigma_bad->Fill(peak->GetSigma());
				h_xdiff_bad->Fill(peak->GetXDiff());
				h_prop_bad->Fill(peak->GetProp());
				h_rot_bad->Fill(peak->GetRotation());
				h_chi2_bad->Fill(peak->GetChi2());
				h_distvar_bad->Fill(peak->GetDistVar());
				h_theta_bad->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
				h_rmsratio_bad->Fill(peak->GetSigmaY()/peak->GetSigmaX());
				h_radius_bad->Fill(peakradius);
				h_rmsdistortion_bad->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius)/(peak->GetSigmaY()/peak->GetSigmaX()));
				h_maxprop_bad->Fill(peak->GetMax(),peak->GetProp());
				h_maxnpes_bad->Fill(peak->GetMax(),peak->GetNPEs());
				h_npesprop_bad->Fill(peak->GetNPEs(),peak->GetProp());
				h_radiustheta_bad->Fill(peakradius,TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
			}
		}
	}
	
	h_max_tot->Add(h_max_good,h_max_bad);
	h_nPE_tot->Add(h_nPE_good,h_nPE_bad);
	h_sigma_tot->Add(h_sigma_good,h_sigma_bad);
	h_xdiff_tot->Add(h_xdiff_good,h_xdiff_bad);
	h_prop_tot->Add(h_prop_good,h_prop_bad);
	h_rot_tot->Add(h_rot_good,h_rot_bad);
	h_chi2_tot->Add(h_chi2_good,h_chi2_bad);
	
	h_distvar_tot->Add(h_distvar_good,h_distvar_bad);
	h_theta_tot->Add(h_theta_good,h_theta_bad);
	h_rmsratio_tot->Add(h_rmsratio_good,h_rmsratio_bad);
	h_radius_tot->Add(h_radius_good,h_radius_bad);
	h_rmsdistortion_tot->Add(h_rmsdistortion_good,h_rmsdistortion_bad);
	h_maxprop_tot->Add(h_maxprop_good,h_maxprop_bad);
	h_maxnpes_tot->Add(h_maxnpes_good,h_maxnpes_bad);
	h_npesprop_tot->Add(h_npesprop_good,h_npesprop_bad);
	h_radiustheta_tot->Add(h_radiustheta_good,h_radiustheta_bad);
	//save canvases
	
	TLegend *leg=new TLegend(0.55,0.55,0.85,0.85);
	leg->AddEntry(h_max_tot,"All Peaks","lp");
	leg->AddEntry(h_max_good,"Single Peaks","lp");
	leg->AddEntry(h_max_bad,"Non-single Peaks","lp");
	leg->SetFillColor(0);
	leg->SetTextAlign(12);
	leg->SetBorderSize(0);
	
	gStyle->SetOptTitle(0);
	TCanvas *can_max=new TCanvas(Form("can_max_hpd%i_col%i",nhpd,ncol),"can_max");
	h_max_tot->SetMarkerStyle(20);
	h_max_good->SetMarkerStyle(21);
	h_max_bad->SetMarkerStyle(22);
	h_max_good->SetMarkerColor(4);
	h_max_bad->SetMarkerColor(2);
	h_max_tot->Draw("hist");
	h_max_good->Draw("sames hist");
	h_max_bad->Draw("sames hist");
	leg->Draw();
	can_max->Write();
	TCanvas *can_nPE=new TCanvas(Form("can_nPE_hpd%i_col%i",nhpd,ncol),"can_nPE");
	h_nPE_tot->SetMarkerStyle(20);
	h_nPE_good->SetMarkerStyle(21);
	h_nPE_bad->SetMarkerStyle(22);
	h_nPE_good->SetMarkerColor(4);
	h_nPE_bad->SetMarkerColor(2);
	h_nPE_tot->Draw("hist");
	h_nPE_good->Draw("sames hist");
	h_nPE_bad->Draw("sames hist");
	leg->Draw();
	can_nPE->Write();
	TCanvas *can_sigma=new TCanvas(Form("can_sigma_hpd%i_col%i",nhpd,ncol),"can_sigma");
	h_sigma_tot->SetMarkerStyle(20);
	h_sigma_good->SetMarkerStyle(21);
	h_sigma_bad->SetMarkerStyle(22);
	h_sigma_good->SetMarkerColor(4);
	h_sigma_bad->SetMarkerColor(2);
	h_sigma_tot->Draw("hist");
	h_sigma_good->Draw("sames hist");
	h_sigma_bad->Draw("sames hist");
	leg->Draw();
	can_sigma->Write();
	TCanvas *can_xdiff=new TCanvas(Form("can_xdiff_hpd%i_col%i",nhpd,ncol),"can_xdiff");
	h_xdiff_tot->SetMarkerStyle(20);
	h_xdiff_good->SetMarkerStyle(21);
	h_xdiff_bad->SetMarkerStyle(22);
	h_xdiff_good->SetMarkerColor(4);
	h_xdiff_bad->SetMarkerColor(2);
	h_xdiff_tot->Draw("hist");
	h_xdiff_good->Draw("sames hist");
	h_xdiff_bad->Draw("sames hist");
	leg->Draw();
	can_xdiff->Write();
	TCanvas *can_prop=new TCanvas(Form("can_prop_hpd%i_col%i",nhpd,ncol),"can_prop");
	h_prop_tot->SetMarkerStyle(20);
	h_prop_good->SetMarkerStyle(21);
	h_prop_bad->SetMarkerStyle(22);
	h_prop_good->SetMarkerColor(4);
	h_prop_bad->SetMarkerColor(2);
	h_prop_tot->Draw("hist");
	h_prop_good->Draw("sames hist");
	h_prop_bad->Draw("sames hist");
	leg->Draw();
	can_prop->Write();
	TCanvas *can_rot=new TCanvas(Form("can_rot_hpd%i_col%i",nhpd,ncol),"can_rot");
	h_rot_tot->SetMarkerStyle(20);
	h_rot_good->SetMarkerStyle(21);
	h_rot_bad->SetMarkerStyle(22);
	h_rot_good->SetMarkerColor(4);
	h_rot_bad->SetMarkerColor(2);
	h_rot_tot->Draw("hist");
	h_rot_good->Draw("sames hist");
	h_rot_bad->Draw("sames hist");
	leg->Draw();
	can_rot->Write();
	TCanvas *can_chi2=new TCanvas(Form("can_chi2_hpd%i_col%i",nhpd,ncol),"can_chi2");
	h_chi2_tot->SetMarkerStyle(20);
	h_chi2_good->SetMarkerStyle(21);
	h_chi2_bad->SetMarkerStyle(22);
	h_chi2_good->SetMarkerColor(4);
	h_chi2_bad->SetMarkerColor(2);
	h_chi2_tot->Draw("hist");
	h_chi2_good->Draw("sames hist");
	h_chi2_bad->Draw("sames hist");
	leg->Draw();
	can_chi2->Write();
	
	TCanvas *can_distvar=new TCanvas(Form("can_distvar_hpd%i_col%i",nhpd,ncol),"can_distvar");
	h_distvar_tot->Draw("hist");
	h_distvar_good->Draw("sames hist");
	h_distvar_bad->Draw("sames hist");
	leg->Draw();
	can_distvar->Write();
	
	TCanvas *can_theta=new TCanvas(Form("can_theta_hpd%i_col%i",nhpd,ncol),"can_theta");
	h_theta_tot->Draw("hist");
	h_theta_good->Draw("sames hist");
	h_theta_bad->Draw("sames hist");
	leg->Draw();
	can_theta->Write();
	
	TCanvas *can_rmsratio=new TCanvas(Form("can_rmsratio_hpd%i_col%i",nhpd,ncol),"can_rmsratio");
	h_rmsratio_tot->Draw("hist");
	h_rmsratio_good->Draw("sames hist");
	h_rmsratio_bad->Draw("sames hist");
	leg->Draw();
	can_rmsratio->Write();
	
	TCanvas *can_radius=new TCanvas(Form("can_radius_hpd%i_col%i",nhpd,ncol),"can_radius");
	h_radius_tot->Draw("hist");
	h_radius_good->Draw("sames hist");
	h_radius_bad->Draw("sames hist");
	leg->Draw();
	can_radius->Write();
	
	TCanvas *can_rmsdistortion=new TCanvas(Form("can_rmsdistortion_hpd%i_col%i",nhpd,ncol),"can_rmsdistortion");
	h_rmsdistortion_tot->Draw("hist");
	h_rmsdistortion_good->Draw("sames hist");
	h_rmsdistortion_bad->Draw("sames hist");
	leg->Draw();
	can_rmsdistortion->Write();
	
	TCanvas *can_maxprop=new TCanvas(Form("can_maxprop_hpd%i_col%i",nhpd,ncol),"can_maxprop");
	h_maxprop_tot->Draw("box");
	h_maxprop_good->Draw("box sames");
	h_maxprop_bad->Draw("box sames");
	leg->Draw();
	can_maxprop->Write();
	
	TCanvas *can_maxnpes=new TCanvas(Form("can_maxnpes_hpd%i_col%i",nhpd,ncol),"can_maxnpes");
	h_maxnpes_tot->Draw("box");
	h_maxnpes_good->Draw("box sames");
	h_maxnpes_bad->Draw("box sames");
	leg->Draw();
	can_maxnpes->Write();
	
	TCanvas *can_npesprop=new TCanvas(Form("can_npesprop_hpd%i_col%i",nhpd,ncol),"can_npesprop");
	h_npesprop_tot->Draw("box");
	h_npesprop_good->Draw("box sames");
	h_npesprop_bad->Draw("box sames");
	leg->Draw();
	can_npesprop->Write();
	
	TCanvas *can_radiustheta=new TCanvas(Form("can_radiustheta_hpd%i_col%i",nhpd,ncol),"can_radiustheta");
	h_radiustheta_tot->Draw("box");
	h_radiustheta_good->Draw("box sames");
	h_radiustheta_bad->Draw("box sames");
	leg->Draw();
	can_radiustheta->Write();
	gStyle->SetOptTitle(1);
}

void CreateDistributions(TList *hpds){
	int i=0,k=0,l=0;
	
	double max_max=140, min_max=0;
	double max_npes=5000, min_npes=0;
	double max_sigma=5, min_sigma=0;
	double max_prop=1, min_prop=0;
	double max_rot=0.2, min_rot=-0.2;
	double max_rms=2, min_rms=0;
	double max_chi2=10, min_chi2=0;
	double max_distvar=20, min_distvar=0;
	double max_newdistvar=20, min_newdistvar=0;
	double max_theta=3.2, min_theta=-3.2;
	double max_rmsratio=5, min_rmsratio=0;
	double max_radius=20, min_radius=0;
	double max_rmsdistortion=5, min_rmsdistortion=-5;
	double max_xdiff=20, min_xdiff=0;
	//good peaks
	TH1D *h_max_good=new TH1D("h_max_good","Good Peaks Max;max;",100,min_max,max_max);
	h_max_good->SetLineColor(4);
	h_max_good->SetLineWidth(2);
	TH1D *h_nPE_good=new TH1D("h_nPE_good","Good Peaks nPE;nPE;",100,min_npes,max_npes);
	h_nPE_good->SetLineColor(4);
	h_nPE_good->SetLineWidth(2);
	TH1D *h_sigma_good=new TH1D("h_sigma_good","Good Peaks #sigma ;#sigma (LHCb pixels);",100,min_sigma,max_sigma);
	h_sigma_good->SetLineColor(4);
	h_sigma_good->SetLineWidth(2);
	TH1D *h_xdiff_good=new TH1D("h_xdiff_good","Good Peaks #xdiff ;#xdiff (LHCb pixels);",100,min_xdiff,max_xdiff);
	h_xdiff_good->SetLineColor(4);
	h_xdiff_good->SetLineWidth(2);
	TH1D *h_prop_good=new TH1D("h_prop_good","Good Peaks Proportion; proportion;",100,min_prop,max_prop);
	h_prop_good->SetLineColor(4);
	h_prop_good->SetLineWidth(2);
	TH1D *h_rot_good=new TH1D("h_rot_good","Good Peaks Rotation; rotation;",800,min_rot,max_rot);
	h_rot_good->SetLineColor(4);
	h_rot_good->SetLineWidth(2);
	TH1D *h_rms_good_X=new TH1D("h_rms_good_X","Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	h_rms_good_X->SetLineColor(4);
	h_rms_good_X->SetLineWidth(2);
	TH1D *h_rms_good_Y_=new TH1D("h_rms_good_Y_","Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	h_rms_good_Y_->SetLineColor(4);
	h_rms_good_Y_->SetLineWidth(2);
	TH1D *h_chi2_good=new TH1D("h_chi2_good","Good Peaks #chi^{2}; #chi^{2};",100,min_chi2,max_chi2);
	h_chi2_good->SetLineColor(4);
	h_chi2_good->SetLineWidth(2);
	TH1D *h_distvar_good=new TH1D("h_distvar_good","Good Peaks dist var; dist var;",400,min_distvar,max_distvar);
	h_distvar_good->SetLineColor(4);
	h_distvar_good->SetLineWidth(2);
	TH1D *h_theta_good=new TH1D("h_theta_good","Good Peaks #theta; #theta (rad);",100,min_theta,max_theta);
	h_theta_good->SetLineColor(4);
	h_theta_good->SetLineWidth(2);
	TH1D *h_rmsratio_good=new TH1D("h_rmsratio_good","Good Peaks rmsratio; rmsy/rmsx;",200,min_rmsratio,max_rmsratio);
	h_rmsratio_good->SetLineColor(4);
	h_rmsratio_good->SetLineWidth(2);
	TH1D *h_radius_good=new TH1D("h_radius_good","Good Peaks radius; radius (LHCb Pixels);",200,min_radius,max_radius);
	h_radius_good->SetLineColor(4);
	h_radius_good->SetLineWidth(2);
	TH1D *h_rmsdistortion_good=new TH1D("h_rmsdistortion_good","Good Peaks rms distortion; rms distortion;",200,min_rmsdistortion,max_rmsdistortion);
	h_rmsdistortion_good->SetLineColor(4);
	h_rmsdistortion_good->SetLineWidth(2);
	TH1D *h_newdistvar_good=new TH1D("h_newdistvar_good","Good Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	h_newdistvar_good->SetLineColor(4);
	h_newdistvar_good->SetLineWidth(2);
	TH2D *h_maxprop_good=new TH2D("h_maxprop_good","Good peaks maxprop",25,min_max,max_max,25,min_prop,max_prop);
	h_maxprop_good->SetMarkerColor(4);
	h_maxprop_good->SetLineColor(4);
// 	h_maxprop_good->SetLineWidth(2);
	TH2D *h_maxnpes_good=new TH2D("h_maxnpes_good","Good peaks maxnpes",25,min_max,max_max,25,min_npes,max_npes);
	h_maxnpes_good->SetMarkerColor(4);
	h_maxnpes_good->SetLineColor(4);
	TH2D *h_npesprop_good=new TH2D("h_npesprop_good","Good peaks npesprop",25,min_npes,max_npes,25,min_prop,max_prop);
	h_npesprop_good->SetMarkerColor(4);
	h_npesprop_good->SetLineColor(4);
	TH2D *h_radiustheta_good=new TH2D("h_radiustheta_good","Good peaks radiustheta",25,min_radius,max_radius,25,min_theta,max_theta);
	h_radiustheta_good->SetMarkerColor(4);
	h_radiustheta_good->SetLineColor(4);
	
	//badpeaks
	TH1D *h_max_bad=new TH1D("h_max_bad","Bad Peaks Max;max;",100,min_max,max_max);
	h_max_bad->SetLineColor(2);
	h_max_bad->SetLineWidth(2);
	TH1D *h_nPE_bad=new TH1D("h_nPE_bad","Bad Peaks nPE;nPE;",100,min_npes,max_npes);
	h_nPE_bad->SetLineColor(2);
	h_nPE_bad->SetLineWidth(2);
	TH1D *h_sigma_bad=new TH1D("h_sigma_bad","Bad Peaks #sigma ;#sigma (LHCb pixels);",100,min_sigma,max_sigma);
	h_sigma_bad->SetLineColor(2);
	h_sigma_bad->SetLineWidth(2);
	TH1D *h_xdiff_bad=new TH1D("h_xdiff_bad","Bad Peaks #xdiff ;#xdiff (LHCb pixels);",100,min_xdiff,max_xdiff);
	h_xdiff_bad->SetLineColor(2);
	h_xdiff_bad->SetLineWidth(2);
	TH1D *h_prop_bad=new TH1D("h_prop_bad","Bad Peaks Proportion; proportion;",100,min_prop,max_prop);
	h_prop_bad->SetLineColor(2);
	h_prop_bad->SetLineWidth(2);
	TH1D *h_rot_bad=new TH1D("h_rot_bad","Bad Peaks Rotation; rotation;",800,min_rot,max_rot);
	h_rot_bad->SetLineColor(2);
	h_rot_bad->SetLineWidth(2);
	TH1D *h_rms_bad_X=new TH1D("h_rms_bad_X","Good Peaks RMS; RMS x (LHCb pixels);",100,min_rms,max_rms);
	h_rms_bad_X->SetLineColor(2);
	h_rms_bad_X->SetLineWidth(2);
	TH1D *h_rms_bad_Y=new TH1D("h_rms_bad_Y","Good Peaks RMS; RMS y (LHCb pixels);",100,min_rms,max_rms);
	h_rms_bad_Y->SetLineColor(2);
	h_rms_bad_Y->SetLineWidth(2);
	TH1D *h_chi2_bad=new TH1D("h_chi2_bad","Bad Peaks #chi^{2}; #chi^{2};",100,min_chi2,max_chi2);
	h_chi2_bad->SetLineColor(2);
	h_chi2_bad->SetLineWidth(2);
	TH1D *h_distvar_bad=new TH1D("h_distvar_bad","Bad Peaks dist var; dist var;",400,min_distvar,max_distvar);
	h_distvar_bad->SetLineColor(2);
	h_distvar_bad->SetLineWidth(2);
	TH1D *h_theta_bad=new TH1D("h_theta_bad","Bad Peaks #theta; #theta (rad);",100,min_theta,max_theta);
	h_theta_bad->SetLineColor(2);
	h_theta_bad->SetLineWidth(2);
	TH1D *h_rmsratio_bad=new TH1D("h_rmsratio_bad","Bad Peaks rmsratio; rmsy/rmsx;",200,min_rmsratio,max_rmsratio);
	h_rmsratio_bad->SetLineColor(2);
	h_rmsratio_bad->SetLineWidth(2);
	TH1D *h_radius_bad=new TH1D("h_radius_bad","Bad Peaks radius; radius (LHCb Pixels);",200,min_radius,max_radius);
	h_radius_bad->SetLineColor(2);
	h_radius_bad->SetLineWidth(2);
	TH1D *h_rmsdistortion_bad=new TH1D("h_rmsdistortion_bad","Bad Peaks rms distortion; rms distortion;",200,min_rmsdistortion,max_rmsdistortion);
	h_rmsdistortion_bad->SetLineColor(2);
	h_rmsdistortion_bad->SetLineWidth(2);
	TH1D *h_newdistvar_bad=new TH1D("h_newdistvar_bad","Bad Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	h_newdistvar_bad->SetLineColor(2);
	h_newdistvar_bad->SetLineWidth(2);
	TH2D *h_maxprop_bad=new TH2D("h_maxprop_bad","Bad peaks maxprop",25,min_max,max_max,25,min_prop,max_prop);
	h_maxprop_bad->SetMarkerColor(2);
	h_maxprop_bad->SetLineColor(2);
	TH2D *h_maxnpes_bad=new TH2D("h_maxnpes_bad","Bad peaks maxnpes",25,min_max,max_max,25,min_npes,max_npes);
	h_maxnpes_bad->SetMarkerColor(2);
	h_maxnpes_bad->SetLineColor(2);
	TH2D *h_npesprop_bad=new TH2D("h_npesprop_bad","Bad peaks npesprop",25,min_npes,max_npes,25,min_prop,max_prop);
	h_npesprop_bad->SetMarkerColor(2);
	h_npesprop_bad->SetLineColor(2);
	TH2D *h_radiustheta_bad=new TH2D("h_radiustheta_bad","Bad peaks radiustheta",25,min_radius,max_radius,25,min_theta,max_theta);
	h_radiustheta_bad->SetMarkerColor(2);
	h_radiustheta_bad->SetLineColor(2);
	
	//totpeaks
	TH1D *h_max_tot=new TH1D("h_max_tot","Tot Peaks Max;M;",100,min_max,max_max);
	TH1D *h_nPE_tot=new TH1D("h_nPE_tot","Tot Peaks nPE;N_{hits};",100,min_npes,max_npes);
	TH1D *h_sigma_tot=new TH1D("h_sigma_tot","Tot Peaks #sigma ;#sigma (LHCb pixels);",100,min_sigma,max_sigma);
	TH1D *h_xdiff_tot=new TH1D("h_xdiff_tot","Tot Peaks #xdiff ;x_{diff} (LHCb pixels);",100,min_xdiff,max_xdiff);
	TH1D *h_prop_tot=new TH1D("h_prop_tot","Tot Peaks Proportion; proportion;",100,min_prop,max_prop);
	TH1D *h_rot_tot=new TH1D("h_rot_tot","Tot Peaks Rotation; rotation;",800,min_rot,max_rot);
	TH1D *h_chi2_tot=new TH1D("h_chi2_tot","Tot Peaks #chi^{2}; #chi^{2};",100,min_chi2,max_chi2);
	TH1D *h_distvar_tot=new TH1D("h_distvar_tot","Tot Peaks dist var; dist var;",400,min_distvar,max_distvar);
	TH1D *h_theta_tot=new TH1D("h_theta_tot","Tot Peaks #theta; #theta (rad);",100,min_theta,max_theta);
	TH1D *h_rmsratio_tot=new TH1D("h_rmsratio_tot","Tot Peaks rmsratio; rmsy/rmsx;",200,min_rmsratio,max_rmsratio);
	TH1D *h_radius_tot=new TH1D("h_radius_tot","Tot Peaks radius; radius (LHCb Pixels);",200,min_radius,max_radius);
	TH1D *h_rmsdistortion_tot=new TH1D("h_rmsdistortion_tot","Tot Peaks rms distortion; rms distortion;",200,min_rmsdistortion,max_rmsdistortion);
	TH1D *h_newdistvar_tot=new TH1D("h_newdistvar_tot","Tot Peaks new dist var; newdistvar;",200,min_newdistvar,max_newdistvar);
	TH2D *h_maxprop_tot=new TH2D("h_maxprop_tot","Tot peaks maxprop",25,min_max,max_max,25,min_prop,max_prop);
	TH2D *h_maxnpes_tot=new TH2D("h_maxnpes_tot","Tot peaks maxnpes",25,min_max,max_max,25,min_npes,max_npes);
	TH2D *h_npesprop_tot=new TH2D("h_npesprop_tot","Tot peaks npesprop",25,min_npes,max_npes,25,min_prop,max_prop);
	TH2D *h_radiustheta_tot=new TH2D("h_radiustheta_tot","Tot peaks radiustheta",25,min_radius,max_radius,25,min_theta,max_theta);
	TH1D *h_rms_tot_X=new TH1D("h_rms_tot_X","Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	TH1D *h_rms_tot_Y=new TH1D("h_rms_tot_Y","Good Peaks RMS; RMS (LHCb pixels);",100,min_rms,max_rms);
	
	
	h_max_tot->SetLineWidth(2);
	h_nPE_tot->SetLineWidth(2);
	h_sigma_tot->SetLineWidth(2);
	h_xdiff_tot->SetLineWidth(2);
	h_prop_tot->SetLineWidth(2);
	h_rot_tot->SetLineWidth(2);
	h_chi2_tot->SetLineWidth(2);
	h_distvar_tot->SetLineWidth(2);
	h_theta_tot->SetLineWidth(2);
	h_rmsratio_tot->SetLineWidth(2);
	h_radius_tot->SetLineWidth(2);
	h_rmsdistortion_tot->SetLineWidth(2);
	h_newdistvar_tot->SetLineWidth(2);
	
	
	//checkdist
	TH1D *h_dist_good=new TH1D("pattern_step_global_good","Good peaks, pattern step; Peak distance (LHCb pixel);",100,0,5);
	h_dist_good->SetLineColor(4);
	TH1D *h_dist_bad=new TH1D("pattern_step_global_bad","Bad peaks, pattern step; Peak distance (LHCb pixel);",100,0,5);
	h_dist_bad->SetLineColor(2);
	TH1D *h_dist_tot=new TH1D("pattern_step_global_tot","Tot peaks, pattern step; Peak distance (LHCb pixel);",100,0,5);
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		CreateDistributions_SingleHPDs(hpd);
		TList *peaklist=(TList*)hpd->GetHitmap();
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			double peakradius=TMath::Sqrt((peak->GetMuX()-hpd->GetMeanX())*(peak->GetMuX()-hpd->GetMeanX())+(peak->GetMuY()-hpd->GetMeanY())*(peak->GetMuY()-hpd->GetMeanY()));
// 			cout << peak->GetMuX() << "\t" << peak->GetMeanX() << endl;
			//get rid of fake peaks
			bool flaggoodpeak=true;
			if(peak->GetMax()<10 || peak->GetSigma()<0.5 || peak->GetChi2()<0.5){
				peak->SetPassLevel("DIS_STEP1");
				flaggoodpeak=false;
			}
			
			//fill histos
			
			if(peak->IsSingle() && hpd->GetPassLevel()!="IFB" && flaggoodpeak){
				h_max_good->Fill(peak->GetMax());
				h_nPE_good->Fill(peak->GetNPEs());
				h_sigma_good->Fill(peak->GetSigma());
				h_xdiff_good->Fill(peak->GetXDiff());
				h_prop_good->Fill(peak->GetProp());
				h_rot_good->Fill(peak->GetRotation());
				h_chi2_good->Fill(peak->GetChi2());
				h_distvar_good->Fill(peak->GetDistVar());
				h_theta_good->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
				h_rmsratio_good->Fill(peak->GetSigmaY()/peak->GetSigmaX());
				h_radius_good->Fill(peakradius);
				h_rmsdistortion_good->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius)/(peak->GetSigmaY()/peak->GetSigmaX()));
				h_maxprop_good->Fill(peak->GetMax(),peak->GetProp());
				h_maxnpes_good->Fill(peak->GetMax(),peak->GetNPEs());
				h_npesprop_good->Fill(peak->GetNPEs(),peak->GetProp());
				h_radiustheta_good->Fill(peakradius,TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
				for(l=k+1; l<hpd->GetNPeaks(); l++){
					CleanPeak *peak_tmp=(CleanPeak*)peaklist->At(l);
					if(abs(peak_tmp->GetCalibStep()-peak->GetCalibStep())==1){
						h_dist_good->Fill(TMath::Sqrt((peak_tmp->GetMuX()-peak->GetMuX())*(peak_tmp->GetMuX()-peak->GetMuX())+(peak_tmp->GetMuY()-peak->GetMuY())*(peak_tmp->GetMuY()-peak->GetMuY())));
					}
				}
			}
			else {
				if(flaggoodpeak){
					h_max_bad->Fill(peak->GetMax());
					h_nPE_bad->Fill(peak->GetNPEs());
					h_sigma_bad->Fill(peak->GetSigma());
					h_xdiff_bad->Fill(peak->GetXDiff());
					h_prop_bad->Fill(peak->GetProp());
					h_rot_bad->Fill(peak->GetRotation());
					h_chi2_bad->Fill(peak->GetChi2());
					h_distvar_bad->Fill(peak->GetDistVar());
					h_theta_bad->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
					h_rmsratio_bad->Fill(peak->GetSigmaY()/peak->GetSigmaX());
					h_radius_bad->Fill(peakradius);
					h_rmsdistortion_bad->Fill(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius)/(peak->GetSigmaY()/peak->GetSigmaX()));
					h_maxprop_bad->Fill(peak->GetMax(),peak->GetProp());
					h_maxnpes_bad->Fill(peak->GetMax(),peak->GetNPEs());
					h_npesprop_bad->Fill(peak->GetNPEs(),peak->GetProp());
					h_radiustheta_bad->Fill(peakradius,TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius));
					for(l=k+1; l<hpd->GetNPeaks(); l++){
						CleanPeak *peak_tmp2=(CleanPeak*)peaklist->At(l);
						if(abs(peak_tmp2->GetCalibStep()-peak->GetCalibStep())==1){
							h_dist_good->Fill(TMath::Sqrt((peak_tmp2->GetMuX()-peak->GetMuX())*(peak_tmp2->GetMuX()-peak->GetMuX())+(peak_tmp2->GetMuY()-peak->GetMuY())*(peak_tmp2->GetMuY()-peak->GetMuY())));
						}
					}
				}
			}
		}
	}
	h_max_tot->Add(h_max_good,h_max_bad);
	h_nPE_tot->Add(h_nPE_good,h_nPE_bad);
	h_sigma_tot->Add(h_sigma_good,h_sigma_bad);
	h_xdiff_tot->Add(h_xdiff_good,h_xdiff_bad);
	h_prop_tot->Add(h_prop_good,h_prop_bad);
	h_rot_tot->Add(h_rot_good,h_rot_bad);
	h_chi2_tot->Add(h_chi2_good,h_chi2_bad);
	h_dist_tot->Add(h_dist_good,h_dist_bad);
	h_distvar_tot->Add(h_distvar_good,h_distvar_bad);
	h_theta_tot->Add(h_theta_good,h_theta_bad);
	h_rmsratio_tot->Add(h_rmsratio_good,h_rmsratio_bad);
	h_radius_tot->Add(h_radius_good,h_radius_bad);
	h_rmsdistortion_tot->Add(h_rmsdistortion_good,h_rmsdistortion_bad);
	h_maxprop_tot->Add(h_maxprop_good,h_maxprop_bad);
	h_maxnpes_tot->Add(h_maxnpes_good,h_maxnpes_bad);
	h_npesprop_tot->Add(h_npesprop_good,h_npesprop_bad);
	h_radiustheta_tot->Add(h_radiustheta_good,h_radiustheta_bad);
	//save canvases
	
	TLegend *leg=new TLegend(0.55,0.55,0.85,0.85);
	leg->AddEntry(h_max_tot,"All Peaks","lp");
	leg->AddEntry(h_max_good,"Single Peaks","lp");
	leg->AddEntry(h_max_bad,"Non-single Peaks","lp");
	leg->SetFillColor(0);
	leg->SetTextAlign(12);
	leg->SetBorderSize(0);
	
	gStyle->SetOptTitle(0);
	TCanvas *can_max=new TCanvas("can_max","can_max");
	h_max_tot->SetMarkerStyle(20);
	h_max_good->SetMarkerStyle(21);
	h_max_bad->SetMarkerStyle(22);
	h_max_good->SetMarkerColor(4);
	h_max_bad->SetMarkerColor(2);
	h_max_tot->Draw();
	h_max_good->Draw("sames");
	h_max_bad->Draw("sames");
	leg->Draw();
	can_max->Write();
	TCanvas *can_nPE=new TCanvas("can_nPE","can_nPE");
	h_nPE_tot->SetMarkerStyle(20);
	h_nPE_good->SetMarkerStyle(21);
	h_nPE_bad->SetMarkerStyle(22);
	h_nPE_good->SetMarkerColor(4);
	h_nPE_bad->SetMarkerColor(2);
	h_nPE_tot->Draw();
	h_nPE_good->Draw("sames");
	h_nPE_bad->Draw("sames");
	leg->Draw();
	can_nPE->Write();
	TCanvas *can_sigma=new TCanvas("can_sigma","can_sigma");
	h_sigma_tot->SetMarkerStyle(20);
	h_sigma_good->SetMarkerStyle(21);
	h_sigma_bad->SetMarkerStyle(22);
	h_sigma_good->SetMarkerColor(4);
	h_sigma_bad->SetMarkerColor(2);
	h_sigma_tot->Draw();
	h_sigma_good->Draw("sames");
	h_sigma_bad->Draw("sames");
	leg->Draw();
	can_sigma->Write();
	TCanvas *can_xdiff=new TCanvas("can_xdiff","can_xdiff");
	h_xdiff_tot->SetMarkerStyle(20);
	h_xdiff_good->SetMarkerStyle(21);
	h_xdiff_bad->SetMarkerStyle(22);
	h_xdiff_good->SetMarkerColor(4);
	h_xdiff_bad->SetMarkerColor(2);
	h_xdiff_tot->Draw();
	h_xdiff_good->Draw("sames");
	h_xdiff_bad->Draw("sames");
	leg->Draw();
	can_xdiff->Write();
	TCanvas *can_prop=new TCanvas("can_prop","can_prop");
	h_prop_tot->SetMarkerStyle(20);
	h_prop_good->SetMarkerStyle(21);
	h_prop_bad->SetMarkerStyle(22);
	h_prop_good->SetMarkerColor(4);
	h_prop_bad->SetMarkerColor(2);
	h_prop_tot->Draw();
	h_prop_good->Draw("sames");
	h_prop_bad->Draw("sames");
	leg->Draw();
	can_prop->Write();
	TCanvas *can_rot=new TCanvas("can_rot","can_rot");
	h_rot_tot->SetMarkerStyle(20);
	h_rot_good->SetMarkerStyle(21);
	h_rot_bad->SetMarkerStyle(22);
	h_rot_good->SetMarkerColor(4);
	h_rot_bad->SetMarkerColor(2);
	h_rot_tot->Draw();
	h_rot_good->Draw("sames");
	h_rot_bad->Draw("sames");
	leg->Draw();
	can_rot->Write();
	TCanvas *can_chi2=new TCanvas("can_chi2","can_chi2");
	h_chi2_tot->SetMarkerStyle(20);
	h_chi2_good->SetMarkerStyle(21);
	h_chi2_bad->SetMarkerStyle(22);
	h_chi2_good->SetMarkerColor(4);
	h_chi2_bad->SetMarkerColor(2);
	h_chi2_tot->Draw();
	h_chi2_good->Draw("sames");
	h_chi2_bad->Draw("sames");
	leg->Draw();
	can_chi2->Write();
	TCanvas *can_dist=new TCanvas("can_dist","can_dist");
	h_dist_tot->SetMarkerStyle(20);
	h_dist_good->SetMarkerStyle(21);
	h_dist_bad->SetMarkerStyle(22);
	h_dist_good->SetMarkerColor(4);
	h_dist_bad->SetMarkerColor(2);
	h_dist_tot->Draw();
	h_dist_good->Draw("sames");
	h_dist_bad->Draw("sames");
	leg->Draw();
	can_dist->Write();
	
	TCanvas *can_distvar=new TCanvas("can_distvar","can_distvar");
	h_distvar_tot->Draw();
	h_distvar_good->Draw("sames");
	h_distvar_bad->Draw("sames");
	leg->Draw();
	can_distvar->Write();
	
	TCanvas *can_theta=new TCanvas("can_theta","can_theta");
	h_theta_tot->Draw();
	h_theta_good->Draw("sames");
	h_theta_bad->Draw("sames");
	leg->Draw();
	can_theta->Write();
	
	TCanvas *can_rmsratio=new TCanvas("can_rmsratio","can_rmsratio");
	h_rmsratio_tot->Draw();
	h_rmsratio_good->Draw("sames");
	h_rmsratio_bad->Draw("sames");
	leg->Draw();
	can_rmsratio->Write();
	
	TCanvas *can_radius=new TCanvas("can_radius","can_radius");
	h_radius_tot->Draw();
	h_radius_good->Draw("sames");
	h_radius_bad->Draw("sames");
	leg->Draw();
	can_radius->Write();
	
	TCanvas *can_rmsdistortion=new TCanvas("can_rmsdistortion","can_rmsdistortion");
	h_rmsdistortion_tot->Draw();
	h_rmsdistortion_good->Draw("sames");
	h_rmsdistortion_bad->Draw("sames");
	leg->Draw();
	can_rmsdistortion->Write();
	
	TCanvas *can_maxprop=new TCanvas("can_maxprop","can_maxprop");
	h_maxprop_tot->Draw("box");
	h_maxprop_good->Draw("box sames");
	h_maxprop_bad->Draw("box sames");
	leg->Draw();
	can_maxprop->Write();
	
	TCanvas *can_maxnpes=new TCanvas("can_maxnpes","can_maxnpes");
	h_maxnpes_tot->Draw("box");
	h_maxnpes_good->Draw("box sames");
	h_maxnpes_bad->Draw("box sames");
	leg->Draw();
	can_maxnpes->Write();
	
	TCanvas *can_npesprop=new TCanvas("can_npesprop","can_npesprop");
	h_npesprop_tot->Draw("box");
	h_npesprop_good->Draw("box sames");
	h_npesprop_bad->Draw("box sames");
	leg->Draw();
	can_npesprop->Write();
	
	TCanvas *can_radiustheta=new TCanvas("can_radiustheta","can_radiustheta");
	h_radiustheta_tot->Draw("box");
	h_radiustheta_good->Draw("box sames");
	h_radiustheta_bad->Draw("box sames");
	leg->Draw();
	can_radiustheta->Write();
	gStyle->SetOptTitle(1);
}

void CreateDistributions_trees(TList *hpds, TString cuts=""){
	int i=0,k=0,l=0;
	
	//totpeaks
	TTree *h_tree=new TTree("h_tree_4fit","Peaks Vars");
	
	double t_max;
	h_tree->Branch("max",&t_max,"max/D");
	double t_nPE;
	h_tree->Branch("npes",&t_nPE,"nPE/D");
	double t_sigma;
	h_tree->Branch("sigma",&t_sigma,"sigma/D");
	double t_xdiff;
	h_tree->Branch("xdiff",&t_xdiff,"xdiff/D");
	double t_prop;
	h_tree->Branch("prop",&t_prop,"prop/D");
	double t_rot;
	h_tree->Branch("rot",&t_rot,"rot/D");
	double t_chi2;
	h_tree->Branch("chi2",&t_chi2,"chi2/D");
	double t_distvar;
	h_tree->Branch("distvar",&t_distvar,"distvar/D");
	double t_theta;
	h_tree->Branch("theta",&t_theta,"theta/D");
	double t_radius;
	h_tree->Branch("radius",&t_radius,"radius/D");
	double t_rmsdistortion;
	h_tree->Branch("rmsdistortion",&t_rmsdistortion,"rmsdistortion/D");
	double t_newdistvar;
	h_tree->Branch("newdistvar",&t_newdistvar,"newdistvar/D");
	
// 	cout << "crashing here"<<endl;
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		TList *peaklist=(TList*)hpd->GetHitmap();
// 		cout << "crashing here loop"<<endl;
		for(k=0; k<peaklist->GetSize(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
// 			cout << hpd->GetMeanX()<< endl;
			double peakradius = TMath::Sqrt(fabs((peak->GetMuX()-hpd->GetMeanX())*(peak->GetMuX()-hpd->GetMeanX())+(peak->GetMuY()-hpd->GetMeanY())*(peak->GetMuY()-hpd->GetMeanY())));
// 			cout << peak->GetMuX() << "\t" << hpd->GetMeanX() << endl;
// 			cout << peak->GetMuY() << "\t" << hpd->GetMeanY() << endl;
			//get rid of fake peaks
			bool flaggoodpeak=true;
			
			//fill histos
			
			if(hpd->GetPassLevel()!="IFB" && flaggoodpeak){
				t_max=peak->GetMax();
				t_nPE=peak->GetNPEs();
				t_sigma=peak->GetSigma();
				t_xdiff=peak->GetXDiff();
				t_prop=peak->GetProp();
				t_rot=peak->GetRotation();
				t_chi2=peak->GetChi2();
				t_distvar=peak->GetDistVar();
				t_theta=TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius);
				t_radius=peakradius;
				t_rmsdistortion=(TMath::ASin((peak->GetMuY()-hpd->GetMeanY())/peakradius)/(peak->GetSigmaY()/peak->GetSigmaX()));
				h_tree->Fill();
			}
		}
// 		peaklist->Clear("nodelete");
	}
// 	cout << "crashing here after"<<endl;
	h_tree->Write();
// 	cout << "crashing here end"<<endl;
// 	//Fit distributions
	FitMax(h_tree);
	FitNPEs(h_tree);
// 	double meandistance=FitInitDist(h_dist_good);
// 	CreateNVDistribution(hpds,meandistance);
	FitProp(h_tree);
	FitXDiff(h_tree);
// 	FitChi2(h_chi2_tot);
// 	FitSigma(h_sigma_tot);
// 	FitRot(h_rot_tot);
// 	FitDistVar(h_distvar_tot);
}



double getLL(double *values, TF1* funcs, double *integrals, int nvalues){
	double value=1;
	int i;
	for(i=0;i<nvalues;i++){
		if(funcs[i].Eval(values[i])>0) value=value*(funcs[i].Eval(values[i])/integrals[i]);
	}
	value=TMath::Log(value);
// 	cout << value << endl;
	return value;
}

void FillFuncList(TString *vars,TF1 *sigfun,double *sigint, TF1 *bkgfun,double *bkgint,int nvars, TFile* infile){
	int i=0;
	for(i=0;i<nvars;i++){
		double t1,t2;
		sigfun[i]=*((TF1*)infile->Get(Form("tf1_%s_sig",(vars[i]).Data())));
		sigfun[i].GetRange(t1,t2);
		sigint[i]=sigfun[i].Integral(t1,t2);
		bkgfun[i]=*((TF1*)infile->Get(Form("tf1_%s_bkg",(vars[i]).Data())));
		bkgfun[i].GetRange(t1,t2);
		bkgint[i]=sigfun[i].Integral(t1,t2);
	}
}

void BuildLL(TList *hpds, TFile *infile){
	
	int n_fun4like=4;
	
	TString* vars4like=new TString[n_fun4like];
	vars4like[0]="max";
	vars4like[1]="prop";
	vars4like[2]="npes";
	vars4like[3]="xdiff";
// 	vars4like[3]="distvar";
	
	double *sig_integrals=new double[n_fun4like];
	double *bkg_integrals=new double[n_fun4like];
	TF1 *Signals=new TF1[n_fun4like];
	TF1 *Backgrounds=new TF1[n_fun4like];
	
	FillFuncList(vars4like, Signals, sig_integrals, Backgrounds, bkg_integrals, n_fun4like, infile);
	
	
	int i,k;
	double *tmpvalues=new double[n_fun4like];
	cout << "Saving LL values"<< endl;
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		TList *peaklist=(TList*)hpd->GetHitmap();
// 		cout << "Processing HPD "<<hpd->GetHPD()<< " on  column "<<hpd->GetColumn()<<endl;
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			tmpvalues[0]=peak->GetMax();
			tmpvalues[1]=peak->GetProp();
			tmpvalues[2]=peak->GetNPEs();
			tmpvalues[3]=peak->GetXDiff();
			peak->SetLike1(getLL(tmpvalues,Signals,sig_integrals,n_fun4like));
			peak->SetLike2(getLL(tmpvalues,Backgrounds,bkg_integrals,n_fun4like));
			peak->SetLike3(getLL(tmpvalues,Signals,sig_integrals,n_fun4like)-getLL(tmpvalues,Backgrounds,bkg_integrals,n_fun4like));
// 			cout << peak->GetLike1() << endl;
		}
	}
	
	
	TCanvas *tmpcan=new TCanvas("pdfs_can","Pdfs");
	tmpcan->Divide(2,2);
	tmpcan->cd(1);
	Signals[0].Draw();
	Backgrounds[0].Draw("same");
	tmpcan->cd(2);
	Signals[1].Draw();
	Backgrounds[1].Draw("same");
	tmpcan->cd(3);
	Signals[2].Draw();
	Backgrounds[2].Draw("same");
	tmpcan->cd(4);
	Signals[3].Draw();
	Backgrounds[3].Draw("same");
	tmpcan->Write();
}

void CreateDLLDistributions(TList *hpds, int thebox){
	int i=0,k=0;
	//good peaks
	
	double dll_ll=-60, dll_ul=20;
	int dll_nbins=100;
	
	cout << "Filling DLL histo..." << endl;
	
	TH1F *h_LL1_good=new TH1F("LL1_good","LL1 good peaks; Signal LL;",dll_nbins,dll_ll,dll_ul);
	h_LL1_good->SetLineColor(4);
	TH1F *h_LL2_good=new TH1F("LL2_good","LL2 good peaks; Background LL;",dll_nbins,dll_ll,dll_ul);
	h_LL2_good->SetLineColor(4);
	TH1F *h_DLL_good=new TH1F("DLL_good","DLL good peaks; DLL;",dll_nbins,dll_ll,dll_ul);
	h_DLL_good->SetLineColor(4);
	TH1F *h_LL1_bad=new TH1F("LL1_bad","LL1 bad peaks; Signal LL;",dll_nbins,dll_ll,dll_ul);
	h_LL1_bad->SetLineColor(2);
	TH1F *h_LL2_bad=new TH1F("LL2_bad","LL2 bad peaks; Background LL;",dll_nbins,dll_ll,dll_ul);
	h_LL2_bad->SetLineColor(2);
	TH1F *h_DLL_bad=new TH1F("DLL_bad","DLL bad peaks; DLL;",dll_nbins,dll_ll,dll_ul);
	h_DLL_bad->SetLineColor(2);
	
	TH1F *h_LL1_tot=new TH1F("LL1_tot","LL1 tot peaks; Signal LL;",dll_nbins,dll_ll,dll_ul);
	TH1F *h_LL2_tot=new TH1F("LL2_tot","LL2 tot peaks; Background LL;",dll_nbins,dll_ll,dll_ul);
	TH1F *h_DLL_tot=new TH1F("DLL_tot","DLL tot peaks; DLL;",dll_nbins,dll_ll,dll_ul);
	
	
	TGraph *g_final=new TGraph();
	double VSpace=30;
	double HTilt=16;
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		TList *peaklist=(TList*)hpd->GetHitmap();
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			//fill histos
			if(peak->GetLike3()>DLL_CUT && peak->GetPassLevel()!="DIS_STEP1"){
				peak->SetPassLevel("OK");
				if(thebox==0) g_final->SetPoint(g_final->GetN(),peak->GetMuY()+32*(13-peak->GetHPD())+HTilt*(((int)peak->GetColumn())%2),peak->GetMuX()+VSpace*(6-peak->GetColumn()));
				if(thebox==1) g_final->SetPoint(g_final->GetN(),peak->GetMuY()+32*peak->GetHPD()+HTilt*((1+(int)peak->GetColumn())%2),peak->GetMuX()+VSpace*peak->GetColumn());
			}
			if(peak->IsSingle() && hpd->GetPassLevel()!="IFB" && peak->GetPassLevel()!="DIS_STEP1"){
				h_LL1_good->Fill(peak->GetLike1());
				h_LL2_good->Fill(peak->GetLike2());
				h_DLL_good->Fill(peak->GetLike3());
			}
			else{
				h_LL1_bad->Fill(peak->GetLike1());
				h_LL2_bad->Fill(peak->GetLike2());
				h_DLL_bad->Fill(peak->GetLike3());
			}
// 			cout << "yeahhh "<< peak->GetLike1() << endl;
			
			//fill histos
		}
	}
	h_LL1_tot->Add(h_LL1_good,h_LL1_bad);
	h_LL2_tot->Add(h_LL1_good,h_LL2_bad);
	h_DLL_tot->Add(h_DLL_good,h_DLL_bad);
	produceplot("Pre_Final", g_final, thebox);
}

void FlagPeakInCalib(TList * List){
	int i=0,k=0;
	
	for(k=0;k<nSTEPs;k++){
		double tmpDLL=DLL_CUT;
		int theindex=0;
		for(i=0;i<List->GetSize();i++){
			CleanPeak *peak=(CleanPeak*)List->At(i);
			if(peak->GetCalibStep()==k && peak->GetLike3()>tmpDLL){
				CleanPeak *peak2=(CleanPeak*)List->At(theindex);
				peak2->SetPassLevel("OK");
				tmpDLL=peak->GetLike3();
				peak->SetPassLevel("OKX2");
			}
		}
	}
}

void StripMultiPeakSteps(TList *hpds, int thebox){
	int i=0,k=0;
	
	cout<< "Flagging final peaks"<<endl;
	
	TGraph *g_final=new TGraph();
	double VSpace=30;
	double HTilt=16;
	int tmpstep=-1;
	int peakcounter=0;
	int currentpeak=0;
	double higherDLL=DLL_CUT-999;
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		TList *peaklist=(TList*)hpd->GetHitmap();
		tmpstep=-1;
		higherDLL=DLL_CUT-999;
		peakcounter=0;
		currentpeak=0;
		TList *MiniList=new TList();
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			if(peak->GetPassLevel()=="OK"){
				MiniList->Add(peak);
			}
		}
		FlagPeakInCalib(MiniList);
// 		MiniList->Clear("nodelete");
// 		peaklist->Clear("nodelete");
	}
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd2=(CleanHPD*)hpds->At(i);
		TList *peaklist2=(TList*)hpd2->GetHitmap();
		tmpstep=-1;
		higherDLL=DLL_CUT-999;
		peakcounter=0;
		currentpeak=0;
		for(k=0; k<hpd2->GetNPeaks(); k++){
			CleanPeak *peak2=(CleanPeak*)peaklist2->At(k);
			peakcounter++;
			if(peak2->GetPassLevel()=="OKX2"){
				if(thebox==0) g_final->SetPoint(g_final->GetN(),peak2->GetMuY()+32*(13-peak2->GetHPD())+HTilt*(((int)peak2->GetColumn())%2),peak2->GetMuX()+VSpace*(6-peak2->GetColumn()));
				if(thebox==1) g_final->SetPoint(g_final->GetN(),peak2->GetMuY()+32*peak2->GetHPD()+HTilt*((1+(int)peak2->GetColumn())%2),peak2->GetMuX()+VSpace*peak2->GetColumn());
			}
		}
// 		peaklist2->Clear("nodelete");
	}
	produceplot("Final", g_final, thebox);
	
}

void OnePeakPerStep(TList *hpds, int thebox){
	int i=0,k=0;
	
	cout<< "Flagging Starting peakmap for Likelihood estimation"<<endl;
	
	TGraph *g_final=new TGraph();
	double VSpace=30;
	double HTilt=16;
	int tmpstep=-1;
	int peakcounter=0;
	int currentpeak=0;
	double higherDLL=DLL_CUT-999;
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd2=(CleanHPD*)hpds->At(i);
// 		cout << hpd2->GetNPeaks() << endl;
// 		cout << hpd2->GetHitmap()->GetSize() <<endl;
		TList *peaklist2=(TList*)hpd2->GetHitmap();
		tmpstep=-1;
		higherDLL=DLL_CUT-999;
		peakcounter=0;
		currentpeak=0;
		for(k=0; k<hpd2->GetNPeaks(); k++){
			CleanPeak *peak2=(CleanPeak*)peaklist2->At(k);
			peakcounter++;
			if(peak2->IsSingle() && peak2->GetPassLevel()!="DIS_STEP1"){
// 				cout << peak2->GetCalibStep() << peak2->GetName() << endl;
				if(thebox==0) g_final->SetPoint(g_final->GetN(),peak2->GetMuY()+32*(13-peak2->GetHPD())+HTilt*(((int)peak2->GetColumn())%2),peak2->GetMuX()+VSpace*(6-peak2->GetColumn()));
				if(thebox==1) g_final->SetPoint(g_final->GetN(),peak2->GetMuY()+32*peak2->GetHPD()+HTilt*((1+(int)peak2->GetColumn())%2),peak2->GetMuX()+VSpace*peak2->GetColumn());
			}
		}
// 		peaklist2->Clear();
	}
	
// 	for(i=0; i<nCOLs*nHPDs; i++){
// 		CleanHPD *hpd2=(CleanHPD*)hpds->At(i);
// 		cout << hpd2->GetNPeaks() << endl;
// 		cout << hpd2->GetHitmap()->GetSize() <<endl;
// 	}
	
	produceplot("One_peak_per_step", g_final, thebox);
	
}

void PlotSingleHPDsPattern(CleanHPD *hpd){
	TGraph *g_initial=new TGraph();
	TGraph *g_single=new TGraph();
	TGraph *g_final=new TGraph();
	
	TList *peaklist=(TList*)hpd->GetHitmap();
	for(int i=0; i<hpd->GetNPeaks();i++){
		CleanPeak *peak=(CleanPeak*)peaklist->At(i);
		g_initial->SetPoint(g_initial->GetN(),peak->GetMuX(),peak->GetMuY());
		if(peak->IsSingle()){
			g_single->SetPoint(g_single->GetN(),peak->GetMuX(),peak->GetMuY());
		}
		if(peak->GetPassLevel()=="OKX2"){
			g_final->SetPoint(g_final->GetN(),peak->GetMuX(),peak->GetMuY());
		}
	}
	
	g_initial->SetName(Form("pattern_initial_COL%i_HPD%i",hpd->GetColumn(),hpd->GetHPD()));
	g_single->SetName(Form("pattern_single_COL%i_HPD%i",hpd->GetColumn(),hpd->GetHPD()));
	g_final->SetName(Form("pattern_final_COL%i_HPD%i",hpd->GetColumn(),hpd->GetHPD()));
	
	g_initial->Write();
	g_single->Write();
	g_final->Write();
}

void PlotHPDsPatterns(TList *hpds){
	for(int i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
		PlotSingleHPDsPattern(hpd);
	}
}

void AllPeaksPlot(TList *hpds, int thebox){
	int i=0,k=0;
	
	TGraph *g_final=new TGraph();
	double VSpace=30;
	double HTilt=16;
	int tmpstep=-1;
	int peakcounter=0;
	int currentpeak=0;
	double higherDLL=DLL_CUT-999;
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd2=(CleanHPD*)hpds->At(i);
// 		cout << hpd2->GetNPeaks() << endl;
// 		cout << hpd2->GetHitmap()->GetSize() <<endl;
		TList *peaklist2=(TList*)hpd2->GetHitmap();
		tmpstep=-1;
		higherDLL=DLL_CUT-999;
		peakcounter=0;
		currentpeak=0;
		for(k=0; k<hpd2->GetNPeaks(); k++){
			CleanPeak *peak2=(CleanPeak*)peaklist2->At(k);
			peakcounter++;
			if(true){
// 				cout << peak2->GetCalibStep() << peak2->GetName() << endl;
				if(thebox==0) g_final->SetPoint(g_final->GetN(),peak2->GetMuY()+32*(13-peak2->GetHPD())+HTilt*(((int)peak2->GetColumn())%2),peak2->GetMuX()+VSpace*(6-peak2->GetColumn()));
				if(thebox==1) g_final->SetPoint(g_final->GetN(),peak2->GetMuY()+32*peak2->GetHPD()+HTilt*((1+(int)peak2->GetColumn())%2),peak2->GetMuX()+VSpace*peak2->GetColumn());
			}
		}
// 		peaklist2->Clear();
	}
	
// 	for(i=0; i<nCOLs*nHPDs; i++){
// 		CleanHPD *hpd2=(CleanHPD*)hpds->At(i);
// 		cout << hpd2->GetNPeaks() << endl;
// 		cout << hpd2->GetHitmap()->GetSize() <<endl;
// 	}
	
	produceplot("All_Peaks", g_final, thebox);
	
}

void FinalTree(TList *hpds, int thebox){
	int i=0,k=0;
	
	cout<< "Flagging final peaks"<<endl;
	
	TGraph *g_final=new TGraph();
	double VSpace=30;
	double HTilt=16;
	int tmpstep=-1;
	int peakcounter=0;
	int currentpeak=0;
	double higherDLL=DLL_CUT-999;
	
	TNtuple *ntuple=new TNtuple("Tree","Tree","col:hpd:oldX:oldY:pattern:max:prop:xdiff:npes:sigma:DLL:is_single");
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd2=(CleanHPD*)hpds->At(i);
		TList *peaklist2=(TList*)hpd2->GetHitmap();
		tmpstep=-1;
		higherDLL=DLL_CUT-999;
		peakcounter=0;
		currentpeak=0;
		for(k=0; k<hpd2->GetNPeaks(); k++){
			CleanPeak *peak2=(CleanPeak*)peaklist2->At(k);
			peakcounter++;
			if(peak2->GetPassLevel()=="OKX2"){
				ntuple->Fill((double)peak2->GetColumn(),(double)peak2->GetHPD(),peak2->GetMuX(),peak2->GetMuY(),(double)peak2->GetCalibStep(),peak2->GetMax(),peak2->GetProp(),peak2->GetXDiff(),peak2->GetNPEs(),peak2->GetSigma(),peak2->GetLike3(),(double)peak2->IsSingle());
			}
		}
	}
// 	ntuple->Write();
}

void NPeaksPerHPD(TList *hpds){
	int i=0,k=0;
	TH1I *nPeaks_HPD_all=new TH1I("nPeaksPerHPD_all","nPeaksPerHPD_all;#peaks/HPD;#HPDs",1000,0,1000);
	TH1I *nPeaks_HPD_single=new TH1I("nPeaksPerHPD_single","nPeaksPerHPD_single;#peaks/HPD;#HPDs",1000,0,1000);
	TH1I *nPeaks_HPD_cleaned=new TH1I("nPeaksPerHPD_cleaned","nPeaksPerHPD_cleaned;#peaks/HPD;#HPDs",1000,0,1000);
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
// 		cout << hpd2->GetNPeaks() << endl;
// 		cout << hpd2->GetHitmap()->GetSize() <<endl;
		TList *peaklist=(TList*)hpd->GetHitmap();
		
		int peakcounter_all=0;
		int peakcounter_single=0;
		int peakcounter_cleaned=0;
		
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			peakcounter_all++;
			if(peak->IsSingle()) peakcounter_single++;
			if(peak->GetPassLevel()=="OKX2") peakcounter_cleaned++;
		}
		if(peakcounter_all>0) nPeaks_HPD_all->Fill(peakcounter_all);
		if(peakcounter_single>0) nPeaks_HPD_single->Fill(peakcounter_single);
		if(peakcounter_cleaned>0) nPeaks_HPD_cleaned->Fill(peakcounter_cleaned);
// 		peaklist2->Clear();
	}
	
	nPeaks_HPD_all->Write();
	nPeaks_HPD_single->Write();
	nPeaks_HPD_cleaned->Write();
}

void LightSpotFeatures(TList *hpds){
	int i=0,k=0;
	
	float QE=0.13;
	float chipE=1;
	
	TH1F *h_sphericity=new TH1F("sphericity","Sphericity;#epsilon;#peaks",100,0,1);
	TH1F *h_sphericity_ratio=new TH1F("sphericity_ratio","Sphericity;#sigma ratio;#peaks",150,0,2);
	TH1F *h_sphericity_ratio_norm=new TH1F("sphericity_ratio_norm","Sphericity;#sigma ratio;#peaks",150,0,1);
	TH1F *h_sphericity_diff=new TH1F("sphericity_diff","Sphericity;#sigma diff;#peaks",150,-2,2);
	TH1F *h_sphericity_diff_norm=new TH1F("sphericity_diff_norm","Sphericity;#sigma diff;#peaks",150,-2,2);
	TH1F *h_spot_ns=new TH1F("spot_ns","spot_ns",1000,0,1);
	TH1F *h_phyield_ns=new TH1F("phyield_ns","LED Photon Yield Distribution;# #gamma / 25ns;Events",1000,0,1/(QE*chipE));
	
	for(i=0; i<nCOLs*nHPDs; i++){
		CleanHPD *hpd=(CleanHPD*)hpds->At(i);
// 		cout << hpd2->GetNPeaks() << endl;
// 		cout << hpd2->GetHitmap()->GetSize() <<endl;
		TList *peaklist=(TList*)hpd->GetHitmap();
		
		for(k=0; k<hpd->GetNPeaks(); k++){
			CleanPeak *peak=(CleanPeak*)peaklist->At(k);
			if(peak->GetPassLevel()=="OKX2"){
				double a=peak->GetRMSX();
				double b=peak->GetRMSY();
				double epsilon=0;
				double ratio=0;
				double ratio_norm=0;
				double diff=0;
				double diff_norm=0;
				if(a>b){
					epsilon=TMath::Sqrt(1.-pow(b/a,2));
					ratio_norm=b/a;
				}
				if(a<b){
					epsilon=TMath::Sqrt(1.-pow(a/b,2));
					ratio_norm=a/b;
				}
				ratio=a/b;
				diff=a-b;
				diff_norm=(a-b)*2/(a+b);
				if(peak->IsSingle()){
					h_sphericity->Fill(epsilon);
					h_sphericity_ratio->Fill(ratio);
					h_sphericity_ratio_norm->Fill(ratio_norm);
					h_sphericity_diff->Fill(diff);
					h_sphericity_diff_norm->Fill(diff_norm);
					h_spot_ns->Fill((double)peak->GetNPEs()/6000.);
					h_phyield_ns->Fill((double)peak->GetNPEs()/(6000.*QE*chipE));
				}
				
			}
		}
// 		peaklist2->Clear();
	}
	
	h_sphericity->Write();
	h_sphericity_ratio->Write();
	h_sphericity_ratio_norm->Write();
	h_sphericity_diff->Write();
	h_sphericity_diff_norm->Write();
	h_spot_ns->Write();
	h_phyield_ns->Write();
}
