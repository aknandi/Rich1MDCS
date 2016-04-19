//HPD class for MDMS cleaning
//implementd by Andrea Contu 6 Sept 2010, acontu@cern.ch
#include "CleanHPD.h"
// #include "TMap.h"
// ClassImp(CleanHPD);
//Private method: m_fill_map
bool CleanHPD::m_fill_map(TTree *tree){
	cout << endl;
	cout << Form("Initialising HPD_%i on Column_%i, BOX %i",m_hpd,m_col,m_box) << endl;
	
	CleanPeak *peak;
	int DIM=1000;
	unsigned int i=0;
	int k=0;
	float *muX=new float[DIM];
	float *muY=new float[DIM];
	float *sigmaX=new float[DIM];
	float *sigmaY=new float[DIM];
	float *meanX=new float[DIM];
	float *meanY=new float[DIM];
	float *rmsX=new float[DIM];
	float *rmsY=new float[DIM];
	float *box=new float[DIM];
	float *col=new float[DIM];
	float *hpd=new float[DIM];
	float *npes=new float[DIM];
	float *peakclass=new float[DIM];
	float *proportion=new float[DIM];
	float *rotation=new float[DIM];
	float *chi2=new float[DIM];
	float *max=new float[DIM];
	int npeaks;
	int calibstep;
	
	int npeaks_this=0;
	
	
	tree->SetBranchAddress("CalibStep",&calibstep);
	tree->SetBranchAddress("max",max);
	tree->SetBranchAddress("muX",muX);
	tree->SetBranchAddress("muY",muY);
	tree->SetBranchAddress("sigmaX",sigmaX);
	tree->SetBranchAddress("sigmaY",sigmaY);
	tree->SetBranchAddress("meanX",meanX);
	tree->SetBranchAddress("meanY",meanY);
	tree->SetBranchAddress("rmsX",rmsX);
	tree->SetBranchAddress("rmsY",rmsY);
	tree->SetBranchAddress("HPDColumn",col);
	tree->SetBranchAddress("HPD",hpd);
	tree->SetBranchAddress("nPeaks",&npeaks);
	tree->SetBranchAddress("Box",box);
	tree->SetBranchAddress("peakClass",peakclass);
	tree->SetBranchAddress("nPEs",npes);
	tree->SetBranchAddress("proportion",proportion);
	tree->SetBranchAddress("rotation",rotation);
	tree->SetBranchAddress("chi2",chi2);
// 	TTreeFormula *formula=new TTreeFormula("cuts",cuts.Data(),tree);
	TString thehpd="";
	
	for(i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);
		for(k=0;k<npeaks;k++){
			if(((int)box[k])==m_box && ((int)hpd[k])==m_hpd && ((int)col[k])==m_col && max[k]>10 && npes[k]>200 && proportion[k]>0.1){
				TString boxstring="U";
// 				cout << max[k] << endl;
				if(box[k]==1) boxstring="D";
				thehpd=Form("peakmap_%07i%07i%02i%02i",(int)calibstep,i,m_hpd,m_col);
				peak=new CleanPeak(calibstep, box[k], col[k], hpd[k], npes[k], peakclass[k], max[k], proportion[k], meanX[k], meanY[k], rmsX[k], rmsY[k], sigmaX[k], sigmaY[k], muX[k], muY[k], rotation[k],chi2[k]);
				peak->SetName(thehpd);
// 				cout << thehpd << endl;
// 				cout << peak->GetCalibStep() << endl;
				npeaks_this++;
				m_hitmap.Add(peak);
			}
		}
	}
	m_hitmap.Sort();
	
	m_npeaks=npeaks_this;
	
	double meanx=0,meany=0;
	int tmpstep=-1;
	int countsteps=0;
	for(i=0; i<m_npeaks; i++){
		CleanPeak *tmp1=(CleanPeak*)m_hitmap.At(i);
// 		thehpd=Form("peakmap_%i02%i02%05i%07",tmp1->GetColumn(),tmp1->GetHPD(),tmp1->GetCalibStep(),i);
// 		tmp1->SetName(thehpd);
// 		cout << thehpd << endl;
		meanx+=tmp1->GetMuX();
		meany+=tmp1->GetMuY();
		if(tmp1->GetCalibStep()!=tmpstep){
			tmpstep=tmp1->GetCalibStep();
			if(countsteps==1){
				CleanPeak *tmp2=(CleanPeak*)m_hitmap.At(i-1);
				tmp2->SetSingle(true);
			}
			if(i==m_npeaks-1) tmp1->SetSingle(true);
			countsteps=0;
		}
		countsteps++;
	}
	
	if(m_npeaks>0) m_meanX=meanx/m_npeaks;
	if(m_npeaks>0) m_meanY=meany/m_npeaks;
	
// 	for(i=0; i<m_npeaks; i++){
// 		CleanPeak *tmp3=(CleanPeak*)m_hitmap.At(i);
// 		cout << tmp3->GetName() << "  "<<tmp3->GetCalibStep()<< "   " << tmp3->IsSingle() << endl;
// 	}
	
	tree->ResetBranchAddresses();
	
	delete muX;
	delete muY;
	delete sigmaX;
	delete sigmaY;
	delete meanX;
	delete meanY;
	delete rmsX;
	delete rmsY;
	delete box;
	delete col;
	delete hpd;
	delete npes;
	delete peakclass;
	delete proportion;
	delete rotation;
	delete max;
	
	cout << "Found "<< m_npeaks << " peaks."<<endl;
	cout << "Estimated HPD centre at X=" << this->GetMeanX() <<" Y="<< this->GetMeanY() <<endl;
	cout << "DONE" << endl;
	cout << endl;
	
	return true;
}

//DrawPeakMap
TGraph* CleanHPD::DrawPeakMap(){
	int i=0;
	double *xvec=new double[m_npeaks];
	double *yvec=new double[m_npeaks];
// 	cout << m_hitmap.GetSize() << endl;
	for(i=0;i<m_npeaks;i++){
		CleanPeak *peak=(CleanPeak*)m_hitmap.At(i);
		xvec[i]=peak->GetMuX();
		yvec[i]=peak->GetMuY();
	}
	
	TGraph *g_map=new TGraph(m_npeaks,xvec,yvec);
	g_map->GetXaxis()->SetTitle("x (LHCb pixels)");
	g_map->GetYaxis()->SetTitle("y (LHCb pixels)");
	g_map->GetXaxis()->SetLimits(0,32);
	g_map->GetYaxis()->SetLimits(0,32);
	g_map->SetMarkerStyle(2);
	g_map->SetMarkerSize(1);
	g_map->SetMarkerColor(4);
	g_map->SetTitle(Form("Peak Map - Box %i,Column %i, HPD %i",m_box,m_col,m_hpd));
	return g_map;
}


bool FitHPDCenter(){
	//fit circle
	return true;
}

//DrawPeakHisto
TH2D* CleanHPD::DrawPeakHisto(){
	int i=0;
	TH2D *h2=new TH2D(Form("HPD2D_HPD%i_Col%i_%i",this->GetHPD(),this->GetColumn(),this->GetBox()),Form("Peak Distribution - HPD%i COL%i BOX%i;x (LHCb pixels);y (LHCb pixels)",this->GetHPD(),this->GetColumn(),this->GetBox()),32,0,32,32,0,32);
	for(i=0;i<m_npeaks;i++){
		CleanPeak *peak=(CleanPeak*)m_hitmap.At(i);
		h2->Fill(peak->GetMuX(),peak->GetMuY());
	}
	
	return h2;
}
