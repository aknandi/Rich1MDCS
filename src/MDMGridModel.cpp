// -----------------------------------------------------------------------------
// Implementation file for class : MDMGridModel
//
// 2010-09-17 : Andrea Contu
// -----------------------------------------------------------------------------

#include "MDMGridModel.h"


//==============================================================================
//Default Constructor
//==============================================================================
MDMGridModel::MDMGridModel(){
	cX=-9999.;
	cY=-9999.;
	step=-9999.;
	phi=-9999.;
	tot_cols=0;
	tot_rows=0;
	minstep=999999;
	maxstep=-999999;
}

//==============================================================================
//Default Constructor
//==============================================================================
MDMGridModel::MDMGridModel(TString name, std::vector<Gaudi::XYZPoint > pts, std::vector<int> steps, int scalef, double angle, double width, int hpd, int col){
	cX=-9999.;
	cY=-9999.;
	step=-9999.;
	phi=-9999.;
	hpd_number=hpd;
	col_number=col;
	tot_cols=0;
	tot_rows=0;
	minstep=999999;
	maxstep=-999999;
	gname=name;
	//define grid parameters
	st=steps;
	points=pts;
	pattern_type.resize(st.size());
	pong_row.resize(st.size());
	pong_col.resize(st.size());
	modules_row_separation=2.5*step;
	modules_col_separation=step;
	scalefactor=scalef;
	this->set_rot(angle);
	this->set_width(width);
	
	this->Initialise();
	this->FlagSteps(points,steps);
	if(pvf) this->PreparePeaksAndGrids();
	
// 	this->FitPoints(points, st, scalefactor);
// 	this->FitPoints(points, st, scalefactor);
// 	this->FitPoints(points, st, scalefactor);
// 	this->FitPoints(points, st, scalefactor);
// 	this->FitPoints(points, st, scalefactor);
}

//==============================================================================
//Destructor
//==============================================================================
MDMGridModel::~MDMGridModel(){
	
}

//==============================================================================
// Initialise
//==============================================================================
void MDMGridModel::Initialise(){
	int theminstep=0;
	LOGIC = this->BuildMDMSLogic();
	int i=0;
	for(i=0;i<st.size();i++){
		if(minstep>st[i]){
			minstep=st[i];
			theminstep=i;
		}
		if(maxstep<st[i]){
			maxstep=st[i];
		}
	}
// 	this->PrintLOGIC(minstep,maxstep);
	PeaksCMS=this->peaksCMS(points);
}

//==============================================================================
// Initialise
//==============================================================================
Gaudi::XYZPoint MDMGridModel::peaksCMS(std::vector< Gaudi::XYZPoint > pts){
	unsigned int i=0;
	double meanx=0, meany=0;
	for(i=0;i<pts.size(); i++){
		meanx+=pts[i].X();
		meany+=pts[i].Y();
	}
	meanx/=pts.size();
	meany/=pts.size();
	
	Gaudi::XYZPoint tmppoint;
	tmppoint.SetXYZ(meanx,meany,0.0);
	return tmppoint;
}
//==============================================================================
//Flag step
//==============================================================================
void MDMGridModel::FlagSteps(std::vector< Gaudi::XYZPoint > datapoints, std::vector<int> steps){
	int i=0,k=0,l=0,m=0, n=0,times=5;
// 	pattern_type.clear();
	double yborder=0;
	
	double meanstep=0;
	
	for(i=0; i<steps.size();i++){
		meanstep+=(double)steps[i];
		if((datapoints[i].Y() - tan(phi)*datapoints[i].X())>yborder) pattern_type[i]="NONE";
		else pattern_type[i]="NONE";
	}
	
	if(steps.size()!=0) meanstep/=steps.size();
	int meanstepint=(int)meanstep;
	
	//check if meanstepint exisists, otherwise take the closest one
	bool exists=false;
	int closest=9999;
	int stepdist=9999;
	int itscol=-1, itsrow=-1;
	for(i=0; i<steps.size();i++){
		if(stepdist>abs(meanstepint-steps[i])){
			stepdist=abs(meanstepint-steps[i]);
			closest=steps[i];
			this->GetRowAndColumnFromLOGIC(steps[i], &itsrow, &itscol);
		}
		if(meanstepint==steps[i]){
			exists=true;
			this->GetRowAndColumnFromLOGIC(steps[i], &itsrow, &itscol);
			break;
		}
// 		std::cout << steps[i]<<std::endl;
	}
	if(!exists){
		meanstepint=closest;
	}
	
	std::cout << "MeanStepC = "<< closest << std::endl;
	std::cout << "MeanStepD = "<< meanstep << std::endl;
	std::cout << "MeanStep = "<< meanstepint << std::endl;
	
	
	//find good pivot
	int pivot_st=-1;
	bool pivot_found=false;
	int count_pivots=0;
	for(i=0; i<steps.size();i++){
// 		count_pivots=0
		int row_i=-1, col_i=-1;
		this->GetRowAndColumnFromLOGIC(steps[i], &row_i, &col_i);
		for(k=0; k<steps.size();k++){
			
			int row_k=-1, col_k=-1;
			this->GetRowAndColumnFromLOGIC(steps[k], &row_k, &col_k);
			double dist1=sqrt((datapoints[i].X()-datapoints[k].X())*(datapoints[i].X()-datapoints[k].X())+(datapoints[i].Y()-datapoints[k].Y())*(datapoints[i].Y()-datapoints[k].Y()));
			
			if(k!=i && fabs(dist1-step)<(step/2) && row_i==row_k && col_i==(col_k-1)){
				
				for(l=0; l<steps.size();l++){
					
					int row_l=-1, col_l=-1;
					this->GetRowAndColumnFromLOGIC(steps[l], &row_l, &col_l);
					double dist2=sqrt((datapoints[i].X()-datapoints[l].X())*(datapoints[i].X()-datapoints[l].X())+(datapoints[i].Y()-datapoints[l].Y())*(datapoints[i].Y()-datapoints[l].Y()));
					
					if(i!=l && l!=k && fabs(dist2-step)<(step/2) && row_i==row_l && col_i==(col_l+1)){
						
						for(m=0; m<steps.size();m++){
							
							int row_m=-1, col_m=-1;
							this->GetRowAndColumnFromLOGIC(steps[m], &row_m, &col_m);
							double dist3=sqrt((datapoints[i].X()-datapoints[m].X())*(datapoints[i].X()-datapoints[m].X())+(datapoints[i].Y()-datapoints[m].Y())*(datapoints[i].Y()-datapoints[m].Y()));
							
							if(m!=i && m!=k && m!=l && fabs(dist3-step)<(step/2) && row_i==(row_m+1) && col_i==col_m){
								
								for(n=0; n<steps.size();n++){
							
									int row_n=-1, col_n=-1;
									this->GetRowAndColumnFromLOGIC(steps[n], &row_n, &col_n);
									double dist4=sqrt((datapoints[i].X()-datapoints[n].X())*(datapoints[i].X()-datapoints[n].X())+(datapoints[i].Y()-datapoints[n].Y())*(datapoints[i].Y()-datapoints[n].Y()));
									
									if(n!=i && n!=k && n!=l && n!=m && fabs(dist4-step)<(step/2) && row_i==(row_n-1) && col_i==col_n){
										
										pivot_st=i;
										pivot_found=true;
										
									}
								}
								
							}
							if(pivot_found) break;
						}
					}
					if(pivot_found) break;
				}
			}
			if(pivot_found) break;
		}
		
		if(pivot_found) break;
	}
	
	if(pivot_found){
		std::cout << "Pivot found at step " <<steps[pivot_st]<< std::endl;
		pattern_type[pivot_st]="PIVOT";
		gridspivot=pivot_st;
// 		this->PrintGridsPivot();
// 		std::cout << "Pivot found at step " <<steps[pivot_st]<< std::endl;
	}
	else{
		std::cout << "Pivot not found! Check HPD!" << std::endl;
		gridspivot=-1;
	}
	
	pvf=pivot_found;
	//step eater
	int row_p=-1, col_p=-1;
	if(pvf) this->GetRowAndColumnFromLOGIC(steps[pivot_st], &row_p, &col_p);
	
	for(i=0; i<steps.size() && pivot_found;i++){
		int tmprow=-1, tmpcol=-1;
		this->GetRowAndColumnFromLOGIC(steps[i], &tmprow, &tmpcol);
		if(i!=pivot_st){
			pattern_type[i]="G1";
			double expectedY = datapoints[pivot_st].Y() + sin(phi)*step*(tmprow-row_p) - cos(phi)*step*(tmpcol-col_p);
			double expectedX = datapoints[pivot_st].X() + cos(phi)*step*(tmprow-row_p) + sin(phi)*step*(tmpcol-col_p);
			
// 			std::cout << tmprow << "\t"<<tmpcol<<std::endl;
// 			std::cout << row_p << "\t"<<col_p<<std::endl;
// 			std::cout << datapoints[i].X() << "\t"<<datapoints[i].Y()<<std::endl;
// 			std::cout << expectedX << "\t"<<expectedY<<std::endl;
// 			std::cout << datapoints[pivot_st].X() << "\t"<<datapoints[pivot_st].Y()<<std::endl;
			
// 			int expsign=(expectedY)/fabs(expectedY);
// 			int realsign=(datapoints[i].Y()-datapoints[pivot_st].Y())/fabs(datapoints[i].Y()-datapoints[pivot_st].Y()-tan(phi)*(datapoints[i].X()-datapoints[pivot_st].X()));
			double dist=sqrt((expectedX-datapoints[i].X())*(expectedX-datapoints[i].X())+(expectedY-datapoints[i].Y())*(expectedY-datapoints[i].Y()));
			
			if(dist<2*step) pattern_type[i]="G1";
			else pattern_type[i]="G2";
		}
	}
	
}

//==============================================================================
//Fit grid from TGraph
//==============================================================================
std::vector<double> MDMGridModel::FitPoints(std::vector< Gaudi::XYZPoint > datapoints, std::vector<int> steps, int scalef=14){
	std::vector<double> pars;
	unsigned int i=0;
	int theminstep=0;
// 	LOGIC = this->BuildMDMSLogic();
// 	int i=0;
	for(i=0;i<st.size();i++){
		if(minstep>st[i]){
			minstep=st[i];
			theminstep=i;
		}
		if(maxstep<st[i]){
			maxstep=st[i];
		}
	}
	
	if(cX==-9999. && cY==9999.){
		cX=(datapoints[theminstep]).X();
		cY=(datapoints[theminstep]).Y();
	}
	
	double meanx=0, meany=0;
	for(i=0;i<st.size();i++){
		Gaudi::XYZPoint tmppoint;
		tmppoint = this->GetClosestGridPoint(datapoints[i],i);
		meanx+=tmppoint.X()-datapoints[i].X();
		meany+=tmppoint.Y()-datapoints[i].Y();
	}
	meanx/=st.size();
	meany/=st.size();
	
	cX-=meanx;
	cY-=meany;
	
// 	std::cout << "Shift: "<< meanx<<"  "<<meany<<std::endl;
	
	tot_rows=(maxstep-minstep)/scalefactor;
	tot_cols=tot_rows;
	
	this->FlagSteps(datapoints,st);
	
	pars.push_back(cX);
	pars.push_back(cY);
	
	return pars;
}

//==============================================================================
//Draw Lines
//==============================================================================
void MDMGridModel::DrawLines(TString options, bool orientation, double Xmin, double Xmax, double Ymin, double Ymax){
	int i=0;
	TLine *l;
	TLine *l2;
	
	double phi_prime=phi;
	double cenX=cX,cenY=cY;
	if(phi<0 && phi>-TMath::PiOver2()) phi_prime+=TMath::PiOver2();
	if(phi>0 && phi>TMath::Pi()) phi_prime-=TMath::Pi();
	
	
	if(orientation){
		cenX=-cY;
		cenY=cX;
	}
	
	//find Y1max	
	double Y1max=0;
	if(fabs(phi_prime)!=TMath::PiOver2() && fabs(phi_prime)!=(TMath::PiOver2()+TMath::Pi())) Y1max=cenY+(Xmin-cenX)*tan(phi_prime);
	else Y1max=9999999;
	
	//find X1min
	double X1min=0;
	if(fabs(phi_prime)!=TMath::Pi() && fabs(phi_prime)!=0) X1min=cenX+(Ymin-cenY)/tan(phi_prime);
	else X1min=-9999999;
	
	std::vector<double> Y1vect, X1vect;
	std::vector< Gaudi::XYZPoint > Y1end, X1end;
	
	//find starting point for Y
	if(Y1max<Ymax && Y1max!=9999999){
		while(Y1max<Ymax){
			Y1max+=step/cos(phi_prime);
		}
		Y1max-=step/cos(phi_prime);
	}
	else{
		if(Y1max>Ymax && Y1max!=9999999){
			while(Y1max>Ymax){
				Y1max-=step/cos(phi_prime);
			}
		}
	}
	
// 	std::cout << Xmin << std::endl;
	
	//find starting point for X
	if(X1min<Xmin && X1min!=-9999999){
		while(X1min<Xmin){
			X1min+=step/sin(phi_prime);
// 			std::cout << X1min << std::endl;
		}
// 		std::cout << X1min << std::endl;
	}
	else{
		if(X1min>Xmin && X1min!=-9999999){
			while(X1min>Xmin){
				X1min-=step/sin(phi_prime);
			}
			X1min+=step/sin(phi_prime);
		}
	}	
	//fill Y1 intercepts
	double ytmp=Y1max;
	while(ytmp>=Ymin && Y1max!=9999999){
		Y1vect.push_back(ytmp);
		ytmp-=step/cos(phi_prime);
	}
	
	//fill X1 intercepts
	double xtmp=X1min;
	while(xtmp<=Xmax && X1min!=-9999999){
		X1vect.push_back(xtmp);
		
		xtmp+=step/sin(phi_prime);
	}
	
	//fill Y2 intercepts
	for(i=0; i<Y1vect.size() && Y1max!=9999999; i++){
		ytmp=Y1vect[i]+(Xmax-Xmin)*tan(phi_prime);
		xtmp=Xmax;
		if(ytmp>Ymax){
			ytmp=Ymax;
			xtmp=Xmin+(Ymax-Y1vect[i])/tan(phi_prime);
		}
		if(ytmp<Ymin){
			ytmp=Ymin;
			xtmp=Xmin+(Ymax-Y1vect[i])/tan(phi_prime);
		}
		Gaudi::XYZPoint ptmp;
		ptmp.SetXYZ(xtmp,ytmp,0.0);
		Y1end.push_back(ptmp);
	}
	
	//fill X2 intercepts
	for(i=0; i<X1vect.size() && X1min!=-9999999; i++){
		xtmp=X1vect[i]+(Ymax-Ymin)/tan(phi_prime);
		ytmp=Ymax;
		if(xtmp>Xmax){
			xtmp=Xmax;
			ytmp=Ymin+(Xmax-X1vect[i])*tan(phi_prime);
		}
		if(xtmp<Xmin){
			xtmp=Xmin;
			ytmp=Ymin+(Xmax-X1vect[i])*tan(phi_prime);
		}
		Gaudi::XYZPoint ptmp;
		ptmp.SetXYZ(xtmp,ytmp,0.0);
		X1end.push_back(ptmp);
	}
	
	
	//Draw Lines
	for(i=0; i<Y1vect.size();i++){
		if(!orientation) l = new TLine(Xmin,Y1vect[i],Y1end[i].X(),Y1end[i].Y());
		else l = new TLine(Y1vect[i],-Xmin,Y1end[i].Y(),-Y1end[i].X());
		l->SetLineStyle(2);
		l->Draw(options.Data());
	}
	
	for(i=0; i<X1vect.size();i++){
		if(!orientation) l = new TLine(X1vect[i],Ymin,X1end[i].X(),X1end[i].Y());
		else l = new TLine(Ymin,-X1vect[i],X1end[i].Y(),-X1end[i].X());
		l->SetLineStyle(2);
		l->Draw(options.Data());
	}
}

//==============================================================================
//Get Relative Point On the Grid
//==============================================================================
Gaudi::XYZPoint MDMGridModel::GetRelativePointOnGrid(int cstep, int cstep_ref, Gaudi::XYZPoint point){
	int col_ref=cstep_ref%scalefactor;
	int row_ref=cstep_ref/scalefactor;
	int col=cstep%scalefactor;
	int row=cstep/scalefactor;
	
	Gaudi::XYZPoint tmppoint;
	double xpoint=point.X()-step*(col-col_ref)*sin(phi)+step*(row-row_ref)*cos(phi);
	double ypoint=point.Y()+step*(col-col_ref)*cos(phi)+step*(row-row_ref)*sin(phi);
	tmppoint.SetXYZ(xpoint,ypoint,0.0);
	return tmppoint;
}

//==============================================================================
// Prepare Peaks and Grids for parameter fitting
//==============================================================================
void MDMGridModel::PreparePeaksAndGrids(){
	unsigned int i=0;
	
	int lpeak_1=99999;
	int lpeak_2=99999;
	
	std::vector<Gaudi::XYZPoint> points_1,points_2;
	
	for(i=0;i<st.size();i++){
		if(pattern_type[i]=="G1" || pattern_type[i]=="PIVOT"){
			points_1.push_back(points[i]);
			if(st[i]<lpeak_1) lpeak_1=st[i];
		}
		if(pattern_type[i]=="G2"){
			points_2.push_back(points[i]);
			if(st[i]<lpeak_2) lpeak_2=st[i];
		}
	}
	
	Gaudi::XYZPoint G1CMS=this->peaksCMS(points_1);
	Gaudi::XYZPoint G2CMS=this->peaksCMS(points_2);
	
	if(G1CMS.Y()>G2CMS.Y()) topgrid=1;
	else topgrid=2;
	
	int row_p=-1,col_p=-1;
	this->GetRowAndColumnFromLOGIC(st[gridspivot], &row_p, &col_p);
	
	for(i=0;i<st.size();i++){
		int tmprow=-1,tmpcol=-1;
		this->GetRowAndColumnFromLOGIC(st[i], &tmprow, &tmpcol);
		double expectedY=-999;
		double expectedX=-999;
		if(pattern_type[i]!="G2"){
			expectedY = points[gridspivot].Y() + sin(phi)*step*(tmprow-row_p) - cos(phi)*step*(tmpcol-col_p);
			expectedX = points[gridspivot].X() + cos(phi)*step*(tmprow-row_p) + sin(phi)*step*(tmpcol-col_p);
		}
		else{
			//69.7 is the difference in Y between 2 modules
			double module_diff=69.7;
			//last module is made of 2 led matrices instead of 2
			if(hpd_number==13 && (col_number%2)==0) module_diff/=2.0;
			
			if(topgrid==1){
				expectedY = points[gridspivot].Y() + sin(phi)*step*(tmprow-row_p) - cos(phi)*step*(tmpcol-col_p) - module_diff;
// 				expectedY = points[gridspivot].Y() + sin(phi)*step*(tmprow-row_p) - cos(phi)*step*(tmpcol-col_p);
				expectedX = points[gridspivot].X() + cos(phi)*step*(tmprow-row_p) + sin(phi)*step*(tmpcol-col_p);
			}
			if(topgrid==2){
				expectedY = points[gridspivot].Y() + sin(phi)*step*(tmprow-row_p) - cos(phi)*step*(tmpcol-col_p) + module_diff;
// 				expectedY = points[gridspivot].Y() + sin(phi)*step*(tmprow-row_p) - cos(phi)*step*(tmpcol-col_p);
				expectedX = points[gridspivot].X() + cos(phi)*step*(tmprow-row_p) + sin(phi)*step*(tmpcol-col_p);
			}
		}
		Gaudi::XYZPoint tmppeak;
		tmppeak.SetXYZ(expectedX,expectedY,0.0);
		pointsOnGrids.push_back(tmppeak);
	}
	
	
	cX_1=-9999;
	cX_2=-9999;
	cY_1=-9999;
	cY_2=-9999;
	lowestpeak_1=lpeak_1;
	lowestpeak_2=lpeak_2;
}

//==============================================================================
// PreFit
//==============================================================================
void MDMGridModel::PreFit(){
	
	gShift_pre.clear();
	gShift_pre_err.clear();
	parameters_pre.clear();
	parameters_pre_err.clear();
	
	gShift_pre.resize(2);
	gShift_pre_err.resize(2);
	parameters_pre.resize(10);
	parameters_pre_err.resize(10);
	
	unsigned int i=0;
	double meanx=0,meany=0;
	int countgoodpeaksfm=0;
	for(i=0;i<st.size();i++){
		if(sqrt((pointsOnGrids[i].X()-points[i].X())*(pointsOnGrids[i].X()-points[i].X())+(pointsOnGrids[i].Y()-points[i].Y())*(pointsOnGrids[i].Y()-points[i].Y()))<4*step){
			meanx+=(pointsOnGrids[i].X()-points[i].X());
			meany+=(pointsOnGrids[i].Y()-points[i].Y());
			countgoodpeaksfm++;
		}
	}
	meanx/=countgoodpeaksfm;
	meany/=countgoodpeaksfm;
// 	meanx=0;
// 	meany=0;
	gShift_pre[0]=(-meanx);
	gShift_pre[1]=(-meany);
	
	//fit good points only
	std::vector<bool> isgood(st.size());
	int countgood=0;
	for(i=0;i<st.size();i++){
		Gaudi::XYZPoint tmppoint;
		tmppoint.SetXYZ(pointsOnGrids[i].X()+gShift_pre[0],pointsOnGrids[i].Y()+gShift_pre[1],0.0);
		double r=points[i].Rho();
		double r_grid=tmppoint.Rho();
		
		double dr=r_grid-r;
		
// 		std::cout << dr << std::endl;
		
		double phi=points[i].Phi();
		double phi_grid=tmppoint.Phi();
		double dphi=phi_grid-phi;
		
		if(fabs(dr)<2*step && fabs(r*dphi)<2*step && fabs(dphi)<1 && fabs(dr)<5){
			countgood++;
			isgood[i]=true;
		}
		else{
			isgood[i]=false;
		}
	}
	
// 	std::cout << "CRASH HERE" << std::endl;
	
	TGraph *g_DRvsR=new TGraph(countgood);
	TGraph *g_DPhivsR=new TGraph(countgood);
	int goodcounter=0;
	
	double rmsx=0,rmsy=0;
	for(i=0;i<st.size();i++){
		Gaudi::XYZPoint tmppoint;
		tmppoint.SetXYZ(pointsOnGrids[i].X()+gShift_pre[0],pointsOnGrids[i].Y()+gShift_pre[1],0.0);
		
		double r=points[i].Rho();
		double r_grid=tmppoint.Rho();
		
		double dr=r_grid-r;
		
		double phi=points[i].Phi();
		
		double phi_grid=tmppoint.Phi();
		double dphi=phi_grid-phi;
		
		if(isgood[i]){
			g_DRvsR->SetPoint(goodcounter,r,dr);
			g_DPhivsR->SetPoint(goodcounter,r,dphi);
			goodcounter++;
		}
		
		rmsx+=((pointsOnGrids[i].X()-points[i].X())-meanx);
		rmsy+=((pointsOnGrids[i].Y()-points[i].Y())-meany);
	}
	rmsx/=st.size();
	rmsy/=st.size();
	rmsx=sqrt(rmsx);
	rmsy=sqrt(rmsy);
	gShift_pre_err[0]=(rmsx);
	gShift_pre_err[1]=(rmsy);
	
	TF1 *fr=new TF1("fr","pol3(0)",0,20);
	fr->FixParameter(0,0.0);
	fr->SetParameter(1,0.0);
	fr->SetParameter(2,0.0);
	fr->FixParameter(3,0.0);
	TF1 *fphi=new TF1("fphi","pol3(0)",0,20);
	fphi->SetParameter(0,0.0);
	fphi->SetParameter(1,0.0);
	fphi->SetParameter(2,0.0);
	fphi->FixParameter(3,0.0);
	
	TCanvas *tmpcan=new TCanvas(Form("PreFit_HPD%i_COL%i",hpd_number,col_number),Form("PreFit_HPD%i_COL%i",hpd_number,col_number));
	tmpcan->Divide(1,2);
	tmpcan->cd(1);
	g_DRvsR->SetMarkerStyle(8);
	g_DRvsR->Draw("AP");
	g_DRvsR->Fit(fr,"MQ");
	tmpcan->cd(2);
	g_DPhivsR->SetMarkerStyle(8);
	g_DPhivsR->Draw("AP");
	g_DPhivsR->Fit(fphi,"MQ");
	tmpcan->Write();
	tmpcan->Close();
	parameters_pre[0]=(fr->GetParameter(0));
	parameters_pre[1]=(fr->GetParameter(1));
	parameters_pre[2]=(fr->GetParameter(2));
	parameters_pre[3]=(fr->GetParameter(3));
	
	parameters_pre[4]=(fphi->GetParameter(0));
	parameters_pre[5]=(fphi->GetParameter(1));
	parameters_pre[6]=(fphi->GetParameter(2));
	parameters_pre[7]=(fphi->GetParameter(3));
	
	parameters_pre_err[0]=(fr->GetParError(0));
	parameters_pre_err[1]=(fr->GetParError(1));
	parameters_pre_err[2]=(fr->GetParError(2));
	parameters_pre_err[3]=(fr->GetParError(3));
	
	parameters_pre_err[4]=(fphi->GetParError(0));
	parameters_pre_err[5]=(fphi->GetParError(1));
	parameters_pre_err[6]=(fphi->GetParError(2));
	parameters_pre_err[7]=(fphi->GetParError(3));
}

//==============================================================================
// PreFit
//==============================================================================
void MDMGridModel::IteratedFit(){
	
	std::vector<Gaudi::XYZPoint> re_eval_p=points;
	std::vector<Gaudi::XYZPoint> re_evalR_p=points;
	std::vector<Gaudi::XYZPoint> re_evalPhi_p=points;
	unsigned int i=0;
	for(i=0;i<st.size();i++){
		double tmpr=sqrt(points[i].X()*points[i].X()+points[i].Y()*points[i].Y());
// 		std::cout << "TEST" <<std::endl;
// 		std::cout << tmpr << std::endl;
// 		std::cout << points[i].Rho() << std::endl;
		double dr=parameters[0]+parameters[1]*tmpr+parameters[2]*tmpr*tmpr+parameters[3]*tmpr*tmpr*tmpr;
		double dphi=parameters[4]+parameters[5]*tmpr+parameters[6]*tmpr*tmpr+parameters[7]*tmpr*tmpr*tmpr;
		
		double signX=points[i].X()/fabs(points[i].X());
		double signY=points[i].X()/fabs(points[i].Y());
		
		re_eval_p[i].SetXYZ(points[i].X()*(1.0+dr/tmpr), points[i].Y()*(1.0+dr/tmpr), 0.0);
		re_eval_p[i].SetXYZ((tmpr+dr)*cos(points[i].Phi()+dphi),(tmpr+dr)*sin(points[i].Phi()+dphi),0.0);
		re_evalR_p[i].SetXYZ(points[i].X(), points[i].Y(), 0.0);
		re_evalR_p[i].SetXYZ((tmpr)*cos(points[i].Phi()+dphi),(tmpr)*sin(points[i].Phi()+dphi),0.0);
		re_evalPhi_p[i].SetXYZ(points[i].X()*(1.0+dr/tmpr), points[i].Y()*(1.0+dr/tmpr), 0.0);
		re_evalPhi_p[i].SetXYZ((tmpr+dr)*cos(points[i].Phi()),(tmpr+dr)*sin(points[i].Phi()),0.0);
	}
	
	double meanx=0,meany=0;
	int countgoodpeaksfm=0;
	for(i=0;i<st.size();i++){
		if(sqrt((pointsOnGrids[i].X()-re_eval_p[i].X())*(pointsOnGrids[i].X()-re_eval_p[i].X())+(pointsOnGrids[i].Y()-re_eval_p[i].Y())*(pointsOnGrids[i].Y()-re_eval_p[i].Y()))<3*step){
			meanx+=(pointsOnGrids[i].X()-re_eval_p[i].X());
			meany+=(pointsOnGrids[i].Y()-re_eval_p[i].Y());
			countgoodpeaksfm++;
		}
	}
	meanx/=countgoodpeaksfm;
	meany/=countgoodpeaksfm;
	gShift[0]=(-meanx);
	gShift[1]=(-meany);
	
	
	//fit good points only
	std::vector<bool> isgood(st.size());
	int countgood=0;
	for(i=0;i<st.size();i++){
		double r=points[i].Rho();
		
		Gaudi::XYZPoint tmppoint;
		tmppoint.SetXYZ(pointsOnGrids[i].X()+gShift[0],pointsOnGrids[i].Y()+gShift[1],0.0);
		double r_grid=tmppoint.Rho();
		
		double dr=r_grid-r;
		
		double phi=points[i].Phi();
		double phi_grid=tmppoint.Phi();
		double dphi=phi_grid-phi;
		
		if(fabs(dr)<2*step && fabs(r*dphi)<2*step && fabs(dphi)<2 && fabs(dr)<5){
			countgood++;
			isgood[i]=true;
		}
		else{
			isgood[i]=false;
		}
	}
	
	TGraph *g_DRvsR=new TGraph(countgood);
	TGraph *g_DPhivsR=new TGraph(countgood);
	int goodcounter=0;
	
	double rmsx=0,rmsy=0;
	for(i=0;i<st.size();i++){
		double r=points[i].Rho();
		Gaudi::XYZPoint tmppoint;
		tmppoint.SetXYZ(pointsOnGrids[i].X()+gShift[0],pointsOnGrids[i].Y()+gShift[1],0.0);
		
		double r_grid=tmppoint.Rho();
		
		double dr=r_grid-r;
		
		double phi=points[i].Phi();
		double phi_grid=tmppoint.Phi();
		double dphi=phi_grid-phi;
		
		if(isgood[i]){
			g_DRvsR->SetPoint(goodcounter,r,dr);
			g_DPhivsR->SetPoint(goodcounter,re_evalPhi_p[i].Rho(),dphi);
			goodcounter++;
		}
		rmsx+=((pointsOnGrids[i].X()-points[i].X())-meanx);
		rmsy+=((pointsOnGrids[i].Y()-points[i].Y())-meany);
	}
	rmsx/=st.size();
	rmsy/=st.size();
	rmsx=sqrt(rmsx);
	rmsy=sqrt(rmsy);
	gShift_err[0]=(rmsx);
	gShift_err[1]=(rmsy);
	
	TF1 *fr=new TF1("fr","pol3(0)",0,20);
	fr->FixParameter(0,parameters[0]);
	fr->SetParameter(1,parameters[1]);
	fr->SetParameter(2,parameters[2]);
	fr->FixParameter(3,parameters[3]);
	TF1 *fphi=new TF1("fphi","pol3(0)",0,20);
	fphi->SetParameter(0,parameters[4]);
	fphi->SetParameter(1,parameters[5]);
	fphi->SetParameter(2,parameters[6]);
	fphi->FixParameter(3,parameters[7]);
	
	TCanvas *tmpcan=new TCanvas(Form("Fit_HPD%i_COL%i",hpd_number,col_number),Form("Fit_HPD%i_COL%i",hpd_number,col_number));
	tmpcan->Divide(1,2);
	tmpcan->cd(1);
	g_DRvsR->SetMarkerStyle(8);
	g_DRvsR->Draw("AP");
	g_DRvsR->Fit(fr,"MQ");
	tmpcan->cd(2);
	g_DPhivsR->SetMarkerStyle(8);
	g_DPhivsR->Draw("AP");
	g_DPhivsR->Fit(fphi,"MQ");
	tmpcan->Close();
	
	parameters[0]=(fr->GetParameter(0));
	parameters[1]=(fr->GetParameter(1));
	parameters[2]=(fr->GetParameter(2));
	parameters[3]=(fr->GetParameter(3));
	
	parameters[4]=(fphi->GetParameter(0));
	parameters[5]=(fphi->GetParameter(1));
	parameters[6]=(fphi->GetParameter(2));
	parameters[7]=(fphi->GetParameter(3));
	
	parameters_err[0]=(fr->GetParError(0));
	parameters_err[1]=(fr->GetParError(1));
	parameters_err[2]=(fr->GetParError(2));
	parameters_err[3]=(fr->GetParError(3));
	
	parameters_err[4]=(fphi->GetParError(0));
	parameters_err[5]=(fphi->GetParError(1));
	parameters_err[6]=(fphi->GetParError(2));
	parameters_err[7]=(fphi->GetParError(3));

	//Add the grid shifts in
	//x shift
	parameters[8] = gShift[0];
	parameters_err[8] = gShift_err[0];
	//y shift
	parameters[9] = gShift[1];
	parameters_err[9] = gShift_err[1];
	
	delete fr;
	delete fphi;
	delete g_DRvsR;
	delete g_DPhivsR;
	delete tmpcan;
}


//==============================================================================
// Fit parameters
//==============================================================================
bool MDMGridModel::FitParameters(){
	
	std::cout << "Initialise parameters"<<std::endl;
	gShift.clear();
	gShift_err.clear();
	parameters.clear();
	parameters_err.clear();
	
	gShift.resize(2);
	gShift_err.resize(2);
	parameters.resize(10);
	parameters_err.resize(10);
	
	
	if(pvf){
		this->PreFit();
		
		gShift=gShift_pre;
		gShift_err=gShift_pre_err;
		parameters=parameters_pre;
		parameters_err=parameters_pre_err;
		
		std::cout << "Fitting..."<<std::endl;
		
		int i=0;
		
		for(i=0;i<50;i++) this->IteratedFit();
	}
	else {
		std::cout << "Too Few or Bad peaks. Fit aborted!"<<std::endl;
		std::cout << "Setting parameters to 0"<<std::endl;
		gShift[0]=0;
		gShift[1]=0;
		gShift_err[0]=0;
		gShift_err[1]=0;
		parameters[0]=0;
		parameters[1]=0;
		parameters[2]=0;
		parameters[3]=0;
		parameters[4]=0;
		parameters[5]=0;
		parameters[6]=0;
		parameters[7]=0;
		parameters_err[0]=0;
		parameters_err[1]=0;
		parameters_err[2]=0;
		parameters_err[3]=0;
		parameters_err[4]=0;
		parameters_err[5]=0;
		parameters_err[6]=0;
		parameters_err[7]=0;
		std::cout << "DONE"<<std::endl;
	}
	return pvf;
	
}

//==============================================================================
//Get Closest Grid Point
//==============================================================================
Gaudi::XYZPoint MDMGridModel::GetClosestGridPoint(Gaudi::XYZPoint point, int k){
	Gaudi::XYZPoint p;
	
	int Xsign=((point.X()-cX) > 0 ) - ((point.X()-cX) < 0 );
	int Ysign=((point.Y()-cY) > 0 ) - ((point.Y()-cY) < 0 );
	
	double dist=sqrt((point.X()-cX)*(point.X()-cX) + (point.Y()-cY)*(point.Y()-cY));
	
	double angle_wrt_local=0;
	if((point.X()-cX)!=0) angle_wrt_local=asin((point.Y()-cY)/dist);
	
	double distlocalY=dist*fabs(sin(angle_wrt_local-phi));
	double distlocalX=dist*fabs(cos(angle_wrt_local-phi));
	
	int distinstepX=0, distinstepY=0;
	
	double tmpdist=0;
	while(tmpdist<distlocalY){
		tmpdist+=step;
		distinstepY++;
	}
	distinstepY--;
	tmpdist=0;
	while(tmpdist<distlocalX){
		tmpdist+=step;
		distinstepX++;
	}
	distinstepX--;
	
// 	std::cout << distinstepX*step <<"\t"<< distlocalX << std::endl;
// 	std::cout << distinstepY*step  <<"\t"<< distlocalY << std::endl;
	
	
	std::vector< Gaudi::XYZPoint > pivots;
	std::vector< int > pgr;
	std::vector< int > pgc;
	
	int col=0, row=0;
	
	int square=10;
	
	for(row=-square; row<=square; row++){
		for(col=-square; col<=square; col++){
			Gaudi::XYZPoint tmppivot;
			tmppivot.SetXYZ(cX+step*(distinstepX+row)*fabs(cos(phi))*Xsign - step*(distinstepY+col)*fabs(sin(phi))*Ysign,cY+step*(distinstepX+row)*fabs(sin(phi))*Xsign + step*(distinstepY+col)*fabs(cos(phi))*Ysign,0.0);
			pivots.push_back(tmppivot);
			pgr.push_back(distinstepX+row);
			pgc.push_back(distinstepY+col);
		}
	}
	
	int i=0;
	
	double pivox=0,pivoy=0;
	
	double fdist=9999999;
	for(i=0;i<pivots.size();i++){
		if(((pivots[i].X()-point.X())*(pivots[i].X()-point.X())+(pivots[i].Y()-point.Y())*(pivots[i].Y()-point.Y()))<fdist){
			fdist=((pivots[i].X()-point.X())*(pivots[i].X()-point.X())+(pivots[i].Y()-point.Y())*(pivots[i].Y()-point.Y()));
			pivox=pivots[i].X();
			pivoy=pivots[i].Y();
			pong_row[k]=pgr[i];
			pong_col[k]=pgc[i];
		}
	}
	
	p.SetXYZ(pivox,pivoy,0.0);
	
	return p;
}

//==============================================================================
//Draw Grid
//==============================================================================
void MDMGridModel::DrawGrid(TString options, std::vector<Gaudi::XYZPoint> points, double Xmin, double Xmax, double Ymin, double Ymax){
	int i=0,k=0;
// 	this->DrawLines(options, false, Xmin, Xmax, Ymin, Ymax);
// 	this->DrawLines(options, true, Xmin, Xmax, Ymin, Ymax);
	
// 	std::vector< TH1D* > h_res;
// 	h_res.push_back(new TH1D("tmphX","tmphX",100,-30,30));
// 	h_res.push_back(new TH1D("tmphY","tmphY",100,-30,30));
	
// 	TMarker *mark=new TMarker(cX,cY,8);
// 	mark->SetMarkerColor(2);
// 	mark->Draw(options.Data());
	
	
	int tmprow=-99999,tmpcol=-99999;
	for(i=0;i<st.size() && pvf;i++){
		
		TLatex latex;
		latex.SetTextSize(0.012);
		latex.SetTextAlign(13);//align at top
		latex.SetTextColor(2);
// 		if(pattern_type[i]=="G1") latex.SetTextColor(1);
// 		if(pattern_type[i]=="G2") latex.SetTextColor(14);
// 		if(pattern_type[i]=="PIVOT") latex.SetTextColor(4);
		latex.DrawLatex(points[i].X(),points[i].Y(),Form("%i",st[i]));
		TLatex latex2;
		latex2.SetTextSize(0.012);
		latex2.SetTextAlign(13);//align at top
		latex2.SetTextColor(1);
		latex2.DrawLatex(pointsOnGrids[i].X()+gShift[0]-2,pointsOnGrids[i].Y()+gShift[1],Form("%i",st[i]));
	}
	
	
	for(i=0;i<st.size() && pvf;i++){
		for(k=0;k<st.size();k++){
			if(i!=k && (pattern_type[i]==pattern_type[k] || (pattern_type[i]=="G1" && pattern_type[k]=="PIVOT") ||(pattern_type[k]=="G1" && pattern_type[i]=="PIVOT"))){
				int tmprow1=-1, tmpcol1=-1;
				this->GetRowAndColumnFromLOGIC(st[i], &tmprow1, &tmpcol1);
				int tmprow2=-1, tmpcol2=-1;
				this->GetRowAndColumnFromLOGIC(st[k], &tmprow2, &tmpcol2);
				if((abs(tmprow1-tmprow2)==1 && abs(tmpcol1-tmpcol2)==0) || (abs(tmprow1-tmprow2)==0 && abs(tmpcol1-tmpcol2)==1) || (abs(tmprow1-tmprow2)==2 && abs(tmpcol1-tmpcol2)==0) || (abs(tmprow1-tmprow2)==0 && abs(tmpcol1-tmpcol2)==2)){
					if(pattern_type[i]=="G1" || pattern_type[i]=="PIVOT"){
						TLine *l=new TLine(pointsOnGrids[i].X()+gShift[0],pointsOnGrids[i].Y()+gShift[1],pointsOnGrids[k].X()+gShift[0],pointsOnGrids[k].Y()+gShift[1]);
						l->SetLineColor(1);
						l->Draw(options.Data());
					}
					if(pattern_type[i]=="G2"){
						TLine *l=new TLine(pointsOnGrids[i].X()+gShift[0],pointsOnGrids[i].Y()+gShift[1],pointsOnGrids[k].X()+gShift[0],pointsOnGrids[k].Y()+gShift[1]);
						l->SetLineColor(14);
						l->Draw(options.Data());
					}
				}
			}
		}	
		
	}
	
}

std::vector<int> SumToVec(std::vector<int> rootvec, int number){
	int i=0;
	std::vector<int> copy=rootvec;
	for(i=0;i<rootvec.size();i++){
		copy[i]+=number;
	}
	return copy;
}

//==============================================================================
// Build MDMS Logic
//==============================================================================
std::vector< std::vector<int> > MDMGridModel::BuildMDMSLogic(){
	std::vector< std::vector<int> > outputmatrix;
	
	int i=0;
	std::vector<int> row_root(scalefactor);
	for(i=0;i<scalefactor;i++){
		row_root[i]=scalefactor-i-1;
	}
	
// 	LOGIC.clear();
	LOGIC_Shift.clear();
	
	int shift=1;
	int angle_shift=0;
	
	for(i=0;i<50;i++){
		if((i+1)%2==1){
			//position #odd
			shift+=3;
			outputmatrix.push_back( SumToVec(row_root,shift) );
			LOGIC_Shift.push_back( angle_shift );
			shift+=15;
			outputmatrix.push_back( SumToVec(row_root,shift) );
			LOGIC_Shift.push_back( angle_shift );
			shift+=15;
			outputmatrix.push_back( SumToVec(row_root,shift) );
			LOGIC_Shift.push_back( angle_shift );
			shift+=15;
		}
		else{
			//position #even
			shift+=3;
			outputmatrix.push_back( SumToVec(row_root,shift) );
			LOGIC_Shift.push_back( angle_shift );
			shift+=15;
			outputmatrix.push_back( SumToVec(row_root,shift) );
			LOGIC_Shift.push_back( angle_shift );
			shift+=30;
			angle_shift++;
		}
	}
	return outputmatrix;
}

//==============================================================================
// Get row and column
//==============================================================================
void MDMGridModel::GetRowAndColumnFromLOGIC(int pattern, int *row, int *col){
	int i=0,k=0;
	
	int tmprow=-999;
	int tmpcol=-999;
	for(i=0;i<LOGIC.size();i++){
		for(k=0;k<LOGIC[i].size();k++){
			if(LOGIC[i][k]==pattern){
				tmprow=i;
				tmpcol=k+LOGIC_Shift[i];
				break;
			}
		}
	}
	
	*col=tmpcol;
	*row=tmprow;
}

//==============================================================================
// Get Pivot Local Expected Coordinate
//==============================================================================
Gaudi::XYZPoint MDMGridModel::GetPivotExpectedCoord(){
	Gaudi::XYZPoint tmppoint; 
	tmppoint.SetXYZ(pointsOnGrids[gridspivot].X()+gShift[0],pointsOnGrids[gridspivot].Y()+gShift[1],0.0);
	return tmppoint;
}

//==============================================================================
// Get Pivot Local Coordinate
//==============================================================================
Gaudi::XYZPoint MDMGridModel::GetPivotCoord(){
	Gaudi::XYZPoint tmppoint; 
	tmppoint.SetXYZ(pointsOnGrids[gridspivot].X()+gShift[0],pointsOnGrids[gridspivot].Y()+gShift[1],0.0);
	return tmppoint;
}

//==============================================================================
// Print MDMS Logic
//==============================================================================
void MDMGridModel::PrintLOGIC(int minstep, int maxstep){
	int i=0,k=0;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << "       MDMS Pattern Logic (same for all lightbars)      " << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << std::endl;
	int printshift=0;
	int minshift=1;
	
	int minrow=-1;
	int maxrow=-1;
	
	if(minstep==-1){
		minstep=0;
		maxstep=2200;
	}
	for(i=(LOGIC.size()-1);i>=0;i--){
		
		for(k=0;k<LOGIC[i].size();k++){
			if(LOGIC[i][k]>=minstep && LOGIC[i][k]<=maxstep){
				if(LOGIC[i][k]==minstep){
					minshift=LOGIC_Shift[i];
					minrow=i;
// 					std::cout << LOGIC_Shift[i] << std::endl;
				}
				if(LOGIC[i][k]==maxstep){
					maxrow=i;
				}
			}
		}
			
		
	}
	if(minstep>=maxstep){
		minstep=0;
		maxstep=2200;
	}
	
	std::cout << "From step: "<<minstep<< " to step "<< maxstep << std::endl;
	std::cout << "Minimum column shift "<<minshift<<std::endl;
	
	for(i=(LOGIC.size()-1);i>=0;i--){
		if(i>=minrow && i<=maxrow){
			std::cout << std::endl;
			std::cout << std::endl;
			for(k=0;k<LOGIC_Shift[i]+printshift-minshift+1;k++){
				std::cout << "\t";
			}
			for(k=0;k<LOGIC[i].size();k++){
				if(LOGIC[i][k]>=minstep && LOGIC[i][k]<=maxstep ) std::cout << "\t" << Form("%04i",LOGIC[i][k]);
				else {
					if(LOGIC[i][k]>=minstep-50 && LOGIC[i][k]<=maxstep+50 ) std::cout << "\t" << "xxxx";
				}
			}
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "********************************************************" << std::endl;
}

//==============================================================================
// Get Resolution
//==============================================================================

bool MDMGridModel::GetResolution(std::vector< Gaudi::XYZPoint > datapoints, TH1D *h, TString Axis){
	int i=0;
	if(pvf){
		if(Axis=="X"){
			for(i=0; i<datapoints.size(); i++){
				h->Fill(datapoints[i].X()-pointsOnGrids[i].X()-gShift[0]);
			}
		}
		if(Axis=="Y"){
			for(i=0; i<datapoints.size(); i++){
				h->Fill(datapoints[i].Y()-pointsOnGrids[i].Y()-gShift[1]);
			}
		}
		if(Axis=="D"){
			for(i=0; i<datapoints.size(); i++){
				h->Fill(sqrt((datapoints[i].X()-pointsOnGrids[i].X()-gShift[0])*(datapoints[i].X()-pointsOnGrids[i].X()-gShift[0])+(datapoints[i].Y()-pointsOnGrids[i].Y()-gShift[1])*(datapoints[i].Y()-pointsOnGrids[i].Y()-gShift[1])));
			}
		}
	}
	return true;
};

//==============================================================================
// Get Intrinsic resolution
//==============================================================================
bool MDMGridModel::GetSelfResolution(TH1D *h, TString Axis){
	int i=0;
	
	std::vector<Gaudi::XYZPoint> datapoints(st.size());
	
	if(pvf){
		for(i=0;i<st.size();i++){
			double tmpr=sqrt(points[i].X()*points[i].X()+points[i].Y()*points[i].Y());
	// 		std::cout << "TEST" <<std::endl;
	// 		std::cout << tmpr << std::endl;
	// 		std::cout << points[i].Rho() << std::endl;
			double dr=parameters[0]+parameters[1]*tmpr+parameters[2]*tmpr*tmpr+parameters[3]*tmpr*tmpr*tmpr;
			double dphi=parameters[4]+parameters[5]*tmpr+parameters[6]*tmpr*tmpr+parameters[7]*tmpr*tmpr*tmpr;
			
			double signX=points[i].X()/fabs(points[i].X());
			double signY=points[i].X()/fabs(points[i].Y());
			
			datapoints[i].SetXYZ(points[i].X()*(1.0+dr/tmpr), points[i].Y()*(1.0+dr/tmpr), 0.0);
			datapoints[i].SetXYZ((tmpr+dr)*cos(points[i].Phi()+dphi),(tmpr+dr)*sin(points[i].Phi()+dphi),0.0);
		}
		
		
		if(Axis=="X"){
			for(i=0; i<datapoints.size(); i++){
				h->Fill(datapoints[i].X()-pointsOnGrids[i].X()-gShift[0]);
			}
		}
		if(Axis=="Y"){
			for(i=0; i<datapoints.size(); i++){
				h->Fill(datapoints[i].Y()-pointsOnGrids[i].Y()-gShift[1]);
			}
		}
		if(Axis=="D"){
			for(i=0; i<datapoints.size(); i++){
				h->Fill(sqrt((datapoints[i].X()-pointsOnGrids[i].X()-gShift[0])*(datapoints[i].X()-pointsOnGrids[i].X()-gShift[0])+(datapoints[i].Y()-pointsOnGrids[i].Y()-gShift[1])*(datapoints[i].Y()-pointsOnGrids[i].Y()-gShift[1])));
			}
		}
	}
	return true;
};
//=============================================================================
