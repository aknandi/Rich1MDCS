// -----------------------------------------------------------------------------
// Implementation file for class : MDMFitParameters
//
// 2010-09-17 : Andrea Contu
// -----------------------------------------------------------------------------

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"
#include "Kernel/RichSide.h"

//local
#include "MDMFitParameters.h"

//ROOT function
#include "TRandom3.h"

// -----------------------------------------------------------------------------

using namespace Rich;
// using namespace RooFit;

DECLARE_ALGORITHM_FACTORY(MDMFitParameters);

//==============================================================================
//Constructor
//==============================================================================
MDMFitParameters::MDMFitParameters(const std::string& name, ISvcLocator* pSvcLocator)
	: Rich::AlgBase ( name, pSvcLocator ),
	m_nEvt(0),
	m_npeaks(0),
	m_npeaks_ref(0),
	m_updMgrSvc(NULL),
	m_magFieldSvc(NULL),
	m_mag_cond(NULL)
{
	declareProperty("SetCleanedNTuple", m_infilepath = "");
	declareProperty("SetRunNumber", m_run = 0);
	declareProperty("SetTree", m_treepath = "Tree");
	declareProperty("Output_MDMS_DB_Xml", m_outfilepath = "out.xml");
	declareProperty("OutputTuple", m_outfiletuple = "MDMSOutput.root");
	declareProperty("Iteration", m_iteration = 0 );
	declareProperty("NHPDs", m_nhpds = 14);
	declareProperty("NColumns", m_ncols = 7);
	declareProperty("Panel", m_panel = 0);
	declareProperty("BField", m_bfield = 0);
	declareProperty("DefaultLinearMagnification", m_defmag = 1./0.18);
	declareProperty("DefaultLinearMagnificationO2", m_defmag_o2 = 0.0);
	declareProperty("DefaultLinearMagnificationO3", m_defmag_o3 = 0.0);
	declareProperty("ReferenceHPD", m_reference_HPD = 49);
	declareProperty("GridSize", m_gstep = 5.08);
	//Set the random seed to be 1,2,3 to make sure the U and D box fits are seeded with the same random values each time
	//Otherwise we can't guarantee that all the HPDs get the same radius change applied in a given fit number
	//Set this value in the python script as a command line argument 
	declareProperty("SetRandomSeed", m_randomseed = 0);
	//Set the radius variation factor for systematics
	declareProperty("SetRadiusFactor", m_radiusfactor = 1.0);
}

//==============================================================================
//Destructor
//==============================================================================
MDMFitParameters::~MDMFitParameters() {
	if(0 != m_magFieldSvc) m_magFieldSvc->release();
	if(0 != m_updMgrSvc) m_updMgrSvc->release();
}

//==============================================================================
//Initialize
//==============================================================================
StatusCode MDMFitParameters::initialize() {
	
	StatusCode sc = Rich::AlgBase::initialize();
	if( sc.isFailure() ) return sc;
	
	sc = sc && this->getMagFieldSvc();
	
	info()<< "Initialising MDMFitParameters..."<<endmsg;
	
// 	SmartDataPtr<TabulatedProperty> HPDdeMag(dataSvc(),"/dd/Materials/RichMaterialTabProperties/HpdDemagnification");
	
	info()<< "Request to process " << m_nhpds <<" HPDs and " << m_ncols << " Columns, Panel: " << m_panel << " (top=0, bottom=1) " << endmsg;
	
	//Fill arrays
	sc = sc && this->ChangeMagneticField(0,0.1);
	
	info() << printMagnetConditions() << endmsg;
	
	//Fill vectors with data in the ntuples
	sc = sc && this->InitializeSmartIDs();
	
	//Open Output File
	TFile *outfile=new TFile(m_outfiletuple.c_str(),"RECREATE");
	
	//Fill peaks graph
	std::vector< TGraph* > hpd_anode_graphs = this->FillTGraphs("peaks_on_anode", 1, m_peaksOnAnode);
	
	
	//Move peaks to the photocathode as with B OFF, false flag means that the old cond DB magnificatio is stored, true for the new ones
	m_peaksOnPhotoCathode = this->MoveToPhotoCathode(m_SmartIDs,m_peaksOnAnode,false);
	//Fill peaks graph on the photocathode
	std::vector< TGraph* > hpd_cathode_graphs = this->FillTGraphs("peaks_on_cathode", 2, m_peaksOnPhotoCathode);
	
	info()<< "Initialising Fit..."<<endmsg;
	
	sc = sc & this->FitGrid(outfile);
// 	sc = sc && this->FitAllHPDS();
	
	if(m_bfield!=0) sc = sc && this->ChangeMagneticField(m_bfield,5850);
	else sc = sc && this->ChangeMagneticField(m_bfield,0.1);
	
// 	sc = sc && this->ChangeMagneticField(m_bfield,5850);
	info() << printMagnetConditions() << endmsg;
	info()<< "Magnifying using old parameters..."<<endmsg;
	
	m_peaksOnPhotoCathode_FIELDON = this->MoveToPhotoCathode(m_SmartIDs,m_peaksOnAnode,true);
	std::vector< TGraph* > hpd_cathode_graphs_FIELDON = this->FillTGraphs("peaks_on_cathode_FIELDON", 4, m_peaksOnPhotoCathode_FIELDON);
	
	info()<< "Updating CondDB..."<<endmsg;
	this->UpdateCondDB();
	
	//so that I can check B=0 stuff (B OFF parameters are in the DDDB... )
	if(m_bfield==0) sc = sc && this->ChangeMagneticField(-1,5850);
	
	info()<< "Magnifying using new parameters..."<<endmsg;
	m_peaksOnPhotoCathode_final = this->MoveToPhotoCathode(m_SmartIDs,m_peaksOnAnode,true);
	std::vector< TGraph* > hpd_cathode_graphs_final = this->FillTGraphs("peaks_on_cathode_final", 9, m_peaksOnPhotoCathode_final);
	
	//peaks in global coordinates
	m_peaksOnPanel_final = this->HPD2Panel(m_SmartIDs,m_peaksOnPhotoCathode_final,true);
	std::vector< TGraph* > hpd_panel_graphs_final = this->FillTGraphs("peaks_on_panel_final", 9, m_peaksOnPanel_final);
// 	sc = sc && this->FillParamHistos();
	
// 	std::vector< std::string > strings;
// 	strings.push_back("Initial");
// 	strings.push_back("FIELDON");
// 	strings.push_back("Final");
	info()<< "Filling resolutions..."<<endmsg;
	this->FillResolutions("Initial",m_peaksOnPhotoCathode);
	this->FillResolutions("OldParams",m_peaksOnPhotoCathode_FIELDON);
	this->FillResolutions("Final",m_peaksOnPhotoCathode_final);
	
	
	info()<< "Making summary plots"<<endmsg;
	this->SummaryPlots(outfile);
	this->GlobalSummary(outfile);
	
	info()<< "Filling HPD Shifts"<<endmsg;
	FillHPDShifts(outfile);
	
	outfile->Close();
	
	return sc;
}

//==============================================================================
//Execute
//==============================================================================
StatusCode MDMFitParameters::execute() {
	StatusCode sc = StatusCode::SUCCESS;
	
	info() << "This algorithm will not process any RAW/DST dataset!" << endmsg;
	
	m_nEvt++;
	
	return sc;
}

//==============================================================================
//Finalize
//==============================================================================
StatusCode MDMFitParameters::finalize() {
	
	const StatusCode sc = Rich::AlgBase::finalize();
	if( sc.isFailure() ) return sc;
	
	info()<< "Parameters fitting is COMPLETED"<<endmsg;
	return sc;
}

//==============================================================================
//Get Magnetic Field Service
//==============================================================================
StatusCode MDMFitParameters::getMagFieldSvc(){
	StatusCode sc = StatusCode::SUCCESS;
	
	
	ISvcLocator* svcLocator = Gaudi::svcLocator();
	if(0 == svcLocator){
		throw GaudiException("ISvcLocator points to NULL!","*DeRichException*",StatusCode::FAILURE);
	}
	
	sc=svcLocator->service("UpdateManagerSvc", m_updMgrSvc);
	
	StatusCode scMag=svcLocator->service("MagneticFieldSvc", m_magFieldSvc);
	if(!scMag.isSuccess()){
		throw GaudiException("Could not locate MagneticFieldSvc!","*DeRichException*",StatusCode::FAILURE);
		sc=scMag;
	}
	
	debug() << "Before playing with OnlineDB" << this->printMagnetConditions() << endreq;
	
	m_updMgrSvc->registerCondition(this, "/dd/Conditions/Online/LHCb/Magnet/Measured",&MDMFitParameters::upd_mag, m_mag_cond);
	
	sc=m_updMgrSvc->update(this);
	
	info() << this->printMagnetConditions() << endmsg;
	
	return sc;
}

//==============================================================================
//return string with magnet conditions
//==============================================================================
std::string MDMFitParameters::printMagnetConditions(){
	std::string output=Form("Magnet Current set to: %f   -   ", m_magFieldSvc->signedRelativeCurrent());
	output+=Form( "Is Magnet Down: %i \n", m_magFieldSvc->isDown());
	return output;
}

//==============================================================================
// Store B OFF MAgnification from DDDB
//==============================================================================
void MDMFitParameters::StoreBOFFMagnification(){
// 	SmartDataPtr<TabulatedProperty> HPDdeMag(dataSvc());
}

//==============================================================================
//Initialize HPDs with MDMS peaks
//==============================================================================
StatusCode MDMFitParameters::InitializeSmartIDs(){
	
	StatusCode sc = StatusCode::SUCCESS;
	
	//Set array size
	m_calibSteps.resize(m_ncols*m_nhpds);
	m_SmartIDs.resize(m_ncols*m_nhpds);
	m_peaksOnAnode.resize(m_ncols*m_nhpds);
	m_peaksOnPhotoCathode.resize(m_ncols*m_nhpds);
	m_hpdCondDB_old.resize(m_ncols*m_nhpds);
	m_hpdCondDB_new.resize(m_ncols*m_nhpds);
	//New parameter error array
	m_hpdCondDB_new_err.resize(m_ncols*m_nhpds);

	m_RadialFunctions.resize(m_ncols*m_nhpds);
	m_AxialFunctions.resize(m_ncols*m_nhpds);
	m_HPDWindowCentresOnPanelRF.resize(m_ncols*m_nhpds);
	m_HPDWindowPivotOnPanelRF.resize(m_ncols*m_nhpds);
	m_HPDWindowPivotOnPanel.resize(m_ncols*m_nhpds);
	m_peaksOnPanel_testing.resize(m_ncols*m_nhpds);
	
	Rich::DetectorType dt=(Rich::DetectorType)0;
	Rich::Side side=(Rich::Side)m_panel;
	
	
	//Load Ntuple
	TFile *infile=new TFile(m_infilepath.c_str());
	info()<< "Opening file \""<< m_infilepath <<"\""<<endmsg;
	
	if(infile){
		info()<< "Searching for tree \""<< m_treepath <<"\""<<endmsg;
		TTree *tree=(TTree*)infile->Get(m_treepath.c_str());
		info()<< "NTuple contains "<<tree->GetEntries()<< " peaks"<<endmsg;
		unsigned int i=0;
		float muX, muY, hpd, col, step;
		tree->SetBranchAddress("oldX",&muX); //uncorrected pixel coordinate
		tree->SetBranchAddress("oldY",&muY); //uncorrected pixel coordinate
		tree->SetBranchAddress("hpd",&hpd);
		tree->SetBranchAddress("col",&col);
		tree->SetBranchAddress("pattern",&step);
		for(i=0; i< tree->GetEntries(); i++){
			tree->GetEntry(i);
			
			int code=this->m_CodeHPD((int)col ,(int)hpd);
// 			debug() << "HPD Code is "<< code << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
			
// 			debug()<< "muX, muY, hpd, column: " << muX << "\t"<<muY<<"\t"<<hpd<<"\t"<<col<<endreq;
// 			debug() << "Array size is "<<m_SmartIDs.size() << "  Current array index is "<< this->m_CodeHPD((int)col,(int)hpd) << endreq;
			
			m_calibSteps[this->m_CodeHPD((int)col,(int)hpd)].push_back((int)step);
			m_SmartIDs[this->m_CodeHPD((int)col,(int)hpd)].push_back(LHCb::RichSmartID(dt, side, (int)hpd, (int)col));
			m_peaksOnAnode[this->m_CodeHPD((int)col,(int)hpd)].push_back(Gaudi::XYZPoint(muX,muY,0.));
		}
	}
	else{
		sc = StatusCode::FAILURE;
	}
	infile->Close();
	
	return sc;
}

//==============================================================================
//Code HPD in Array
//==============================================================================
int MDMFitParameters::m_CodeHPD(int col, int hpd){
	unsigned int array_pos=0;
	array_pos=col*m_nhpds+hpd;
	return array_pos;
}

//==============================================================================
//Code HPD in Array
//==============================================================================
void MDMFitParameters::m_DecodeHPD(int code, int* col, int* hpd){
	*col=code/m_nhpds;
	*hpd=code%m_nhpds;
}

//==============================================================================
//Move To Photocathode
//==============================================================================
std::vector< std::vector< Gaudi::XYZPoint > > MDMFitParameters::MoveToPhotoCathode(std::vector< LHCb::RichSmartID::Vector > SmartIDs, std::vector< std::vector< Gaudi::XYZPoint > > peaksOnAnode, bool iscorrected){
	
	info() << "Moving to Photocathode" << endmsg;


	//Sum of x and y for each HPD (used to work out average x and y)
	double sum_x = 0;
	double sum_y = 0;
	//Average x and y for a HPD
	double mean_x = 0;
	double mean_y = 0;

	//Unshifted x and y
	double x = 0;
	double y = 0;

	//Shifted x and y
	double x_new = 0;
	double y_new = 0;

	//hpd number and column number
	int hpd = 0;
	int col = 0;

	//Peak number index
	int peak_num = 0;

	
	std::vector< std::vector< Gaudi::XYZPoint > > output_vector(m_ncols*m_nhpds);


	//Give each point for a given HPD the same alpha and shift
	//Each HPD gets a different alpha and shift
	//Set random seed for systematcis use

	//Random seed comes from the command line, to ensure that fit i always gets the same seed (so U and D boxes seeded the same)
	int seed = m_randomseed;
	TRandom3 r3(seed);
	//Random radius scaling factor for magnification (can make radius bigger or smaller, but on average keep it the same)
	//Radius change amount specified in percent on the command line. Divide it by 100 here
	float radius_factor = m_radiusfactor/100;
	//double alpha = r3.Gaus(1.0,radius_factor);
	//No radius change
	double alpha = 1.0;
	//Random x and y shift values (on average no shift, but 2mm within 1 sigma) 
	double rand_x = r3.Gaus(0.0,0.01);
	double rand_y = r3.Gaus(0.0,0.01);

	
// 	panel->initialize();
	unsigned int i=0, k=0;
	for(i=0;i<m_ncols*m_nhpds;i++){
	  hpd=-1;
	  col=-1;
 
		this->m_DecodeHPD(i, &col ,&hpd);
		debug() << "(Decoded) HPD Code is "<< i << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		int code=this->m_CodeHPD(col ,hpd);
		debug() << "(Code) HPD Code is "<< code << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		
		const std::string hpdlocation=Form("/dd/Structure/LHCb/BeforeMagnetRegion/Rich1/HPDPanel%i/HPD:%i",m_panel,i+98*m_panel);

		//Quick loop to calculate the average x and y for the HPD in question
		//Loop over all the points, add together the x and y's
		//Divide by the number of points at the end to get the mean position
		for(int j=0; j < SmartIDs[i].size(); j++){

		  sum_x += peaksOnAnode[i][j].x();
		  sum_y += peaksOnAnode[i][j].y();

		}
		mean_x = sum_x/(SmartIDs[i].size());
		mean_y = sum_y/(SmartIDs[i].size());

		//Reset the sum values to zero before the next HPD is done, otherwise the centres are not calculated peroperly
		sum_x = 0;
		sum_y = 0;

// 		IDataProviderSvc* m_detSvc=detSvc();
		
// 		m_updMgrSvc->update(m_detSvc);
		
		DeRichHPD *tmp_hpd= getDet<DeRichHPD>(detSvc(),hpdlocation);
// 		tmp_hpd->initialize();
// 		tmp_hpd->updateMagParInMDMSAlgos();
// 		detSvc()->updateObject(tmp_hpd);
		for(k=0; k < SmartIDs[i].size(); k++){
			Gaudi::XYZPoint tmppoint;

			peak_num = k;


			//Place where the peak positions on the anode are defined
			//Can move them around by some random amount at this point
			//Add a number on average 0, but with sigma 2 mm

			//x with shifting and magnification changes
			//When alpha = 1, no magnification, x = regular x
			//When x = mean_x, nothing added to it as there is no magnification at the centre
			//rand_x and rand_y are random shifts
			x = peaksOnAnode[i][k].x();
			y = peaksOnAnode[i][k].y();
			info() << "Mean x for peaks on HPD "<<i<<" = "<<mean_x<<endmsg;
			info() << "Mean y for peaks on HPD "<<i<<" = "<<mean_y<<endmsg;
			info() << "Normal x = "<<x<<endmsg;
			info() << "Normal y = "<<y<<endmsg;
			info() << "Radius change alpha = "<<alpha<<endmsg;
			info() << "Random x shift = "<<rand_x<<endmsg;
		        info() << "Random y shift = "<<rand_y<<endmsg;
			x_new = alpha*(x - mean_x) + mean_x ;//+ rand_x; 
			y_new = alpha*(y - mean_y) + mean_y ;//+ rand_y; 
			info() << "Shifted x = "<<x_new<<endmsg;
			info() << "Shifted y = "<<y_new<<endmsg;
			//Make x and y the new shifted positions
			x = x_new;
			y = y_new;
			
			tmp_hpd->detectionPoint(x,y,tmppoint,false); //set to false to include refraction, set to true if want to be wrong
			output_vector[i].push_back(tmppoint);

		       
// 			if(i==0 && k==0){
// 				std::cout<<"Test"<<std::endl;
// 				std::cout<<x<<"\t"<<y<<std::endl;
// 				std::cout<<tmppoint.X()<<"\t"<<tmppoint.Y()<<std::endl;
// 			}
		}
       
		
		//Stuff for debugging
		if(SmartIDs[i].size()>0 || 1){
// 			const LHCb::RichSmartID thehpd(m_SmartIDs[i][0]);
// 			const DeRichHPD* thehpdid=panel->deHPD( thehpd );
// 			thehpdid->initialize();
		  const std::string parlist=Form("hpd%i_rec",i+98*m_panel);
		  SmartRef<Condition> hpd_demag_conditions = tmp_hpd->condition("DemagParametersFieldNegative");
		  std::vector<double> params = hpd_demag_conditions->paramAsDoubleVect(parlist);
		  if(i==0) info() << "Mag0 from DeRichHPD\t HPD" <<hpd<<"_COL"<<col << " ->  Demag parameters \t" << params[0] << "\t"<< params[1] << "\t"<< params[2] << "\t"<< params[3] << "\t"<< params[4] << "\t"<< params[5] << "\t"<< params[6] << "\t"<< params[7] << endmsg;
				
			SmartRef<Condition> hpd_demag_conditions_down = tmp_hpd->condition("DemagParametersFieldNegative");
			std::vector<double> params_down = hpd_demag_conditions_down->paramAsDoubleVect(parlist);
			if(i==0) info() << "Mag-1 from DeRichHPD\t HPD" <<hpd<<"_COL"<<col << " ->  Demag parameters \t" << params_down[0] << "\t"<< params_down[1] << "\t"<< params_down[2] << "\t"<< params_down[3] << "\t"<< params_down[4] << "\t"<< params_down[5] << "\t"<< params_down[6] << "\t"<< params_down[7] << endmsg;
				
			SmartRef<Condition> hpd_demag_conditions_up = tmp_hpd->condition("DemagParametersFieldPositive");
			std::vector<double> params_up = hpd_demag_conditions_up->paramAsDoubleVect(parlist);
			if(i==0) info() << "Mag1 from DeRichHPD\t HPD" <<hpd<<"_COL"<<col << " ->  Demag parameters \t" << params_up[0] << "\t"<< params_up[1] << "\t"<< params_up[2] << "\t"<< params_up[3] << "\t"<< params_up[4] << "\t"<< params_up[5] << "\t"<< params_up[6] << "\t"<< params_up[7] << endmsg;
			
			debug() << "New Parameters" << endreq;
			if(m_hpdCondDB_new[i].size()==8){
				debug() << "Mag from fit\t HPD" <<hpd<<"_COL"<<col << " ->  Demag parameters \t" << m_hpdCondDB_new[i][0] << "\t"<< m_hpdCondDB_new[i][1] << "\t"<< m_hpdCondDB_new[i][2] << "\t"<< m_hpdCondDB_new[i][3] << "\t"<< m_hpdCondDB_new[i][4] << "\t"<< m_hpdCondDB_new[i][5] << "\t"<< m_hpdCondDB_new[i][6] << "\t"<< m_hpdCondDB_new[i][7] << endreq;
			}
			
			if(!iscorrected){
				
				if(m_bfield==-1) m_hpdCondDB_old[i]=params_down;
				if(m_bfield==1) m_hpdCondDB_old[i]=params_up;
				if(m_bfield==0) m_hpdCondDB_old[i]=params;
				
// 				std::cout << "old "<< i << std::endl; 
			}
// 			delete thehpdid;
		}
		
// 		delete tmp_hpd;
	}
// 	delete panel;
	return output_vector;
}

//==============================================================================
//Move To Photocathod
//==============================================================================
std::vector< std::vector< Gaudi::XYZPoint > > MDMFitParameters::HPD2Panel(std::vector< LHCb::RichSmartID::Vector > SmartIDs, std::vector< std::vector< Gaudi::XYZPoint > > peaksOnWindow, bool iscorrected){
	
	info() << "Moving to Photocathode" << endmsg;
	
	std::vector< std::vector< Gaudi::XYZPoint > > output_vector(m_ncols*m_nhpds);
	
// 	panel->initialize();
	unsigned int i=0, k=0;
	for(i=0;i<m_ncols*m_nhpds;i++){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
		debug() << "(Decoded) HPD Code is "<< i << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		int code=this->m_CodeHPD(col ,hpd);
		debug() << "(Code) HPD Code is "<< code << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		
		const std::string hpdlocation=Form("/dd/Structure/LHCb/BeforeMagnetRegion/Rich1/HPDPanel%i/HPD:%i",m_panel,i+98*m_panel);
		
		DeRichHPD *tmp_hpd= getDet<DeRichHPD>(detSvc(),hpdlocation);
// 		tmp_hpd->updateMagParInMDMSAlgos();
		
		if(fitted[i])m_HPDWindowPivotOnPanel[i]=(tmp_hpd->fromHPDToPanel())*MDMGrids[i]->GetPivotCoord();
// 		detSvc()->updateObject(tmp_hpd);
		for(k=0; k < SmartIDs[i].size(); k++){
			Gaudi::XYZPoint tmppoint;
			tmppoint=(tmp_hpd->fromHPDToPanel())*peaksOnWindow[i][k];
// 			kop
 			//double x= peaksOnWindow[i][k].x();
 			//double y= peaksOnWindow[i][k].y();
 			//tmp_hpd->detectionPoint(x,y,tmppoint,false); //set to false to include refraction, set to true if want to be wrong
			
			output_vector[i].push_back(tmppoint);
		}
		
		Gaudi::XYZPoint testpoint;
		testpoint.SetXYZ(4,0,0);
		m_peaksOnPanel_testing[i]=(tmp_hpd->fromHPDToPanel())*testpoint;
	}
// 	delete panel;
	return output_vector;
}

//==============================================================================
//Fill a TGraph with the peaks
//==============================================================================
std::vector< TGraph* > MDMFitParameters::FillTGraphs(std::string Prefix, int color, std::vector< std::vector< Gaudi::XYZPoint > > peaks){
	std::vector< TGraph* > output_graph;
	unsigned int i=0,k=0;
	
	for(i=0; i<peaks.size(); i++){
		
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
		
		unsigned int hpd_npeaks=peaks[i].size();
		double *xcoord = new double[hpd_npeaks];
		double *ycoord = new double[hpd_npeaks];
		for(k=0; k<hpd_npeaks; k++){
			xcoord[k]=peaks[i][k].X();
			ycoord[k]=peaks[i][k].Y();
		}
		
		TGraph *g=new TGraph(hpd_npeaks,xcoord,ycoord);
		g->SetName(Form("%s_HPD%i_COL%i",Prefix.c_str(),hpd,col));
		g->SetTitle(Form("%s_HPD%i_COL%i - Peak Map, Outside Quartz Window",Prefix.c_str(),hpd,col));
		g->GetXaxis()->SetTitle("X (mm)");
		g->GetYaxis()->SetTitle("Y (mm)");
		g->SetMarkerColor(color);
		g->SetMarkerStyle(2);
		g->SetMarkerSize(1);
		output_graph.push_back(g);
		g->Write();
	}
	return output_graph;
}


//==============================================================================
//Update Magnetic Field
//==============================================================================
StatusCode MDMFitParameters::ChangeMagneticField(int polarity, double current){
	StatusCode sc = StatusCode::SUCCESS;
	
	info() << "Acting on Magnet Settings  -  Current is: " << m_mag_cond->paramAsDouble("Current") << ",  Polarity is: " << m_mag_cond->paramAsInt("Polarity") << endmsg;
	
	
	Condition newMagnetConditions;
	newMagnetConditions.addParam<double>("Current",current);
	newMagnetConditions.addParam<int>("Polarity",polarity);
	
	m_mag_cond->update(newMagnetConditions);
	m_updMgrSvc->invalidate(m_mag_cond);
	
	sc = sc && m_updMgrSvc->update(m_mag_cond);
	
// 	m_updMgrSvc->registerCondition(this, m_magFieldSvc,&MDMFitParameters::upd_mag);
	m_updMgrSvc->update(m_magFieldSvc);
	info() << "Magnet Settings Changed -  Current is: " << m_mag_cond->paramAsDouble("Current") << ",  Polarity is: " << m_mag_cond->paramAsInt("Polarity") << endmsg;
	
	return sc;
}

//==============================================================================
//Fill the vector of new parameters given the fitted TF1
//==============================================================================
StatusCode MDMFitParameters::fromTF1toDoubleVector(int hpdnumber, TF1 *radialF, TF1* axialF, bool includemag){
	StatusCode sc = StatusCode::SUCCESS;
	std::vector<double> tmpcond(8); //8 MDMS parameters
	double r0=radialF->GetParameter(0);
	double r1=radialF->GetParameter(1);
	double r2=radialF->GetParameter(2);
	double r3=radialF->GetParameter(3);
	double a0=axialF->GetParameter(0);
	double a1=axialF->GetParameter(1);
	double a2=axialF->GetParameter(2);
	double a3=axialF->GetParameter(3);

	std::vector<double> tmpcond_err(8); //8 MDMS parameter errors
	//Parameter errors for MDMS parameters
	double r0_err=radialF->GetParError(0);
	double r1_err=radialF->GetParError(1);
	double r2_err=radialF->GetParError(2);
	double r3_err=radialF->GetParError(3);
	double a0_err=axialF->GetParError(0);
	double a1_err=axialF->GetParError(1);
	double a2_err=axialF->GetParError(2);
	double a3_err=axialF->GetParError(3);
	
	
	info() << "HPD #"<<hpdnumber+98*m_panel<<endmsg;
	info() << "Radial correction fitted is "<<r0<<"\t"<<r1<<"\t"<<r2<<"\t"<<r3<<endmsg;
	info() << "Axial correction fitted is "<<a0<<"\t"<<a1<<"\t"<<a2<<"\t"<<a3<<endmsg;
	
	if(!includemag){
		tmpcond[0]=r0;
		tmpcond[1]=r1;
		tmpcond[2]=r2;
		tmpcond[3]=r3;
		tmpcond[4]=a0;
		tmpcond[5]=a1;
		tmpcond[6]=a2;
		tmpcond[7]=a3;

		tmpcond_err[0]=r0_err;
		tmpcond_err[1]=r1_err;
		tmpcond_err[2]=r2_err;
		tmpcond_err[3]=r3_err;
		tmpcond_err[4]=a0_err;
		tmpcond_err[5]=a1_err;
		tmpcond_err[6]=a2_err;
		tmpcond_err[7]=a3_err;
	}
	else{
		tmpcond[0]=r0;
		tmpcond[1]= (r1+1)*m_defmag;
		tmpcond[2]= r2*m_defmag*m_defmag + (r1+1)*m_defmag_o2;
		tmpcond[3]= r3*m_defmag*m_defmag*m_defmag + 2*r2*m_defmag_o2*m_defmag + (r1+1)*m_defmag_o3;
		tmpcond[4]=a0;
		tmpcond[5]=a1*tmpcond[1];
		tmpcond[6]=a2*tmpcond[1]*tmpcond[1]+a1*tmpcond[2];
		tmpcond[7]=a3*tmpcond[1]*tmpcond[1]*tmpcond[1]+2*a2*tmpcond[1]*tmpcond[2]+a1*tmpcond[3];

		tmpcond_err[0]=r0_err;
		tmpcond_err[1]=r1_err;
		tmpcond_err[2]=r2_err;
		tmpcond_err[3]=r3_err;
		tmpcond_err[4]=a0_err;
		tmpcond_err[5]=a1_err;
		tmpcond_err[6]=a2_err;
		tmpcond_err[7]=a3_err;


		info() << "HPD #"<<hpdnumber+98*m_panel<<endmsg;
		info() << "Radial correction is "<<tmpcond[0]<<"\t"<<tmpcond[1]<<"\t"<<tmpcond[2]<<"\t"<<tmpcond[3]<<endmsg;
		info() << "Axial correction is "<<tmpcond[4]<<"\t"<<tmpcond[5]<<"\t"<<tmpcond[6]<<"\t"<<tmpcond[7]<<endmsg;
	}
	m_hpdCondDB_new[hpdnumber]=tmpcond;
	m_hpdCondDB_new_err[hpdnumber]=tmpcond_err;
	return sc;
}

//==============================================================================
//Update Conditions Database
//==============================================================================
StatusCode MDMFitParameters::UpdateCondDB(){
	StatusCode sc = StatusCode::SUCCESS;
	
	info() << "Updating parameters in the CondDB..." << endmsg;
	
	unsigned int i=0;
	
	
	std::string panel_location;
	if(m_panel==0) panel_location="/dd/Structure/LHCb/BeforeMagnetRegion/Rich1/HPDPanel0";
	else panel_location="/dd/Structure/LHCb/BeforeMagnetRegion/Rich1/HPDPanel1";
	
	std::string fieldstring="DemagParameters";
	if(m_bfield==-1) fieldstring="DemagParametersFieldNegative";
	if(m_bfield==1) fieldstring="DemagParametersFieldPositive";
	if(m_bfield==0) fieldstring="DemagParametersFieldNegative";
// 	fieldstring="DemagParametersFieldNegative";
	
	DeRichHPDPanel *panel=getDet<DeRichHPDPanel>(panel_location);
	panel->initialize();
	
// 	info() << parInCondDB->printParams() << endmsg;
// 	const LHCb::RichSmartID hpd_id(m_SmartIDs[0][0]); //just need a DeRichHPD to get the condition
// 	const DeRichHPD* hpd=panel->deHPD( 10 ); //just need a DeRichHPD to get the condition
// 	info() << panel->nPDColumns() << "\t"<< panel->nPDsPerCol() << endmsg;
// 	info() << panel->nPDs() << endmsg;
// 	const DeRichHPD* hpd=NULL;
	const std::string hpdlocation=Form("/dd/Structure/LHCb/BeforeMagnetRegion/Rich1/HPDPanel%i/HPD:%i",m_panel,10+98*m_panel);
// 		DeRichHPD* tmp_hpd=panel->deHPD( i );
	DeRichHPD *tmp_hpd= getDet<DeRichHPD>(detSvc(),hpdlocation);
	Condition* parInCondDB = tmp_hpd->condition(fieldstring.c_str());
	debug() << "Previous conditions"<< endreq;
	debug() << parInCondDB->printParams() << endreq;
	Condition newPars;
	
	for(i=0;i<m_nhpds*m_ncols;i++){
		
		int nhpd=-1,ncol=-1;
		debug() << "Update demag conditions for HPD"<< i << endreq;
		this->m_DecodeHPD(i, &ncol ,&nhpd);
		
		std::vector<double> params = parInCondDB->paramAsDoubleVect(Form("hpd%i_rec",i+98*m_panel));
		
		debug() << "Mag0\t HPD" <<nhpd<<"_COL"<<ncol << " ->  Demag parameters \t" << params[0] << "\t"<< params[1] << "\t"<< params[2] << "\t"<< params[3] << "\t"<< params[4] << "\t"<< params[5] << "\t"<< params[6] << "\t"<< params[7] << endreq;
		
		
// 		newPars.addParam< std::vector<double> >(Form("hpd%i_rec",i+98*(m_panel-1)),m_hpdCondDB_new[i]);
		newPars.addParam< std::vector<double> >(Form("hpd%i_rec",i+98*m_panel),m_hpdCondDB_new[i]);
		
		debug() << "New conditions"<< endreq;
		params = parInCondDB->paramAsDoubleVect(Form("hpd%i_rec",i+98*m_panel));
		
		debug() << "Mag0\t HPD" <<nhpd<<"_COL"<<ncol << " ->  Demag parameters \t" << params[0] << "\t"<< params[1] << "\t"<< params[2] << "\t"<< params[3] << "\t"<< params[4] << "\t"<< params[5] << "\t"<< params[6] << "\t"<< params[7] << endreq;
		
		//adding dummy sim hpds so that the DB does not get upset
		std::vector<double> params_sim = parInCondDB->paramAsDoubleVect(Form("hpd%i_sim",i+98*m_panel));
		newPars.addParam< std::vector<double> >(Form("hpd%i_sim",i+98*m_panel),params_sim);
	}
	
	std::string condition_name="";
	std::string condname_up=Form("demagParsR1P%i_FieldPositive",m_panel);
	std::string condname_down=Form("demagParsR1P%i_FieldNegative",m_panel);
	std::string condname_off=Form("demagParsR1P%i_FieldNegative",m_panel);
	if(m_bfield==-1) condition_name=condname_down;
	if(m_bfield==1) condition_name=condname_up;
	if(m_bfield==0) condition_name=condname_off;
	
	info() << "Before"<<endmsg;
	std::string file_name=condition_name;
	if(m_bfield==-1) file_name=Form("demagParsR1P%i_FieldNegative_Run%i",m_panel,m_run);
	if(m_bfield==1) file_name=Form("demagParsR1P%i_FieldPositive_Run%i",m_panel,m_run);
	if(m_bfield==0) file_name=Form("demagParsR1P%i_FieldOff_Run%i",m_panel,m_run);
// 	info() << parInCondDB->toXml() << endmsg;
	std::ofstream outbefore((m_outfilepath+".previous").c_str());
	outbefore << parInCondDB->toXml(condition_name,true,6);
	outbefore.close();
	
// 	m_updMgrSvc->purge();
	parInCondDB->update(newPars);
	m_updMgrSvc->invalidate(parInCondDB);
	sc = sc && m_updMgrSvc->update(parInCondDB);
	
// 	m_updMgrSvc->registerCondition(condition_name);
// 	runUpdate();
// 	sc = sc && m_updMgrSvc->update(detSvc());
	debug() << parInCondDB->printParams() << endreq;
	info() << "After"<<endmsg;
// 	info() << parInCondDB->toXml() << endmsg;
	std::ofstream outafter((m_outfilepath).c_str());
	outafter << parInCondDB->toXml(condition_name,true,6);
	outafter.close();
	info() << m_outfilepath<<endmsg;
	info() << "xml outfiles written"<<endmsg;
	for(i=0;i<m_ncols*m_nhpds;i++){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
		debug() << "(Decoded) HPD Code is "<< i << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		int code=this->m_CodeHPD(col ,hpd);
		debug() << "(Code) HPD Code is "<< code << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		
		const std::string hpdlocation=Form("/dd/Structure/LHCb/BeforeMagnetRegion/Rich1/HPDPanel%i/HPD:%i",m_panel,i+98*m_panel);
// 		DeRichHPD* tmp_hpd=panel->deHPD( i );
		DeRichHPD *tmp_hpd= getDet<DeRichHPD>(detSvc(),hpdlocation);
		info() << hpdlocation << endmsg;
		m_updMgrSvc->update(tmp_hpd);
	}
	//should not be hard
// 	Condition* parInCondDB = 
	
	info() << "DONE" << endmsg;
	
	return sc;
}


//==============================================================================
// Fit Grid
//==============================================================================
StatusCode MDMFitParameters::FitGrid(TFile *infile){
	StatusCode sc = StatusCode::SUCCESS;
	unsigned int i=0,k=0;
	info() << "Fitting..." << endmsg;
	for(i=0;i<m_nhpds*m_ncols;i++){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
		
		info() << Form("Initialising Grid_HPD%i_COL%i",hpd,col) << endmsg;
		
		MDMGridModel *grid=new MDMGridModel(Form("Grid_HPD%i_COL%i",hpd,col),m_peaksOnPhotoCathode[i], m_calibSteps[i], 14, 0.1974, m_gstep,hpd,col);
		
		bool fitdone=grid->FitParameters();
		fitted.push_back(fitdone);
		
		MDMGrids.push_back(grid);
		
		std::vector<double> pars=grid->GetParameters();
		std::vector<double> pars_err=grid->GetParameters_err();
		
		info() << "Creating functions" << endmsg;
		
		TF1* radialF = new TF1(Form("Radial_TF1_HPD%i_COL%i",hpd,col),"pol3(0)",0.0,30.0);
		for(k=0;k<4;k++){
			radialF->SetParameter(k,pars[k]);
			radialF->SetParError(k,pars_err[k]);
		}
		TF1* axialF = new TF1(Form("Radial_TF1_HPD%i_COL%i",hpd,col),"pol3(0)",0.0,30.0);
		for(k=0;k<4;k++){
			axialF->SetParameter(k,pars[k+4]);
			//info() << "Axial parameter "<<k<<" error is "<<pars_err[k+4]<< endmsg;
			axialF->SetParError(k,pars_err[k+4]);
		}
		
		m_RadialFunctions.push_back(radialF);
		m_AxialFunctions.push_back(axialF);
		fromTF1toDoubleVector(i, radialF, axialF, true);
// 		tmpcanvas2->cd(2);
// 		h[0]->SetName(Form("resX_HPD%i_COL%i",hpd,col));
// 		h[1]->SetName(Form("resY_HPD%i_COL%i",hpd,col));
// 		h[0]->SetLineColor(2);
// 		h[0]->Draw();
// 		h[1]->Draw("SAMES");
// 		tmpcanvas2->Write();
	}
	
	info() << "Fit COMPLETED" << endmsg;
	return sc;
}
//==============================================================================
// Fill Resolutions
//==============================================================================
StatusCode MDMFitParameters::FillResolutions(TString Prefix, std::vector< std::vector< Gaudi::XYZPoint > > peaks){
	StatusCode sc = StatusCode::SUCCESS;
	
	unsigned int i=0,k=0;
	for(i=0;i<m_nhpds*m_ncols;i++){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
		TH1D *htmp = new TH1D(Form("%s_Resolution_HPD%i_COL%i",Prefix.Data(),hpd,col),Form("%s_Resolution_HPD%i_COL%i; Residual (mm);",Prefix.Data(),hpd,col),100,0,15);
		TH1D *htmpX = new TH1D(Form("%s_ResolutionX_HPD%i_COL%i",Prefix.Data(),hpd,col),Form("%s_ResolutionX_HPD%i_COL%i; Residual X (mm);",Prefix.Data(),hpd,col),100,-10,10);
		TH1D *htmpY = new TH1D(Form("%s_ResolutionY_HPD%i_COL%i",Prefix.Data(),hpd,col),Form("%s_ResolutionY_HPD%i_COL%i; Residual Y (mm);",Prefix.Data(),hpd,col),100,-10,10);
		
		MDMGrids[i]->GetResolution(peaks[i],htmp,"D");
		MDMGrids[i]->GetResolution(peaks[i],htmpX,"X");
		MDMGrids[i]->GetResolution(peaks[i],htmpY,"Y");
		
		htmp->Write();
		htmpX->Write();
		htmpY->Write();
		
		if(Prefix=="Final"){
			TH1D *htmp_self = new TH1D(Form("%s_Resolution_HPD%i_COL%i","Self",hpd,col),Form("%s_Resolution_HPD%i_COL%i; Residual (mm);",Prefix.Data(),hpd,col),100,0,15);
			TH1D *htmpX_self = new TH1D(Form("%s_ResolutionX_HPD%i_COL%i","Self",hpd,col),Form("%s_ResolutionX_HPD%i_COL%i; Residual X (mm);",Prefix.Data(),hpd,col),100,-10,10);
			TH1D *htmpY_self = new TH1D(Form("%s_ResolutionY_HPD%i_COL%i","Self",hpd,col),Form("%s_ResolutionY_HPD%i_COL%i; Residual Y (mm);",Prefix.Data(),hpd,col),100,-10,10);
		
			MDMGrids[i]->GetSelfResolution(htmp_self,"D");
			MDMGrids[i]->GetSelfResolution(htmpX_self,"X");
			MDMGrids[i]->GetSelfResolution(htmpY_self,"Y");
			
			htmp_self->Write();
			htmpX_self->Write();
			htmpY_self->Write();
		
		}
	}
	
	
	return sc;
}

//==============================================================================
// Summary Plots
//==============================================================================
StatusCode MDMFitParameters::SummaryPlots(TFile *infile){
	StatusCode sc = StatusCode::SUCCESS;
	
	unsigned int i=0,k=0;
	for(i=0;i<m_nhpds*m_ncols;i++){
		info()<< "Filling Summary for HPD "<< i <<endmsg;
		if(fitted[i]){
			int hpd=-1,col=-1;
			this->m_DecodeHPD(i, &col ,&hpd);
			TH1D *htmp_init =(TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","Initial",hpd,col));
			TH1D *htmpX_init =(TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","Initial",hpd,col));
			TH1D *htmpY_init =(TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","Initial",hpd,col));
			htmp_init->SetLineColor(1);
			htmpX_init->SetLineColor(1);
			htmpY_init->SetLineColor(1);
			TH1D *htmp_old =(TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","OldParams",hpd,col));
			TH1D *htmpX_old =(TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","OldParams",hpd,col));
			TH1D *htmpY_old =(TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","OldParams",hpd,col));
			htmp_old->SetLineColor(2);
			htmpX_old->SetLineColor(2);
			htmpY_old->SetLineColor(2);
			TH1D *htmp_final =(TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","Final",hpd,col));
			TH1D *htmpX_final =(TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","Final",hpd,col));
			TH1D *htmpY_final =(TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","Final",hpd,col));
			htmp_final->SetLineColor(4);
			htmpX_final->SetLineColor(4);
			htmpY_final->SetLineColor(4);
			htmp_final->SetLineWidth(2);
			htmpX_final->SetLineWidth(2);
			htmpY_final->SetLineWidth(2);
			
			
			TLegend *legend=new TLegend(0.7,0.2,0.99,0.5);
			legend->AddEntry(htmp_init,"No Distortion Correction","l");
			legend->AddEntry(htmp_old,"Current Parameters","l");
			legend->AddEntry(htmp_final,"New Parameters","l");
			
			TGraph *g_prev = (TGraph*)infile->Get(Form("%s_HPD%i_COL%i","peaks_on_cathode",hpd,col));
			TGraph *g_old = (TGraph*)infile->Get(Form("%s_HPD%i_COL%i","peaks_on_cathode_FIELDON",hpd,col));
			TGraph *g_final = (TGraph*)infile->Get(Form("%s_HPD%i_COL%i","peaks_on_cathode_final",hpd,col));
			
			TLatex latex;
			latex.SetTextSize(0.027);
			latex.SetTextAlign(13);//align at top
			TLatex latex1;
			latex1.SetTextSize(0.035);
			latex1.SetTextAlign(13);//align at top
			latex1.SetTextColor(1);
			TLatex latex2;
			latex2.SetTextSize(0.035);
			latex2.SetTextAlign(13);//align at top
			latex2.SetTextColor(14);
			
			TCanvas *sum_canvas=new TCanvas(Form("Summary_HPD%i_COL%i",hpd,col),Form("Summary_HPD%i_COL%i",hpd,col),1000,1000);
			sum_canvas->Divide(3,2);
			sum_canvas->cd(1);
			g_prev->Draw("AP");
			g_prev->GetYaxis()->SetRangeUser(-50.,50.);
			g_prev->GetXaxis()->SetLimits(-50.,50.);
			g_prev->Draw("AP");
			latex.DrawLatex(-45,48,("B OFF magnification"));
			latex.DrawLatex(-45,44,Form("Linear coeff: %0.3g ,  quadratic coeff: %0.3g",m_defmag,m_defmag_o2));
			MDMGrids[i]->DrawGrid("SAME",m_peaksOnPhotoCathode[i],-50.,50.,-50.,50.);
			latex1.DrawLatex(33,-38,("Grid 1"));
			latex2.DrawLatex(33,-43,("Grid 2"));
			sum_canvas->cd(2);
			g_old->Draw("AP");
			g_old->GetYaxis()->SetRangeUser(-50.,50.);
			g_old->GetXaxis()->SetLimits(-50.,50.);
			g_old->Draw("AP");
			latex.DrawLatex(-45,48,("Old parameters"));
			latex.DrawLatex(-45,44,Form("Radial params: %0.3g+%0.3g*R+%0.3g*R^{2}+%0.3g*R^{3}",m_hpdCondDB_old[i][0],m_hpdCondDB_old[i][1],m_hpdCondDB_old[i][2],m_hpdCondDB_old[i][3]));
			latex.DrawLatex(-45,40,Form("Axial params: %0.3g+%0.3g*R+%0.3g*R^{2}+%0.3g*R^{3}",m_hpdCondDB_old[i][4],m_hpdCondDB_old[i][5],m_hpdCondDB_old[i][6],m_hpdCondDB_old[i][7]));
			MDMGrids[i]->DrawGrid("SAME",m_peaksOnPhotoCathode_FIELDON[i],-50.,50.,-50.,50.);
			latex1.DrawLatex(33,-38,("Grid 1"));
			latex2.DrawLatex(33,-43,("Grid 2"));
			sum_canvas->cd(3);
			g_final->Draw("AP");
			g_final->GetYaxis()->SetRangeUser(-50.,50.);
			g_final->GetXaxis()->SetLimits(-50.,50.);
			g_final->Draw("AP");
			latex.DrawLatex(-45,48,("New parameters"));
			latex.DrawLatex(-45,44,Form("Radial params: %0.3g+%0.3g*R+%0.3g*R^{2}+%0.3g*R^{3}",m_hpdCondDB_new[i][0],m_hpdCondDB_new[i][1],m_hpdCondDB_new[i][2],m_hpdCondDB_new[i][3]));
			latex.DrawLatex(-45,40,Form("Axial params: %0.3g+%f*R+%0.3g*R^{2}+%0.3g*R^{3}",m_hpdCondDB_new[i][4],m_hpdCondDB_new[i][5],m_hpdCondDB_new[i][6],m_hpdCondDB_new[i][7]));
			MDMGrids[i]->DrawGrid("SAME",m_peaksOnPhotoCathode_final[i],-50.,50.,-50.,50.);
			latex1.DrawLatex(33,-38,("Grid 1"));
			latex2.DrawLatex(33,-43,("Grid 2"));
			sum_canvas->cd(4);
			htmp_final->SetMaximum(TMath::Max(htmp_final->GetMaximum(),TMath::Max(htmp_old->GetMaximum(),htmp_init->GetMaximum())) * 1.1);
			htmp_final->SetMinimum(0);
			htmp_final->Draw();
			htmp_old->Draw("SAMES");
			htmp_init->Draw("SAMES");
			legend->Draw("SAME");
			sum_canvas->cd(5);
			htmpX_final->SetMaximum(TMath::Max(htmpX_final->GetMaximum(),TMath::Max(htmpX_old->GetMaximum(),htmpX_init->GetMaximum())) * 1.1);
			htmpX_final->SetMinimum(0);
			htmpX_final->Draw();
			htmpX_old->Draw("SAMES");
			htmpX_init->Draw("SAMES");
			legend->Draw("SAME");
			sum_canvas->cd(6);
			htmpY_final->SetMaximum(TMath::Max(htmpY_final->GetMaximum(),TMath::Max(htmpY_old->GetMaximum(),htmpY_init->GetMaximum())) * 1.1);
			htmpY_final->SetMinimum(0);
			htmpY_final->Draw();
			htmpY_old->Draw("SAMES");
			htmpY_init->Draw("SAMES");
			legend->Draw("SAME");
			sum_canvas->Write();
		}
	}
	
	
	return sc;
}

//==============================================================================
// Global Summary
//==============================================================================
StatusCode MDMFitParameters::GlobalSummary(TFile *infile){
	StatusCode sc = StatusCode::SUCCESS;

	info() << "Filling Global Summary"<<endmsg;
	
	TH1D *htmp_init = new TH1D("Initial_Resolution","Initial_Resolution;Residual (mm);",100,0,15);
	TH1D *htmpX_init = new TH1D("Initial_Resolution_X","Initial_Resolution_X;Residual (mm);",100,-10,10);
	TH1D *htmpY_init = new TH1D("Initial_Resolution_Y","Initial_Resolution_Y;Residual (mm);",100,-10,10);
	TH1D *htmp_old = new TH1D("Old_Resolution","Old_Resolution;Residual (mm);",100,0,15);
	TH1D *htmpX_old = new TH1D("Old_Resolution_X","Old_Resolution_X;Residual (mm);",100,-10,10);
	TH1D *htmpY_old = new TH1D("Old_Resolution_Y","Old_Resolution_Y;Residual (mm);",100,-10,10);
	TH1D *htmp_final = new TH1D("final_Resolution","final_Resolution;Residual (mm);",100,0,15);
	TH1D *htmpX_final = new TH1D("final_Resolution_X","final_Resolution_X;Residual (mm);",100,-10,10);
	TH1D *htmpY_final = new TH1D("final_Resolution_Y","final_Resolution_Y;Residual (mm);",100,-10,10);
	
	TH1D *htmp_self = new TH1D("self_Resolution","self_Resolution;Residual (mm);",100,0,15);
	TH1D *htmpX_self = new TH1D("self_Resolution_X","self_Resolution_X;Residual (mm);",100,-10,10);
	TH1D *htmpY_self = new TH1D("self_Resolution_Y","self_Resolution_Y;Residual (mm);",100,-10,10);
	
	htmp_init->SetLineColor(1);
	htmpX_init->SetLineColor(1);
	htmpY_init->SetLineColor(1);
	htmp_old->SetLineColor(2);
	htmpX_old->SetLineColor(2);
	htmpY_old->SetLineColor(2);
	htmp_final->SetLineColor(4);
	htmpX_final->SetLineColor(4);
	htmpY_final->SetLineColor(4);
	htmp_final->SetLineWidth(2);
	htmpX_final->SetLineWidth(2);
	htmpY_final->SetLineWidth(2);
	
	htmp_self->SetLineColor(1);
	htmpX_self->SetLineColor(2);
	htmpY_self->SetLineColor(4);
	
	TLegend *legend=new TLegend(0.6,0.2,0.99,0.5);
	legend->AddEntry(htmp_init,"No Magnetic Correction","l");
	legend->AddEntry(htmp_old,"Current Parameters","l");
	legend->AddEntry(htmp_final,"New Parameters","l");
	
	TLegend *legend2=new TLegend(0.7,0.2,0.99,0.5);
	legend2->AddEntry(htmp_init,"Self resolution","l");
	legend2->AddEntry(htmp_old,"Self Resolution X","l");
	legend2->AddEntry(htmp_final,"Self resolution Y","l");
	
	TH1D* h_r0=new TH1D("r0_dist","r0_dist;r0;",100,-1,1);
	TH1D* h_r1=new TH1D("r1_dist","r1_dist;r1;",100,4.5,6);
	TH1D* h_r2=new TH1D("r2_dist","r2_dist;r2;",100,-0.3,0.3);
	TH1D* h_r3=new TH1D("r3_dist","r3_dist;r3;",100,-0.1,0.1);
	
	TH1D* h_a0=new TH1D("a0_dist","a0_dist;a0;",100,-1,1);
	TH1D* h_a1=new TH1D("a1_dist","a1_dist;a1;",100,-0.4,0.4);
	TH1D* h_a2=new TH1D("a2_dist","a2_dist;a2;",100,-0.02,0.02);
	TH1D* h_a3=new TH1D("a3_dist","a3_dist;a3;",100,-0.002,0.002);
	
	TH2D* h_r0_2D=new TH2D("r0_2d","r0_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_r1_2D=new TH2D("r1_2d","r1_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_r2_2D=new TH2D("r2_2d","r2_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_r3_2D=new TH2D("r3_2d","r3_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_a0_2D=new TH2D("a0_2d","a0_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_a1_2D=new TH2D("a1_2d","a1_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_a2_2D=new TH2D("a2_2d","a2_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_a3_2D=new TH2D("a3_2d","a3_2d;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	
	TH2D* h_Res_2D=new TH2D("Res_2D","Res_2D;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_ResX_2D=new TH2D("ResX_2D","ResX_2D;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	TH2D* h_ResY_2D=new TH2D("ResY_2D","ResY_2D;HPD;COL",m_nhpds,0,m_nhpds,m_ncols,0,m_ncols);
	
	info() << "Opening output Ntuple"<<endmsg;
	
	TTree *outtree=new TTree("parameters","parameters");
	int hpd_t,col_t,isfit;
	double pr0,pr1,pr2,pr3;
	double pa0,pa1,pa2,pa3;
	//Param errors
	double pr0_err,pr1_err,pr2_err,pr3_err;
	double pa0_err,pa1_err,pa2_err,pa3_err;
	//Grid shifts
	double gs0,gs1;
	double gs0_err,gs1_err;
	outtree->Branch("hpd",&hpd_t,"hpd/I");
	outtree->Branch("col",&col_t,"col/I");
	outtree->Branch("fitted",&isfit,"fitted/I");
	outtree->Branch("r0",&pr0,"r0/D");
	outtree->Branch("r1",&pr1,"r1/D");
	outtree->Branch("r2",&pr2,"r2/D");
	outtree->Branch("r3",&pr3,"r3/D");
	outtree->Branch("a0",&pa0,"a0/D");
	outtree->Branch("a1",&pa1,"a1/D");
	outtree->Branch("a2",&pa2,"a2/D");
	outtree->Branch("a3",&pa3,"a3/D");

	//Fill errors for the MDMS params, so I can do pull distributions
	outtree->Branch("r0_err",&pr0_err,"r0_err/D");
	outtree->Branch("r1_err",&pr1_err,"r1_err/D");
	outtree->Branch("r2_err",&pr2_err,"r2_err/D");
	outtree->Branch("r3_err",&pr3_err,"r3_err/D");
	outtree->Branch("a0_err",&pa0_err,"a0_err/D");
	outtree->Branch("a1_err",&pa1_err,"a1_err/D");
	outtree->Branch("a2_err",&pa2_err,"a2_err/D");
	outtree->Branch("a3_err",&pa3_err,"a3_err/D");

	//Add the grid x and y shifts
	outtree->Branch("grid_x_shift",&gs0,"gs0/D");
	outtree->Branch("grid_y_shift",&gs1,"gs1/D");
	outtree->Branch("grid_x_shift_err",&gs0_err,"gs0_err/D");
	outtree->Branch("grid_y_shift_err",&gs1_err,"gs1_err/D");
	
	unsigned int i=0,k=0;
	for(i=0;i<m_nhpds*m_ncols;i++){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
		
		hpd_t=hpd;
		col_t=col;
		
// 		std::cout << "Filling: " << fitted[i] << std::endl;
		
		if(fitted[i]){
			htmp_init->Add((TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","Initial",hpd,col)));
			htmpX_init->Add((TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","Initial",hpd,col)));
			htmpY_init->Add((TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","Initial",hpd,col)));
			
			htmp_old->Add((TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","OldParams",hpd,col)));
			htmpX_old->Add((TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","OldParams",hpd,col)));
			htmpY_old->Add((TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","OldParams",hpd,col)));
			
			htmp_final->Add((TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","Final",hpd,col)));
			htmpX_final->Add((TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","Final",hpd,col)));
			htmpY_final->Add((TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","Final",hpd,col)));
			
			htmp_self->Add((TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","Self",hpd,col)));
			htmpX_self->Add((TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","Self",hpd,col)));
			htmpY_self->Add((TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","Self",hpd,col)));
			
			h_r0->Fill(m_hpdCondDB_new[i][0]);
			h_r1->Fill(m_hpdCondDB_new[i][1]);
			h_r2->Fill(m_hpdCondDB_new[i][2]);
			h_r3->Fill(m_hpdCondDB_new[i][3]);
			h_a0->Fill(m_hpdCondDB_new[i][4]);
			h_a1->Fill(m_hpdCondDB_new[i][5]);
			h_a2->Fill(m_hpdCondDB_new[i][6]);
			h_a3->Fill(m_hpdCondDB_new[i][7]);
			h_r0_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][0]);
			h_r1_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][1]);
			h_r2_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][2]);
			h_r3_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][3]);
			h_a0_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][4]);
			h_a1_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][5]);
			h_a2_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][6]);
			h_a3_2D->SetBinContent(hpd+1,col+1,m_hpdCondDB_new[i][7]);
			
			h_Res_2D->SetBinContent(hpd+1,col+1,((TH1D*)infile->Get(Form("%s_Resolution_HPD%i_COL%i","Final",hpd,col)))->GetMean());
			h_ResX_2D->SetBinContent(hpd+1,col+1,((TH1D*)infile->Get(Form("%s_ResolutionX_HPD%i_COL%i","Final",hpd,col)))->GetRMS());
			h_ResY_2D->SetBinContent(hpd+1,col+1,((TH1D*)infile->Get(Form("%s_ResolutionY_HPD%i_COL%i","Final",hpd,col)))->GetRMS());
			
			isfit=1;
			pr0=m_hpdCondDB_new[i][0];
			pr1=m_hpdCondDB_new[i][1];
			pr2=m_hpdCondDB_new[i][2];
			pr3=m_hpdCondDB_new[i][3];
			pa0=m_hpdCondDB_new[i][4];
			pa1=m_hpdCondDB_new[i][5];
			pa2=m_hpdCondDB_new[i][6];
			pa3=m_hpdCondDB_new[i][7];
			
			pr0_err=m_hpdCondDB_new_err[i][0];
			pr1_err=m_hpdCondDB_new_err[i][1];
			pr2_err=m_hpdCondDB_new_err[i][2];
			pr3_err=m_hpdCondDB_new_err[i][3];
			pa0_err=m_hpdCondDB_new_err[i][4];
			pa1_err=m_hpdCondDB_new_err[i][5];
			pa2_err=m_hpdCondDB_new_err[i][6];
			pa3_err=m_hpdCondDB_new_err[i][7];

			//Grid shifts
			gs0=m_hpdCondDB_new[i][8];
			gs1=m_hpdCondDB_new[i][9];
			gs0_err=m_hpdCondDB_new_err[i][8];
			gs1_err=m_hpdCondDB_new_err[i][9];
			
		}
		else{
			
			isfit=0;
			pr0=m_hpdCondDB_old[i][0];
			pr1=m_hpdCondDB_old[i][1];
			pr2=m_hpdCondDB_old[i][2];
			pr3=m_hpdCondDB_old[i][3];
			pa0=m_hpdCondDB_old[i][4];
			pa1=m_hpdCondDB_old[i][5];
			pa2=m_hpdCondDB_old[i][6];
			pa3=m_hpdCondDB_old[i][7];

			pr0_err=m_hpdCondDB_new_err[i][0];
			pr1_err=m_hpdCondDB_new_err[i][1];
			pr2_err=m_hpdCondDB_new_err[i][2];
			pr3_err=m_hpdCondDB_new_err[i][3];
			pa0_err=m_hpdCondDB_new_err[i][4];
			pa1_err=m_hpdCondDB_new_err[i][5];
			pa2_err=m_hpdCondDB_new_err[i][6];
			pa3_err=m_hpdCondDB_new_err[i][7];

			//Grid shifts
			gs0=m_hpdCondDB_old[i][8];
			gs1=m_hpdCondDB_old[i][9];
			gs0_err=m_hpdCondDB_new_err[i][8];
			gs1_err=m_hpdCondDB_new_err[i][9];

// 			pr0=0;
// 			pr1=0;
// 			pr2=0;
// 			pr3=0;
// 			pa0=0;
// 			pa1=0;
// 			pa2=0;
// 			pa3=0;
		}
		
		outtree->Fill();
	}
	
	outtree->Write();
	
	info() << "NTuple Written!"<<endmsg;
	
	TCanvas *sum_canvas_1=new TCanvas("GlobalSummary_1","GlobalSummary_1",1000,1000);
	sum_canvas_1->Divide(2,2);
	sum_canvas_1->cd(1);
	htmp_final->SetMaximum(TMath::Max(htmp_final->GetMaximum(),TMath::Max(htmp_old->GetMaximum(),htmp_init->GetMaximum())) * 1.1);
	htmp_final->SetMinimum(0);
	htmp_final->Draw();
	htmp_old->Draw("SAMES");
	htmp_init->Draw("SAMES");
	legend->Draw("SAME");
	sum_canvas_1->cd(2);
	htmpX_final->SetMaximum(TMath::Max(htmpX_final->GetMaximum(),TMath::Max(htmpX_old->GetMaximum(),htmpX_init->GetMaximum())) * 1.1);
	htmpX_final->SetMinimum(0);
	htmpX_final->Draw();
	htmpX_old->Draw("SAMES");
	htmpX_init->Draw("SAMES");
	legend->Draw("SAME");
	sum_canvas_1->cd(3);
	htmpY_final->SetMaximum(TMath::Max(htmpY_final->GetMaximum(),TMath::Max(htmpY_old->GetMaximum(),htmpY_init->GetMaximum())) * 1.1);
	htmpY_final->SetMinimum(0);
	htmpY_final->Draw();
	htmpY_old->Draw("SAMES");
	htmpY_init->Draw("SAMES");
	legend->Draw("SAME");
	sum_canvas_1->cd(4);
	htmpY_self->SetMaximum(TMath::Max(htmp_self->GetMaximum(),TMath::Max(htmpX_self->GetMaximum(),htmpY_self->GetMaximum())) * 1.1);
	htmpY_self->SetMinimum(0);
	htmpY_self->Draw();
	htmpX_self->Draw("SAMES");
	htmp_self->Draw("SAMES");
	legend2->Draw("SAME");
	
	sum_canvas_1->Write();
	
	
	TCanvas *sum_canvas_2=new TCanvas("GlobalSummary_2","GlobalSummary_2",1000,1000);
	sum_canvas_2->Divide(2,2);
	sum_canvas_2->cd(1);
	h_Res_2D->Draw("colz");
	sum_canvas_2->cd(2);
	h_ResX_2D->Draw("colz");
	sum_canvas_2->cd(3);
	h_ResY_2D->Draw("colz");
	sum_canvas_2->cd(4);
	
	sum_canvas_2->Write();
	
	TCanvas *sum_canvas_radial_1=new TCanvas("GlobalSummary_radial_1","GlobalSummary_radial_1",1000,1000);
	sum_canvas_radial_1->Divide(2,2);
	sum_canvas_radial_1->cd(1);
	h_r0->Draw();
	sum_canvas_radial_1->cd(2);
	h_r1->Draw();
	sum_canvas_radial_1->cd(3);
	h_r2->Draw();
	sum_canvas_radial_1->cd(4);
	h_r3->Draw();
	sum_canvas_radial_1->Write();
	
	TCanvas *sum_canvas_radial_2=new TCanvas("GlobalSummary_radial_2","GlobalSummary_radial_2",1000,1000);
	sum_canvas_radial_2->Divide(2,2);
	sum_canvas_radial_2->cd(1);
	h_r0_2D->Draw("colz");
	sum_canvas_radial_2->cd(2);
	h_r1_2D->Draw("colz");
	sum_canvas_radial_2->cd(3);
	h_r2_2D->Draw("colz");
	sum_canvas_radial_2->cd(4);
	h_r3_2D->Draw("colz");
	
	sum_canvas_radial_2->Write();
	
	
	TCanvas *sum_canvas_axial_1=new TCanvas("GlobalSummary_axial_1","GlobalSummary_axial_1",1000,1000);
	sum_canvas_axial_1->Divide(2,2);
	sum_canvas_axial_1->cd(1);
	h_a0->Draw();
	sum_canvas_axial_1->cd(2);
	h_a1->Draw();
	sum_canvas_axial_1->cd(3);
	h_a2->Draw();
	sum_canvas_axial_1->cd(4);
	h_a3->Draw();
	sum_canvas_axial_1->Write();
	
	TCanvas *sum_canvas_axial_2=new TCanvas("GlobalSummary_axial_2","GlobalSummary_axial_2",1000,1000);
	sum_canvas_axial_2->Divide(2,2);
	sum_canvas_axial_2->cd(1);
	h_a0_2D->Draw("colz");
	sum_canvas_axial_2->cd(2);
	h_a1_2D->Draw("colz");
	sum_canvas_axial_2->cd(3);
	h_a2_2D->Draw("colz");
	sum_canvas_axial_2->cd(4);
	h_a3_2D->Draw("colz");
	
	sum_canvas_axial_2->Write();
	
	return sc;
}

//==============================================================================
// Fill HPD positions in Panel Ref Frame
//==============================================================================

StatusCode MDMFitParameters::FillHPDWinCentresOnPF(){
	
	StatusCode sc = StatusCode::SUCCESS;
	info() << "Filling HPD windows centres on HPDPlane" << endmsg;
	
// 	std::vector< std::vector< Gaudi::XYZPoint > > output_vector(m_ncols*m_nhpds);
	
// 	panel->initialize();
	unsigned int i=0, k=0;
	for(i=0;i<m_ncols*m_nhpds;i++){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
		debug() << "(Decoded) HPD Code is "<< i << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		int code=this->m_CodeHPD(col ,hpd);
		debug() << "(Code) HPD Code is "<< code << "  HPD index is "<< hpd << "   Column index is " << col << endreq;
		
		const std::string hpdlocation=Form("/dd/Structure/LHCb/BeforeMagnetRegion/Rich1/HPDPanel%i/HPD:%i",m_panel,i+98*m_panel);
		
		DeRichHPD *tmp_hpd= getDet<DeRichHPD>(detSvc(),hpdlocation);
// 		tmp_hpd->updateMagParInMDMSAlgos();
// 		detSvc()->updateObject(tmp_hpd);
// 		for(k=0; k < SmartIDs[i].size(); k++){
			Gaudi::XYZPoint tmppoint;
			tmppoint.SetXYZ(0.0,0.0,0.0);
// 			double x= peaksOnAnode[i][k].x();
// 			double y= peaksOnAnode[i][k].y();
			tmppoint=(tmp_hpd->fromHPDToPanel())*tmppoint;
			m_HPDWindowCentresOnPanelRF[i]=(tmppoint);
// 			std::cout << "i: "<< i << std::endl;
			if(fitted[i])m_HPDWindowPivotOnPanelRF[i]=(tmp_hpd->fromHPDToPanel())*MDMGrids[i]->GetPivotExpectedCoord();
			else m_HPDWindowPivotOnPanelRF[i]=tmppoint;
// 		}
		
		
// 		delete tmp_hpd;
	}
// 	delete panel;
	return sc;
}

//==============================================================================
// Fill HPD Shifts
//==============================================================================
StatusCode MDMFitParameters::FillHPDShifts(TFile *infile){
	StatusCode sc = StatusCode::SUCCESS;
	int i=0,k=0;
	
	this->FillHPDWinCentresOnPF();
	
	TH1F *h_res_X=new TH1F("HPD_Shifts_Distribution_X","HPD Shifts wrt CondDB, X (Panel Coordinates);X_{CondDB}-X_{Fit} (mm);# HPDs",50,-10,10);
	TH1F *h_res_Y=new TH1F("HPD_Shifts_Distribution_Y","HPD Shifts wrt CondDB, Y (Panel Coordinates);Y_{CondDB}-Y_{Fit} (mm);# HPDs",50,-10,10);
	
	//check reference hpd
	if(fitted[m_reference_HPD]){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(m_reference_HPD, &col ,&hpd);
		std::cout << "Calculating shifts wrt HPD" << hpd << " COL"<<col<<std::endl;
	}
	else{
		std::cout << "No fit data available for reference HPD, changing it.." <<std::endl;
		for(i=0;i<m_nhpds*m_ncols;i++){
			if(fitted[i]){
				m_reference_HPD=i;
				int hpd=-1,col=-1;
				this->m_DecodeHPD(m_reference_HPD, &col ,&hpd);
				std::cout << "Calculating shifts wrt HPD" << hpd << " COL"<<col<<std::endl;
				break;
			}
		}
	}
	
	//difference in Y between modules
	double module_diff=69.7;
	double step=m_gstep;
	double phi=MDMGrids[m_reference_HPD]->rot();
	std::cout << "Light Bar Angle: " << phi << std::endl;
	//reference hpd position, col and row
	Gaudi::XYZPoint ref_pos=MDMGrids[m_reference_HPD]->GetPivotExpectedCoord();
	int ref_col=-1;
	int ref_row=-1;
	MDMGrids[m_reference_HPD]->GetPivotRowCol(&ref_row,&ref_col);
	
	TGraph *gRef=new TGraph(m_nhpds*m_ncols);
	gRef->SetName("Fitted_HPD_New_Centres_Positions");
	gRef->SetMarkerStyle(2);
	gRef->SetMarkerColor(2);
	gRef->SetMarkerSize(2);
	int hpd_ref=-1,col_ref=-1;
	this->m_DecodeHPD(m_reference_HPD, &col_ref ,&hpd_ref);
	
	for(i=0;i<m_nhpds*m_ncols;i++){
		if(fitted[i]){
			
			double getexpTube_X=m_HPDWindowCentresOnPanelRF[i].X()-m_HPDWindowCentresOnPanelRF[m_reference_HPD].X();
			double getexpTube_Y=m_HPDWindowCentresOnPanelRF[i].Y()-m_HPDWindowCentresOnPanelRF[m_reference_HPD].Y();
			
			int hpd=-1,col=-1;
			this->m_DecodeHPD(i, &col ,&hpd);
			
			int tmpcol=-1;
			int tmprow=-1;
			MDMGrids[i]->GetPivotRowCol(&tmprow,&tmpcol);
			
// 			MDMGrids[i]->PrintGridsPivot();
// 			std::cout << MDMGrids[i]->GetGridsPivot() << std::endl;
// 			std::cout << tmprow << "\t" << tmpcol << std::endl;
			
			double x_ref=m_HPDWindowPivotOnPanelRF[m_reference_HPD].X();
			double y_ref=m_HPDWindowPivotOnPanelRF[m_reference_HPD].Y();
			
			double bdist=999;
			int thek=0;
			for(k=0;k<200; k++){
				double EX1=((x_ref - module_diff/2*(hpd_ref-hpd+(k-100)) + sin(phi)*step*(tmprow-ref_row) - cos(phi)*step*(tmpcol-ref_col))) - m_HPDWindowPivotOnPanelRF[i].X();
				if(m_panel==0) EX1=((x_ref + module_diff/2*(hpd_ref-hpd+(k-100)) - sin(phi)*step*(tmprow-ref_row) + cos(phi)*step*(tmpcol-ref_col))) - m_HPDWindowPivotOnPanelRF[i].X();
				if(fabs(bdist)>fabs(EX1)) {
					thek=k;
					bdist=EX1;
				}
			}
			double expectedX;
			expectedX=x_ref+sin(phi)*step*(tmprow-ref_row)-cos(phi)*step*(tmpcol-ref_col) - module_diff/2*(hpd_ref-hpd+(thek-100));
			if(m_panel==0) expectedX=x_ref-sin(phi)*step*(tmprow-ref_row)+cos(phi)*step*(tmpcol-ref_col) + module_diff/2*(hpd_ref-hpd+(thek-100));
			double expectedY;
// 			double expectedY=((m_HPDWindowPivotOnPanelRF[m_reference_HPD].Y() + cos(phi)*step*(tmprow-ref_row) - sin(phi)*step*(tmpcol-ref_col)));
			
			expectedY=y_ref+cos(phi)*step*(tmprow-ref_row) +sin(phi)*step*(tmpcol-ref_col);
			if(m_panel==0)expectedY=y_ref-cos(phi)*step*(tmprow-ref_row) -sin(phi)*step*(tmpcol-ref_col);
			
			
// 			if(m_panel==0){
// 				expectedX=-expectedX;
// 				expectedY=-expectedY;
// 			}
			
			gRef->SetPoint(i,(m_HPDWindowPivotOnPanelRF[i].X()-expectedX)+m_HPDWindowCentresOnPanelRF[i].X(),(m_HPDWindowPivotOnPanelRF[i].Y()-expectedY)+m_HPDWindowCentresOnPanelRF[i].Y());
			
			h_res_X->Fill(m_HPDWindowPivotOnPanelRF[i].X()-expectedX);
			h_res_Y->Fill(m_HPDWindowPivotOnPanelRF[i].Y()-expectedY);
// 			std::cout << col << "\t" << hpd << std::endl;
// 			std::cout << expectedX-m_HPDWindowPivotOnPanelRF[i].X() << "\t" << expectedY-m_HPDWindowPivotOnPanelRF[i].Y() << std::endl;
			
// 			gRef->Print();
		}
	}
	
	
	
	TGraph *g=new TGraph(m_nhpds*m_ncols);
	g->SetName("HPD_CondDB_Positions");
	g->SetMarkerStyle(2);
	g->SetMarkerSize(2);
	TGraph *gP=new TGraph(m_nhpds*m_ncols);
	gP->SetName("HPD_Pivot_Positions");
	gP->SetMarkerStyle(2);
	gP->SetMarkerColor(3);
	gP->SetMarkerSize(2);
	
	TGraph *gRefData=new TGraph(m_nhpds*m_ncols);
	gRefData->SetName("TestPoint");
	gRefData->SetMarkerStyle(2);
	gRefData->SetMarkerColor(3);
	gRefData->SetMarkerSize(2);
	
	for(i=0;i<m_nhpds*m_ncols;i++){
		g->SetPoint(i,m_HPDWindowCentresOnPanelRF[i].X(),m_HPDWindowCentresOnPanelRF[i].Y());
		gP->SetPoint(i,m_HPDWindowPivotOnPanelRF[i].X(),m_HPDWindowPivotOnPanelRF[i].Y());
		gRefData->SetPoint(i,m_peaksOnPanel_testing[i].X(),m_peaksOnPanel_testing[i].Y());
	}
	
	TLatex latex;
	latex.SetTextSize(0.015);
	latex.SetTextAlign(13);//align at top
	latex.SetTextColor(2);
	
	
	TCanvas *can_shift=new TCanvas("HPD_Shift_summary","HPD_Shift_summary");
	g->Draw("AP");
// 	gP->Draw("P SAME");
	gRef->Draw("P SAME");
	double *xvec=gRef->GetX();
	double *yvec=gRef->GetY();
	latex.DrawLatex(m_HPDWindowCentresOnPanelRF[m_reference_HPD].X()+1,m_HPDWindowCentresOnPanelRF[m_reference_HPD].Y(),"Ref HPD");
	
	
	
	for(i=0;i<m_nhpds*m_ncols;i++){
		int hpd=-1,col=-1;
		this->m_DecodeHPD(i, &col ,&hpd);
// 		latex.SetTextColor(4);
// 		latex.DrawLatex(m_HPDWindowCentresOnPanelRF[i].X(),m_HPDWindowCentresOnPanelRF[i].Y()+1,Form("C%i,H%i",col,hpd));
// 		latex.SetTextColor(2);
// 		latex.DrawLatex(xvec[i],yvec[i],Form("C%i,H%i",col,hpd));
// 		if(fitted[i]) latex.DrawLatex(xvec[i],yvec[i],Form("%i",MDMGrids[i]->GetPivotStep()));
		TArc *arc=new TArc(m_HPDWindowCentresOnPanelRF[i].X(),m_HPDWindowCentresOnPanelRF[i].Y(),42);
		arc->SetFillStyle(0);
		arc->Draw("SAME");
		TArc *arc2=new TArc(xvec[i],yvec[i],42);
		arc2->SetFillStyle(0);
		arc2->SetLineColor(2);
		arc2->SetLineStyle(2);
		if(fitted[i] && false){
			TGraph *g_panel = (TGraph*)infile->Get(Form("%s_HPD%i_COL%i","peaks_on_panel_final",hpd,col));
			g_panel->SetMarkerSize(0.1);
			g_panel->Draw("P SAME");
			
			
		}
		arc2->Draw("SAME");
	}
	g->Draw("P SAME");
	gRef->Draw("P SAME");
// 	gRefData->Draw("P SAME");
	can_shift->Update();
	can_shift->Write();
	
	TCanvas *can_shift_2=new TCanvas("HPD_Shift_summary_2","HPD_Shift_summary_2");
	can_shift_2->Divide(2,1);
	can_shift_2->cd(1);
	h_res_X->Draw();
	can_shift_2->cd(2);
	h_res_Y->Draw();
	
	can_shift_2->Write();
	
	g->Write();
	gRef->Write();
	
	return sc;
}

