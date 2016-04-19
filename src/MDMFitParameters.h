#ifndef MDMFitParameters_H
#define MDMFitParameters_H 1

#include <cmath>
#include <iostream>
#include <fstream>
// #include "boost/multi_array.hpp"

#include "RichKernel/RichAlgBase.h"
#include "RichKernel/BoostArray.h"

#include "GaudiAlg/GaudiTupleAlg.h"

#include "Kernel/RichSmartID.h"
#include "GaudiKernel/IUpdateManagerSvc.h"
// #include "GaudiKernel/IMagneticFieldSvc.h"
#include "Kernel/ILHCbMagnetSvc.h"
// #include "Kernel/MagnetCondLocations.h"
#include "Event/RichDigit.h"

//RichDet
#include "RichDet/DeRich.h"
#include "RichDet/DeRichSystem.h"
#include "RichDet/DeRichHPD.h"
#include "RichDet/DeRichHPDPanel.h"
#include "RichDet/DeRichRadiator.h"

#include "LHCbMath/LHCbMath.h"
#include "MDMGridModel.h"

using namespace Rich;

class MDMFitParameters : public Rich::AlgBase{
	public:
		MDMFitParameters(const std::string& name, ISvcLocator* pSvcLocator);
		virtual ~MDMFitParameters();
		
		virtual StatusCode initialize();
		virtual StatusCode execute();
		virtual StatusCode finalize();
		
	private:
		unsigned long int m_nEvt;
		unsigned long int m_iteration;
		unsigned int m_nhpds;
		unsigned int m_ncols;
		unsigned int m_npeaks;
		unsigned int m_npeaks_ref;
		int m_bfield;
		int m_bfield_ref;
		int m_panel;
		int m_run;
		double m_defmag;
		double m_defmag_o2;
		double m_defmag_o3;
		double m_gstep;
		std::vector<bool> fitted;
		
		std::vector< TF1* > m_RadialFunctions;
		std::vector< TF1* > m_AxialFunctions;
		std::vector< TH1D* > m_InitResolutions;
		std::vector< TH1D* > m_FinalResolutions;
		std::vector< MDMGridModel* > MDMGrids;
// 		std::vector< TF2* > m_ParFunctions;
		std::vector< LHCb::RichSmartID::Vector > m_SmartIDs;
		std::vector< std::vector< int > > m_calibSteps;
		std::vector< Gaudi::XYZPoint > m_HPDShifts;
		int m_reference_HPD;
		std::vector< Gaudi::XYZPoint > m_HPDCentres;
		std::vector< Gaudi::XYZPoint > m_HPDWindowCentresOnPanelRF;
		std::vector< Gaudi::XYZPoint > m_HPDWindowPivotOnPanelRF;
		std::vector< Gaudi::XYZPoint > m_HPDCentresOnWindow;
		std::vector< Gaudi::XYZPoint > m_HPDWindowPivotOnPanel;
		std::vector< std::vector< Gaudi::XYZPoint > > m_peaksOnAnode;
		std::vector< std::vector< Gaudi::XYZPoint > > m_peaksOnPhotoCathode;
		std::vector< std::vector< Gaudi::XYZPoint > > m_peaksOnPhotoCathode_FIELDON;
		std::vector< std::vector< Gaudi::XYZPoint > > m_peaksOnPhotoCathode_final;
		std::vector< std::vector< Gaudi::XYZPoint > > m_peaksOnPanel_final;
		std::vector< Gaudi::XYZPoint > m_peaksOnPanel_testing;
// 		std::vector< DeRichHPD > m_HPD_list;
		std::string m_infilepath;
		std::string m_outfilepath;
		std::string m_outfiletuple;
		std::string m_treepath;
		//New random seed for systematic jobs (guarantess that U and D boxes get the same value in a given fit)
		int m_randomseed;
		//New radius change factor for systematics jobs
		float m_radiusfactor;
		ILHCbMagnetSvc* m_magFieldSvc;
		// The Update Manager Service
		IUpdateManagerSvc* m_updMgrSvc;
		Condition* m_mag_cond;
		std::vector< std::vector<double> > m_hpdCondDB_old;
		std::vector< std::vector<double> > m_hpdCondDB_new;
		//New MDMS parameter errors vector
		std::vector< std::vector<double> > m_hpdCondDB_new_err;
		
	private: //methods
		int m_CodeHPD(int col, int hpd);
		void m_DecodeHPD(int code, int* col, int* hpd);
		StatusCode getMagFieldSvc();
		void StoreBOFFMagnification();
		StatusCode InitializeSmartIDs();
		std::vector< std::vector< Gaudi::XYZPoint > > MoveToPhotoCathode(std::vector< LHCb::RichSmartID::Vector > SmartIDs, std::vector< std::vector< Gaudi::XYZPoint > > peaksOnAnode, bool iscorrected);
		std::vector< std::vector< Gaudi::XYZPoint > > HPD2Panel(std::vector< LHCb::RichSmartID::Vector > SmartIDs, std::vector< std::vector< Gaudi::XYZPoint > > peaksOnWindow, bool iscorrected);
		
		std::vector< TGraph* > FillTGraphs(std::string Prefix, int color, std::vector< std::vector< Gaudi::XYZPoint > > peaks);
// 		StatusCode CreateHPDList();
		StatusCode UpdateCondDB();
		StatusCode ChangeMagneticField(int polarity, double current);
		StatusCode upd_mag(){return StatusCode::SUCCESS;};
		std::string printMagnetConditions();
		StatusCode fromTF1toDoubleVector(int hpdnumber, TF1 *radialF, TF1* axialF, bool includemag);
		StatusCode FillParamHistos();
		StatusCode FitGrid(TFile *infile);
		StatusCode FillResolutions(TString Prefix, std::vector< std::vector< Gaudi::XYZPoint > > peaks);
		StatusCode SummaryPlots(TFile *infile);
		StatusCode GlobalSummary(TFile *infile);
		StatusCode FillHPDShifts(TFile *infile);
		StatusCode FillHPDWinCentresOnPF();
};
#endif
