#ifndef MDMGridModel_H
#define MDMGridModel_H 1

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "gsl/gsl_math.h"
#include "GaudiKernel/Point3DTypes.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TLine.h"
#include "TProfile.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TVirtualFitter.h"
//#include "TFitterMinuit.h"
#include "TArc.h"
//minuit
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FCNBase.h"
// using namespace std;

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag>::Scalar geoScalar;

class MDMGridModel {
	
	public:
		MDMGridModel();
		MDMGridModel(TString name, std::vector< Gaudi::XYZPoint > pts, std::vector<int> steps, int scalef, double angle, double width, int hpd, int column);
		~MDMGridModel();
		
		double centreX(){return cX;};
		double centreY(){return cY;};
		double width(){return step;};
		double rot(){return phi;};
		double hpd(){return hpd_number;};
		double col(){return col_number;};
		
		int pivotClass(){return minstep%(scalefactor+1);};
		
		void set_centreX(double val){cX=val;}
		void set_centreY(double val){cY=val;}
		void set_width(double val){step=val;}
		void set_rot(double val){phi=val;}
		Gaudi::XYZPoint PeaksCMS;
		bool GetResolution(std::vector< Gaudi::XYZPoint > datapoints,TH1D* h, TString Axis);
		bool GetSelfResolution(TH1D* h, TString Axis);
		std::vector<double> FitPoints(std::vector< Gaudi::XYZPoint > datapoints, std::vector<int> steps, int scalef);
		void DrawGrid(TString options, std::vector< Gaudi::XYZPoint > datapoints, double Xmin, double Xmax, double Ymin, double Ymax);
		bool FitParameters();
		std::vector<double> GetGridShift(){return gShift;};
		std::vector<double> GetGridShift_err(){return gShift_err;};
		std::vector<double> GetParameters(){return parameters;};
		std::vector<double> GetParameters_err(){return parameters_err;};
		double modules_row_separation;
		double modules_col_separation;
		int GetGridsPivot(){return gridspivot;};
		void PrintGridsPivot(){std::cout << gridspivot<< std::endl;};
		int GetPivotStep(){return st[gridspivot];};
		Gaudi::XYZPoint GetPivotCoord();
		Gaudi::XYZPoint GetPivotExpectedCoord();
		void GetPivotRowCol(int *row, int *col){this->GetRowAndColumnFromLOGIC(st[gridspivot], row, col);};
		int GetTopGrid(){return topgrid;};
		
	private: //data
		double cX;
		double cY;
		double step;
		double phi;
		int hpd_number;
		int col_number;
		bool pvf;
		TString gname;
		std::vector<int> st;
		
		std::vector<int> pong_row;
		std::vector<int> pong_col;
		std::vector< Gaudi::XYZPoint > points;
		std::vector< Gaudi::XYZPoint > pointsOnGrids;
		std::vector<TString> pattern_type;
		int scalefactor;
		int tot_cols;
		int tot_rows;
		int minstep;
		int maxstep;
		std::vector< std::vector<int> > LOGIC;
		std::vector<int> LOGIC_Shift;
		//grid parameters
		double cX_1; //grid 1 pivot X
		double cY_1; //grid 1 pivot Y
		double cX_2; //grid 2 pivot X
		double cY_2; //grid 2 pivot Y
		int lowestpeak_1;
		int lowestpeak_2;
		int gridspivot;
		//classified peaks
		int topgrid;
		
		std::vector<double> gShift;
		std::vector<double> gShift_err;
		std::vector<double> parameters;
		std::vector<double> parameters_err;
		
		std::vector<double> gShift_pre;
		std::vector<double> gShift_pre_err;
		std::vector<double> parameters_pre;
		std::vector<double> parameters_pre_err;
// 		std::vector< Gaudi::XYZPoint > points_1;
// 		std::vector< Gaudi::XYZPoint > points_2;
// 		std::vector< Gaudi::XYZPoint > pointsOnGrid_1;
// 		std::vector< Gaudi::XYZPoint > pointsOnGrid_2;
		
	private: //methods
		void Initialise();
		void FlagSteps(std::vector< Gaudi::XYZPoint > datapoints, std::vector<int> steps);
		void DrawLines(TString options, bool orientation, double Xmin, double Xmax, double Ymin, double Ymax);
		void PreparePeaksAndGrids();
		Gaudi::XYZPoint GetRelativePointOnGrid(int cstep, int cstep_ref, Gaudi::XYZPoint point);
		Gaudi::XYZPoint GetClosestGridPoint(Gaudi::XYZPoint point, int k);
		std::vector< std::vector<int> > BuildMDMSLogic(); //this is insane, think of cooler ways of doing it
		void PrintLOGIC(int minstep=-1, int maxstep=-1);
		void GetRowAndColumnFromLOGIC(int pattern, int *row, int *col);
		Gaudi::XYZPoint peaksCMS(std::vector< Gaudi::XYZPoint > pts);
		void PreFit();
		void IteratedFit();
};

#endif
