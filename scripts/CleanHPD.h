//HPD class for MDMS cleaning
//implementd by Andrea Contu 6 Sept 2010, acontu@cern.ch

#ifndef __CleanHPD__
#define __CleanHPD__

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "TStyle.h"
#include "TROOT.h"
#include "TObject.h"
#include "TNamed.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TRegexp.h"
#include "TString.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TList.h"
#include "TMap.h"
#include "TObjArray.h"

using namespace std;

class CleanPeak : public TNamed{
private:
	int m_step; //calibration step
	int m_box; //box
	int m_column; //column
	int m_hpd; //hpd
	int m_npes; //number of hits in the peak
	int m_peak_class;
	double m_max; //peak height
	double m_prop; //peak proportion wrt total hits
	double m_meanX; //mean peak position from hits, X
	double m_meanY; //mean peak position from hits, Y
	double m_rmsX; //rms X
	double m_rmsY; //rms Y
	double m_sigmaX; //2D gaussian fit rms X
	double m_sigmaY; //2D gaussian fit rms Y
	double m_muX; //2D fitted peak position X
	double m_muY; //2D fitted peak position Y
	double m_rot; //peak rotation
	double m_chi2; //fit chi2
	TString m_pass;
	bool m_is_single;
	double m_like1;
	double m_like2;
	double m_like3;
	double m_distvar;
	double m_xdiff;
	
public:
	CleanPeak(){};
	CleanPeak(int step, int box, int column, int hpd, int npes=-1, int peak_class=-1, double max=-1, double prop=-1, double meanX=-1, double meanY=-1, double rmsX=-1, double rmsY=-1, double sigmaX=-1, double sigmaY=-1, double muX=-1, double muY=-1, double rot=-1, double chi2=-1){this->SetPeak(step, box, column, hpd, npes, peak_class, max, prop, meanX, meanY, rmsX, rmsY, sigmaX, sigmaY, muX, muY, rot, chi2);};
	bool SetPeak(int step, int box, int column, int hpd, int npes, int peak_class, double max, double prop, double meanX, double meanY, double rmsX, double rmsY, double sigmaX, double sigmaY, double muX, double muY, double rot, double chi2){
		m_step=step;
		m_box=box;
		m_column=column;
		m_hpd=hpd;
		m_npes=npes;
		m_peak_class=peak_class;
		m_max=max;
		m_prop=prop;
		m_meanX=meanX;
		m_meanY=meanY;
		m_rmsX=rmsX;
		m_rmsY=rmsY;
		m_sigmaX=sigmaX;
		m_sigmaY=sigmaY;
		m_muX=muX;
		m_muY=muY;
		m_rot=rot;
		m_chi2=chi2;
		m_is_single=false;
		m_distvar=0.;
		return true;
	};
	void SetDistVar(double distvar){m_distvar=distvar;};
	double GetDistVar(){return m_distvar;};
	void SetXDiff(double xdiff){m_xdiff=xdiff;};
	double GetXDiff(){return m_xdiff;};
	void SetSingle(bool single){m_is_single=single;};
	bool IsSingle(){return m_is_single;};
	void SetPassLevel(TString passstring){m_pass=passstring;};
	void SetLike1(double like){m_like1=like;};
	void SetLike2(double like){m_like2=like;};
	void SetLike3(double like){m_like3=like;};
	TString GetPassLevel(){return m_pass;};
	double GetLike1(){return m_like1;};
	double GetLike2(){return m_like2;};
	double GetLike3(){return m_like3;};
	int GetCalibStep(){return m_step;};
	int GetBox(){return m_box;};
	int GetColumn(){return m_column;};
	int GetHPD(){return m_hpd;};
	int GetNPEs(){return m_npes;};
	int GetClass(){return m_peak_class;};
	double GetMax(){return m_max;};
	double GetProp(){return m_prop;};
	double GetMeanX(){return m_meanX;};
	double GetMeanY(){return m_meanY;};
	double GetRMSX(){return m_rmsX;};
	double GetRMSY(){return m_rmsY;};
	double GetSigmaX(){return m_sigmaX;};
	double GetSigmaY(){return m_sigmaY;};
	double GetSigma(){return TMath::Sqrt(m_sigmaX*m_sigmaX/2+m_sigmaY*m_sigmaY/2);};
	double GetMuX(){return m_muX;};
	double GetMuY(){return m_muY;};
	double GetRotation(){return m_rot;};
	double GetChi2(){return m_chi2;};
	
// 	ClassDef(CleanPeak,1)
};


class CleanHPD : public TObject{
private:
	TList m_hitmap;
	int m_col;
	int m_hpd;
	int m_box;
	int m_npeaks;
	double m_centerX;
	double m_centerY;
	double m_meanX;
	double m_meanY;
	double m_fitradius;
	bool m_fill_map(TTree *tree);
	TString m_pass;
public:
	CleanHPD(){m_box=-1; m_col=-1; m_hpd=-1;};
	CleanHPD(int box, int column, int hpd){m_box=box; m_col=column; m_hpd=hpd;};
	CleanHPD(int box, int column, int hpd, TTree *tree){m_box=box; m_col=column; m_hpd=hpd; this->m_fill_map(tree);};
	~CleanHPD(){};
	TList* GetHitmap(){return &m_hitmap;};
	bool FitHPDCenter();
	void SetPassLevel(TString pass){m_pass=pass;};
	TString GetPassLevel(){return m_pass;};
	double GetCenterX(){return m_centerX;}
	double GetCenterY(){return m_centerY;}
	double GetMeanX(){return m_meanX;}
	double GetMeanY(){return m_meanY;}
	double GetFittedRadius(){return m_fitradius;}
	int GetBox(){return m_box;};
	int GetColumn(){return m_col;};
	int GetHPD(){return m_hpd;};
	int GetNPeaks(){return m_npeaks;};
	bool FillMap(TTree *tree){return this->m_fill_map(tree);};
	TGraph* DrawPeakMap();
	TH2D* DrawPeakHisto();
// 	ClassDef(CleanHPD,1)
};

#endif
