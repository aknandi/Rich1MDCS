#include "TStyle.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TText.h"
#include "TColor.h"
#include "TLatex.h"

#include <iomanip>
#include <iostream>
using namespace std;

TPaveText *myName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
TPaveText *myUni = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
TText *myLabel = new TText();
TLatex *myLatex = new TLatex();

void mystyle(){
	cout << "**************************" << endl
	     << "*  ROOT Scrilli's Style  *" << endl
	     << "**************************" << endl;
	
// 	TStyle *myStyle= new TStyle("myStyle","Scrilli's official plots style");
	TStyle *myStyle=(TStyle*)gROOT->GetStyle("Plain");
	myStyle->SetName("myStyle");
	myStyle->SetTitle("Scrilli's official plots style");
	// use helvetica-bold-r-normal, precision 2 (rotatable)
	Int_t myFont = 62;
	// line thickness
	Int_t myWidth = 1;

	// use plain black on white colors
	myStyle->SetFrameBorderMode(0);
	myStyle->SetCanvasBorderMode(0);
	myStyle->SetPadBorderMode(0);
	myStyle->SetPadColor(0);
	myStyle->SetCanvasColor(0);
	myStyle->SetStatColor(0);
	myStyle->SetPalette(1);
	myStyle->SetTitleColor(1);
// 	myStyle->SetFillColor(0);

	// set the paper & margin sizes
// 	myStyle->SetPaperSize(20,26);
// 	myStyle->SetPadTopMargin(0.05);
	myStyle->SetPadRightMargin(0.12); // increase for colz plots!!
	myStyle->SetPadBottomMargin(0.12);
	myStyle->SetPadLeftMargin(0.12);

	// use large fonts
	myStyle->SetTextFont(myFont);
	myStyle->SetTextSize(0.05);
	myStyle->SetLabelFont(myFont,"x");
	myStyle->SetLabelFont(myFont,"y");
	myStyle->SetLabelFont(myFont,"z");
	myStyle->SetLabelSize(0.05,"x");
	myStyle->SetLabelSize(0.05,"y");
	myStyle->SetLabelSize(0.05,"z");
	myStyle->SetTitleFont(myFont,"x");
	myStyle->SetTitleFont(myFont,"y");
	myStyle->SetTitleFont(myFont,"z");
	myStyle->SetTitleSize(0.05,"x");
	myStyle->SetTitleSize(0.05,"y");
	myStyle->SetTitleSize(0.05,"z");
	myStyle->SetTitleOffset(1.1,"x");
	myStyle->SetTitleOffset(1.1,"y");
	myStyle->SetTitleOffset(1.1,"z");

	// use bold lines and markers
	myStyle->SetLineWidth(myWidth);
	myStyle->SetFrameLineWidth(myWidth);
	myStyle->SetHistLineWidth(myWidth);
	myStyle->SetFuncWidth(myWidth+1);
	myStyle->SetGridWidth(myWidth);
	myStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	myStyle->SetMarkerStyle(1);
	myStyle->SetMarkerSize(1);

	// by default, do not display histogram decorations:
	myStyle->SetOptStat(0);  
// 	myStyle->SetOptStat(1111);
	myStyle->SetOptTitle(0);
// 	myStyle->SetOptFit(0);
	myStyle->SetOptFit(1111); // show probability, parameters and errors

	// look of the statistics box:
	myStyle->SetStatColor(0);
// 	myStyle->SetStatTextColor(4);
	myStyle->SetStatBorderSize(1);
	myStyle->SetStatFont(myFont);
	myStyle->SetStatFontSize(0.03);
	myStyle->SetStatX(0.9);
	myStyle->SetStatY(0.9);
	myStyle->SetStatW(0.20);
	myStyle->SetStatH(0.15);

	// put tick marks on top and RHS of plots
	myStyle->SetPadTickX(1);
	myStyle->SetPadTickY(1);
	myStyle->SetTickLength(0.02,"X");
	myStyle->SetTickLength(0.02,"Y");
	myStyle->SetTickLength(0.02,"Z");

	// histogram divisions: only 5 in x to avoid label overlaps
	myStyle->SetNdivisions(510,"x");
	myStyle->SetNdivisions(510,"y");

// 	gROOT->ForceStyle();
	//TPaveText *zName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
	myName->SetFillColor(0);
// 	myName->SetTextAlign(12);
	myName->SetBorderSize(0);
	myName->SetTextSize(0.04);
	myName->AddText("Andrea Contu");
	myUni->SetFillColor(0);
// 	myName->SetTextAlign(12);
	myUni->SetBorderSize(0);
	myUni->SetTextSize(0.04);
	myUni->AddText("University of Oxford");

	//TText *zLabel = new TText();
	myLabel->SetTextFont(myFont);
	myLabel->SetTextColor(1);
	myLabel->SetTextSize(0.04);
	myLabel->SetTextAlign(12);

	//TLatex *zLatex = new TLatex();
	myLatex->SetTextFont(myFont);
	myLatex->SetTextColor(1);
	myLatex->SetTextSize(0.04);
	myLatex->SetTextAlign(12);
	gROOT->SetStyle("myStyle");
// 	gROOT->SetStyle("Plain");
	//continous colour palette
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	myStyle->SetNumberContours(NCont);
//	myStyle->SetPaintTextFormat("5.2f");
	myStyle->SetPalette(1);
	gROOT->SetStyle("myStyle");
}
