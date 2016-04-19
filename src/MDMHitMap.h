// $Id: MDMHitMap.h,v 1.1.1.1 2009/04/14 16:06:30 mjohn Exp $
#ifndef MDMHITMAP_H 
#define MDMHITMAP_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"

#include "Kernel/RichSmartID.h"
#include "Event/RichDigit.h"

// RichDet
#include "RichDet/DeRich.h"
#include "RichDet/DeRichHPD.h"
#include "RichDet/DeRichHPDPanel.h"

#include "LHCbMath/LHCbMath.h"
#include "gsl/gsl_math.h"

// Root
#include "TH2F.h"
#include "TCanvas.h"

const unsigned int NHPDPANELS       =  2;
const unsigned int NHPDCOLUMNS      =  7;
const unsigned int NHPDSPERCOLUMN   = 14;
const unsigned int NPIXELCOLUMNS    = 32;
const unsigned int NPIXELSPERCOLUMN =256;

typedef std::pair<double,double> PixelCoordinate;
typedef std::pair<int,PixelCoordinate> PixelHit;

typedef std::vector< std::pair<int,int> > PixelList;

class Peak{
public:
  Peak(double px,double py,double sx,double sy,int n, double p,int cl=1);
  void storeFitResult(double,double,double,double,double,double,double);
  double maximum(){return m_max;}
  double meanX(){return m_meanX;}
  double meanY(){return m_meanY;}
  double rmsX(){return m_rmsX;}
  double rmsY(){return m_rmsY;}
  int nHits(){return m_N;}
  double proportion(){return m_prop;}
  int    classOfPeak(){return m_class;}
  double muX(){return m_muX;}
  double muY(){return m_muY;}
  double sigmaX(){return m_sigmaX;}
  double sigmaY(){return m_sigmaY;}
  double rotation(){return m_rotation;}
  void closest(Peak* other){m_closest=other;}//std::cout<<this<<" is registering "<<m_closest<<std::endl;}
  Peak* closest(){/*std::cout<<this<<" is returning "<<m_closest<<std::endl;*/return m_closest;}
  double chi2(){return m_chi2;}
  bool change(double);
private:
  double m_meanX;
  double m_meanY;
  double m_muX;
  double m_muY;
  double m_sigmaX;
  double m_sigmaY;
  double m_rmsX;
  double m_rmsY;
  int m_N;
  int m_class;
  double m_prop;
  double m_rotation;
  double m_chi2;
  double m_max;
  Peak* m_closest;
};
typedef std::vector<Peak> Peaks;

inline Peak::Peak(double px,double py,double sx,double sy,int n, double p,int cl)
{
  m_meanX=px;
  m_meanY=py;
  m_rmsX=sx;
  m_rmsY=sy;
  m_N=n;
  m_prop=p;
  m_class=cl;
  m_muX=-999;
  m_muY=-999;
  m_sigmaX=-999;
  m_sigmaY=-999;
  m_rotation=-999;
  m_closest=NULL;
  m_max=0;
}

inline void Peak::storeFitResult(double maximum, double mX,double mY,double sX,double sY,double r,double c)
{
  m_muX=mX;
  m_muY=mY;
  m_sigmaX=sX;
  m_sigmaY=sY;
  m_rotation=r;
  m_chi2=c;
  m_max=maximum;
}

inline bool Peak::change(double evidenceContour)
{
  if( fabs(m_closest->meanX()-m_meanX)>m_rmsX * evidenceContour ||
      fabs(m_closest->meanY()-m_meanY)>m_rmsY * evidenceContour ){
    return true;
  }
  return false;
}

class MDMHitMap{
public:
  MDMHitMap();
  ~MDMHitMap(){}
  void addHit(unsigned int,unsigned int);
  void suppressNoisyPixels(int);
  int nHits(){return (int)map->GetEntries();}
  Peaks peaks(bool=false);
  void doClustering(bool=false);
  void setName(const char*, const char*);
  const char* getName(){return map->GetName();}
  const char* getTitle(){return map->GetTitle();}
  int nXBins(){return map->GetNbinsX();}
  int nYBins(){return map->GetNbinsY();}
  int binContent(int i, int j){return (int)map->GetBinContent(i,j);}
  bool badPixel(unsigned int c,unsigned int r);
  void clear();
  void tidy();
  void draw(bool b=true){DRAW=b;}
  void howManyEvents(int n){m_nEvents=n;}
  int howManyEvents(){return m_nEvents;}
  PixelList hotPixels(){return m_hotPixelList;}
private:
  double piBy4;
  TCanvas *can;
  TH2F *map;
  TH2F *noisyMap;
  Peaks m_peaks;
  bool m_doneFitting;
  bool m_doneClustering;
  unsigned int m_nEvents;
  bool DRAW;
  PixelList m_hotPixelList;
};

#endif // MDMHITMAP.H
