#include "TF2.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMarker.h"
#include <boost/lexical_cast.hpp>
#include "MDMHitMap.h"

#include "RichKernel/RichDAQDefinitions.h"

#include "GaudiKernel/PhysicalConstants.h"

double gauss2D(double *x, double *par) {
  double a = par[5];
  double xUnrot = x[0]-par[1];
  double yUnrot = x[1]-par[3];
  double xRotated = xUnrot*cos(a)+yUnrot*sin(a);
  double yRotated =-xUnrot*sin(a)+yUnrot*cos(a);
  double z1 = (0==par[2]) ? 0 : xRotated/par[2];
  double z2 = (0==par[4]) ? 0 : yRotated/par[4];
  return par[0]*exp(-0.5*(z1*z1+z2*z2));
}

MDMHitMap::MDMHitMap()
{
  piBy4=Gaudi::Units::halfpi/4.;
  m_doneClustering=false;
  m_doneFitting=false;
  DRAW=false;
  can=0;
}

void MDMHitMap::setName(const char *name,const char *title="Hitmap")
{
  map = new TH2F(name,title,32,-0.5,31.5,256,-0.5,255.5);
  map->SetStats(false);
  std::string newname = boost::lexical_cast<std::string>(map->GetName())+"_Noisy";
  noisyMap = new TH2F(newname.c_str(),title,32,-0.5,31.5,256,-0.5,255.5);
}

void MDMHitMap::addHit(unsigned int c,unsigned int r)
{
  m_doneFitting = false;
  m_doneClustering = false;
  if(badPixel(c,r)) return;
  map->Fill(c,r);
}

void MDMHitMap::suppressNoisyPixels(int nmax)
{
  for(int count=0;count<nmax;count++){
    int i,j,k;
    map->GetBinXYZ(map->GetMaximumBin(),i,j,k);
    int thisBinContent = (int)map->GetBinContent(i,j);
    if(thisBinContent<0.1*m_nEvents) continue;
    if(((0==map->GetBinContent(i+1,j  )) ||
        (0==map->GetBinContent(i+1,j+1)) ||
        (0==map->GetBinContent(i  ,j+1)) ||
        (0==map->GetBinContent(i-1,j+1)) ||
        (0==map->GetBinContent(i-1,j  )) ||
        (0==map->GetBinContent(i-1,j-1)) ||
        (0==map->GetBinContent(i-1,j  )) ||
        (0==map->GetBinContent(i-1,j+1)))&&
       (thisBinContent>100)){
      map->SetBinContent(i,j,0);
      map->SetEntries(map->GetEntries()-thisBinContent-1);
      noisyMap->SetBinContent(i,j,1);
      m_hotPixelList.push_back(std::pair<int,int>(i-1,j-1));
    }else{
      count=nmax;
    }
  }
}

bool MDMHitMap::badPixel(unsigned int c,unsigned int r)
{
  return (bool)noisyMap->GetBinContent(c+1,r+1);
}

Peaks MDMHitMap::peaks(bool fit)
{
  if(fit){
    if(!m_doneFitting) doClustering(fit);
  }else{
    if(!m_doneClustering) doClustering();
  }
  return m_peaks;
}

void MDMHitMap::clear(){
  m_hotPixelList.clear();
  m_doneClustering=false;
  m_doneFitting=false;
  noisyMap->Reset();
  m_peaks.clear();
  map->Reset();
}

void MDMHitMap::tidy(){
  delete map; 
  delete noisyMap;
  if(can) delete can; 
}

void MDMHitMap::doClustering(bool fit)
{
  m_peaks.clear();
  suppressNoisyPixels(16);
  if(nHits()<200 || nHits()<0.2*m_nEvents){
    return;
  }

  //Rebin to LHCb mode for peak-finding
  std::string newname = boost::lexical_cast<std::string>(map->GetName())+"_Rebinned";
  TH2F *rebinnedMap = (TH2F*)map->RebinY(Rich::DAQ::NumAlicePixelsPerLHCbPixel,newname.c_str());

  //Loop over all bins searching for pixels that have a good proportion of the total hits
  for(int iBin=1;iBin<=rebinnedMap->GetNbinsX();iBin++){
    for(int jBin=1;jBin<=rebinnedMap->GetNbinsY();jBin++){
      int thisBinContent = (int)rebinnedMap->GetBinContent(iBin,jBin);
      if(thisBinContent<(0.01*map->GetEntries())) continue;

      bool isAPeak=true;
      int  peakClass=1;

      //Check that there is no, already-found peak next door
      for(Peaks::iterator ipeak=m_peaks.begin();ipeak!=m_peaks.end();ipeak++){
        if((fabs((iBin-1)-(*ipeak).meanX())<1.5)&&
           (fabs((jBin-1)-(*ipeak).meanY())<1.5)){isAPeak=false;}
      }
      if(!isAPeak) continue;
      
      //Figure out which pixels are close by
      std::vector<PixelHit> nearestPixels;
      std::vector<PixelHit> nextToNearestPixels;
      std::vector<PixelHit> nextToNextToNearestPixels;
      for(int i=iBin-3;i<=iBin+3;i++){
        for(int j=jBin-3;j<=jBin+3;j++){
          if(i==iBin&&j==jBin) continue;
          if(i<1||i>32||j<1||j>32) continue;
          float x = fabs(i-iBin);
          float y = fabs(j-jBin);
          float r = sqrt(x*x+y*y);
          PixelHit p((int)rebinnedMap->GetBinContent(i,j),PixelCoordinate(x,y));
          if(r<2){
            nearestPixels.push_back(p);
          }else if(r<3){
            nextToNearestPixels.push_back(p);
          }else if(r<4){
            nextToNextToNearestPixels.push_back(p);
          }
        }
      }

      //Check the neighbouring pixels to see if they are higher
      for(unsigned int k=0;k<nearestPixels.size();k++){
        if(nearestPixels[k].first>thisBinContent) isAPeak=false;
      }
      if(!isAPeak) continue;

      for(unsigned int k=0;k<nextToNearestPixels.size();k++){
        if(nextToNearestPixels[k].first>thisBinContent) isAPeak=false;
      }
      if(!isAPeak) continue;

      for(unsigned int k=0;k<nextToNextToNearestPixels.size();k++){
        if(nextToNextToNearestPixels[k].first>thisBinContent) peakClass=2;
      }
      if(!isAPeak) continue;

      //Calculate the approx. RMS of the peak
      double xRms=0; double yRms=0; int N=thisBinContent;
      for(unsigned int k=0;k<nearestPixels.size();k++){
        float x = nearestPixels[k].second.first;
        float y = nearestPixels[k].second.second;
        int n = nearestPixels[k].first;
          xRms += x*x*n;
          yRms += y*y*n;
          N += n;
      }
      for(unsigned int k=0;k<nearestPixels.size();k++){
        float x = nextToNearestPixels[k].second.first;
        float y = nextToNearestPixels[k].second.second;
        int n = nextToNearestPixels[k].first;
        xRms += x*x*n;
        yRms += y*y*n;
        N += n;
      }
      if(peakClass==1){
        for(unsigned int k=0;k<nearestPixels.size();k++){
          float x = nextToNextToNearestPixels[k].second.first;
          float y = nextToNextToNearestPixels[k].second.second;
          int n = nextToNextToNearestPixels[k].first;
          xRms += x*x*n;
          yRms += y*y*n;
          N += n;
        }
      }
      xRms = sqrt(xRms/N);
      yRms = sqrt(yRms/N);

      double xMean=iBin-1;
      double yMean=jBin-1;

      //Find the highest ALICE pixel
      int nZeroBins=0;
      int highestBin=0;      
      float highestValue=0;
      for(int j=(jBin-1)*8;j<=jBin*8+1;j++){//Look across 10 ALICE pixels (1+8+1)
        float aliceBinValue = map->GetBinContent(iBin,j);
        if(aliceBinValue>highestValue){
          highestValue=aliceBinValue;
          highestBin=j;
        }else if(0==aliceBinValue){
          nZeroBins++;
        }
      }        
      double yMeanALICE=highestBin-1;

      if(nZeroBins>2) continue;//Good protection against a high number of evenly distributed noisy channels

      Peak peak(xMean,yMean,xRms,yRms,N,float(N)/map->GetEntries(),peakClass);

      if(fit){

        double xSigma=1.0;
        double ySigma=8.0;
        double xRangeLo=xMean-1.5*xSigma;
        double xRangeHi=xMean+1.5*xSigma;
        double yRangeLo=yMeanALICE-1.5*ySigma;
        double yRangeHi=yMeanALICE+1.5*ySigma;

        double iniParams[6]= {   100,    xMean,    xSigma, yMeanALICE,    ySigma,   0  };
        double hiParams[6] = {100000, xRangeHi,  5*xSigma, yRangeHi  ,  5*ySigma, piBy4};
        double loParams[6] = {     0, xRangeLo,.01*xSigma, yRangeLo  ,.01*ySigma,-piBy4};
        TF2 *func = new TF2("func",gauss2D,xRangeLo,xRangeHi,yRangeLo,yRangeHi,6);
        for(int i=0; i<6;i++){
          func->SetParameter(i,iniParams[i]);
          func->SetParLimits(i,loParams[i],hiParams[i]);
        }
        map->Fit(func,"R0NQ");
        
        float maximum = func->GetParameter(0);
        float posX = func->GetParameter(1);
        float posY = (func->GetParameter(3)+0.5)/8.-0.5;
        float sigX = fabs(func->GetParameter(2));
        float sigY = fabs(fabs(func->GetParameter(4))/8.);
        float aRot = func->GetParameter(5);
        float chi2 = func->GetChisquare();
        float nDF  = func->GetNDF();
        peak.storeFitResult(maximum,posX,posY,sigX,sigY,aRot,chi2/nDF);
        
        delete func;
      }
      m_peaks.push_back(peak);
    }
  }
  if(fit) m_doneFitting=true;
  m_doneClustering=true;
  delete rebinnedMap;

  if(m_peaks.size()>0 && DRAW){
    if(!can) can = new TCanvas(map->GetName(),map->GetTitle(),500,500);
    can->cd();
    map->Draw("colz");
    for(Peaks::iterator ipeak=m_peaks.begin();ipeak!=m_peaks.end();ipeak++){
      Peak peak = (*ipeak);
      double textX =(peak.meanX()<22?peak.meanX()+1:22);
      double textY =(peak.meanY()+0.5)*8-0.5;
      TMarker *p1 = new TMarker(peak.meanX(),(peak.meanY()+0.5)*8-0.5,4);
      p1->SetMarkerSize(2*p1->GetMarkerSize());
      p1->SetMarkerColor(10);
      p1->Draw();
      TMarker *p2 = new TMarker(peak.meanX(),(peak.meanY()+0.5)*8-0.5,4);
      p2->SetMarkerSize(1.7*p2->GetMarkerSize());
      p2->Draw();
      TLatex *t1 = new TLatex(textX,textY-5,Form("rms=%3.1f %i %2.0f%% %s",
                                                 (peak.rmsX()+peak.rmsY())/2.,peak.nHits(),peak.proportion()*100,(peak.classOfPeak()==2?"2nd":"1st")));
      t1->SetTextSize(0.5*t1->GetTextSize());
      t1->Draw();
      if(fit){
        TMarker *p = new TMarker(peak.muX(),(peak.muY()+0.5)*8-0.5,20);
        p->SetMarkerColor(10);
        p->Draw();
        TMarker *p0 = new TMarker(peak.muX(),(peak.muY()+0.5)*8-0.5,20);
        p0->SetMarkerSize(0.7*p0->GetMarkerSize());
        p0->Draw();
        TLatex *t0 = new TLatex(textX,textY+5,Form("#sigma: %3.1f pixels, #theta: %2.1f^{#circ}",
                                                            (peak.sigmaX()+peak.sigmaY())/2.,peak.rotation()*45./piBy4));
        t0->SetTextSize(0.5*t0->GetTextSize());
        t0->Draw();  
      }
    }
    can->Update();
  }
}
