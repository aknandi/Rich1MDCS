#include "GaudiKernel/ToolFactory.h"
#include <boost/lexical_cast.hpp>

// Gaudi
#include "GaudiUtils/Aida2ROOT.h"

// ROOT
#include "TApplication.h"

// local
#include "MDMHitMapManager.h"

// LHCb
#include "Event/ODIN.h"

// Rich
#include "RichUtils/RichDAQDefinitions.h"

DECLARE_TOOL_FACTORY( MDMHitMapManager )

// Standard constructor
MDMHitMapManager::MDMHitMapManager( const std::string& type,
                  const std::string& name,
                  const IInterface* parent )
  : TupleToolBase ( type, name, parent )
  , m_eventsSinceEndOfPreviousPattern(0)
  , m_hitMap()
  , m_cache()
  , m_pattern(0)
  , DEBUG(false)
  , FIT(true)
{
  // Define interface
  declareInterface<IHitMapManagerTool>(this);
}

StatusCode MDMHitMapManager::initialize()
{
  // Sets up various tools and services
  const StatusCode sc = TupleToolBase::initialize();
  return sc;
}

void MDMHitMapManager::setup(int n, bool b, bool storehistos)
{
  m_nEvtsPerStep=n;
  DEBUG=b;
	m_StoreHistos=storehistos;
  for(unsigned int a=0;a<NHPDPANELS;a++){
    for(unsigned int b=0;b<NHPDCOLUMNS;b++){
      for(unsigned int c=0;c<NHPDSPERCOLUMN;c++){
        std::string s_b = boost::lexical_cast<std::string>(b);
        std::string s_c = boost::lexical_cast<std::string>(c);
        if(c<10) s_c = "0"+s_c;
        std::string name="UpperBox_Col"+s_b+"_HPD"+s_c;
        if((bool)a) name="LowerBox_Col"+s_b+"_HPD"+s_c;
        std::string mapname = "MainMap_"+name;
        m_hitMap[a][b][c].setName(mapname.c_str(),mapname.c_str());
        if(DEBUG){
          m_hitMap[a][b][c].draw();
        }
      }
    }
  }
//   std::cout << "SETTING stuff: " << m_nEvtsPerStep << std::endl;
  m_cache.setMaximumLength(m_nEvtsPerStep);
  m_currentCalibrationStep=-1;
  return;
}

StatusCode MDMHitMapManager::finalize()
{
	const StatusCode sc = TupleToolBase::finalize();
  return sc;
//   return StatusCode::SUCCESS;
}

void MDMHitMapManager::tidy()
{
  while(m_cache.firstKey()>0){
    m_currentCalibrationStep=m_cache.firstKey();
    fillMaps(m_hitMap,m_currentCalibrationStep);
    dumpToNtuple();
    m_cache.removeSmartIDs(m_currentCalibrationStep);
    resetMaps(m_hitMap);
  }

  for(unsigned int a=0;a<NHPDPANELS;a++){
    for(unsigned int b=0;b<NHPDCOLUMNS;b++){
      for(unsigned int c=0;c<NHPDSPERCOLUMN;c++){
        m_hitMap[a][b][c].tidy();
      }
    }
  }
  m_cache.tidy();
}

void MDMHitMapManager::addHits( LHCb::RichSmartID::Vector & smartIDs , int cStep , unsigned int nEvt ) 
{
  //add the hits from this event to a cache
  m_nEvent = nEvt;
  m_cache.addSmartIDs(cStep,smartIDs);
  m_eventsSinceEndOfPreviousPattern++;
  if(m_cache.aStepHasCompleted()){
    //info() << "Calibration Step has completed" << endmsg;
    m_currentCalibrationStep=m_cache.theCompletedStepNumber();
    fillMaps(m_hitMap,m_currentCalibrationStep);
    dumpToNtuple();
    m_cache.removeSmartIDs(m_currentCalibrationStep);
    m_cache.resetCompletedStep();
    resetMaps(m_hitMap);
  }
}

bool MDMHitMapManager::containsOnlyNoise(Rich1HitMap& map)
{
  int nHPDsWithPeaks=0;
  int nPeaksInMainHitmap=0;
  for(unsigned int a=0;a<NHPDPANELS;a++){
    for(unsigned int b=0;b<NHPDCOLUMNS;b++){
      for(unsigned int c=0;c<NHPDSPERCOLUMN;c++){
        bool seePeaks=false;
        Peaks hitMapPeaks = map[a][b][c].peaks();
        for(Peaks::iterator ipeak=hitMapPeaks.begin();ipeak!=hitMapPeaks.end();ipeak++){
          if((*ipeak).proportion()<0.1) continue;
          nPeaksInMainHitmap++;
          seePeaks=true;
        }
        if(seePeaks) nHPDsWithPeaks++;
      }
    }
  }
  if(nPeaksInMainHitmap<6 || nHPDsWithPeaks<5) return true;
  return false;
}

void MDMHitMapManager::fillMaps(Rich1HitMap& map, int step )
{
  SmartIDVector largeSmartIDVector=m_cache.smartIDsData(step);
  for(SmartIDVector::const_iterator iter=largeSmartIDVector.begin();iter!=largeSmartIDVector.end();iter++){
    LHCb::RichSmartID::Vector smallSmartIDVector = (*iter);
    for ( LHCb::RichSmartID::Vector::const_iterator iS = smallSmartIDVector.begin();iS != smallSmartIDVector.end(); ++iS ){
      int col = (*iS).pixelCol();
      int row = Rich::DAQ::NumAlicePixelsPerLHCbPixel*(*iS).pixelRow()+(*iS).pixelSubRow();
      map[(*iS).panel()][(*iS).pdCol()][(*iS).pdNumInCol()].addHit(col,row);
    }
  }
}

void MDMHitMapManager::resetMaps(Rich1HitMap& map)
{
  // reset to zero the hit counters
  for(unsigned int a=0;a<NHPDPANELS;a++){
    for(unsigned int b=0;b<NHPDCOLUMNS;b++){
      for(unsigned int c=0;c<NHPDSPERCOLUMN;c++){
        map[a][b][c].clear();
      }
    }
  }
}

void MDMHitMapManager::dumpToNtuple()
{
  std::string s_currentCalibrationStep = boost::lexical_cast<std::string>(m_currentCalibrationStep);
  if(m_currentCalibrationStep<1000) s_currentCalibrationStep="0"+s_currentCalibrationStep;
  if(m_currentCalibrationStep<100) s_currentCalibrationStep="0"+s_currentCalibrationStep;
  if(m_currentCalibrationStep<10) s_currentCalibrationStep="0"+s_currentCalibrationStep;
  int nU_HPDs=0;
  int nU_Peaks=0;
  int nD_HPDs=0;
  int nD_Peaks=0;
  std::vector<float> maximum;
  std::vector<float> Box;
  std::vector<float> HPD;
  std::vector<float> muX;
  std::vector<float> muY;
  std::vector<float> chi2;
  std::vector<float> nPEs;
  std::vector<float> rmsX;
  std::vector<float> rmsY;
  std::vector<float> meanX;
  std::vector<float> meanY;
  std::vector<float> sigmaX;
  std::vector<float> sigmaY;
  std::vector<float> rotation;
  std::vector<float> HPDColumn;
  std::vector<float> proportion;
  std::vector<float> classOfPeak;

  for(unsigned int a=0;a<NHPDPANELS;a++){
    for(unsigned int b=0;b<NHPDCOLUMNS;b++){
      for(unsigned int c=0;c<NHPDSPERCOLUMN;c++){

        m_hitMap[a][b][c].howManyEvents(m_cache.length(m_currentCalibrationStep));
        Peaks hitMapPeaks = m_hitMap[a][b][c].peaks(FIT);
        if(hitMapPeaks.size()==0) continue;

        for(Peaks::iterator ipeak=hitMapPeaks.begin();ipeak!=hitMapPeaks.end();ipeak++){
          Box.push_back(a);
          HPD.push_back(c);
          HPDColumn.push_back(b);
          maximum.push_back((*ipeak).maximum());
          muX.push_back((*ipeak).muX());
          muY.push_back((*ipeak).muY());
          rmsX.push_back((*ipeak).rmsX());
          rmsY.push_back((*ipeak).rmsY());
          rmsX.push_back((*ipeak).rmsX());
          rmsY.push_back((*ipeak).rmsY());
          chi2.push_back((*ipeak).chi2());
          nPEs.push_back((*ipeak).nHits());
          meanX.push_back((*ipeak).meanX());
          meanY.push_back((*ipeak).meanY());
          sigmaX.push_back((*ipeak).sigmaX());
          sigmaY.push_back((*ipeak).sigmaY());
          rotation.push_back((*ipeak).rotation());
          proportion.push_back((*ipeak).proportion());
          classOfPeak.push_back((*ipeak).classOfPeak());	  
        }
	if(a){
	  nD_Peaks+=hitMapPeaks.size();
	  nD_HPDs++;
	}else{
	  nU_Peaks+=hitMapPeaks.size();
	  nU_HPDs++;
	}
      }
    }
  }

  bool histStored=false;
  if(!containsOnlyNoise(m_hitMap)) histStored=true;

  Tuple ntuple = nTuple( "Map", "Result of MDM peak finding", CLID_ColumnWiseTuple );
  ntuple->column( "HistStored",histStored                                          );
  ntuple->column( "CalibStep", m_currentCalibrationStep                            );
  ntuple->column( "nEvents",   m_cache.length(m_currentCalibrationStep)            );
  ntuple->column( "eventCount",m_eventsSinceEndOfPreviousPattern                   );
  ntuple->column( "nU_Peaks",  nU_Peaks                                            );
  ntuple->column( "nU_HPDs",   nU_HPDs                                             );
  ntuple->column( "nD_Peaks",  nD_Peaks                                            );
  ntuple->column( "nD_HPDs",   nD_HPDs                                             );
  ntuple->farray( "max",       maximum.begin(),    maximum.end(),     "nPeaks",1000);
  ntuple->farray( "Box",       Box.begin(),        Box.end(),         "nPeaks",1000);
  ntuple->farray( "HPDColumn", HPDColumn.begin(),  HPDColumn.end(),   "nPeaks",1000);
  ntuple->farray( "HPD",       HPD.begin(),        HPD.end(),         "nPeaks",1000);
  ntuple->farray( "nPEs",      nPEs.begin(),       nPEs.end(),        "nPeaks",1000);
  ntuple->farray( "chi2",      chi2.begin(),       chi2.end(),        "nPeaks",1000);
  ntuple->farray( "proportion",proportion.begin(), proportion.end(),  "nPeaks",1000);
  ntuple->farray( "meanX",     meanX.begin(),      meanX.end(),       "nPeaks",1000);
  ntuple->farray( "meanY",     meanY.begin(),      meanY.end(),       "nPeaks",1000);
  ntuple->farray( "rmsX",      rmsX.begin(),       rmsX.end(),        "nPeaks",1000);
  ntuple->farray( "rmsY",      rmsY.begin(),       rmsY.end(),        "nPeaks",1000);
  ntuple->farray( "muX",       muX.begin(),        muX.end(),         "nPeaks",1000);
  ntuple->farray( "muY",       muY.begin(),        muY.end(),         "nPeaks",1000);
  ntuple->farray( "sigmaX",    sigmaX.begin(),     sigmaX.end(),      "nPeaks",1000);
  ntuple->farray( "sigmaY",    sigmaY.begin(),     sigmaY.end(),      "nPeaks",1000);
  ntuple->farray( "rotation",  rotation.begin(),   rotation.end(),    "nPeaks",1000);
  ntuple->farray( "peakClass", classOfPeak.begin(),classOfPeak.end(), "nPeaks",1000);
  ntuple->write();
	
  if(histStored && m_StoreHistos){
    for(unsigned int a=0;a<NHPDPANELS;a++){
      for(unsigned int b=0;b<NHPDCOLUMNS;b++){
	for(unsigned int c=0;c<NHPDSPERCOLUMN;c++){
	  std::string hName = "Step"+s_currentCalibrationStep+"_"+m_hitMap[a][b][c].getName();
	  for(int i = 1; i<=m_hitMap[a][b][c].nXBins();i++){
	    for(int j = 1; j<=m_hitMap[a][b][c].nYBins();j++){
	      float col = i-1;
	      float row =(j-1)/8.-0.5;
	      float content = m_hitMap[a][b][c].binContent(i,j);
	      AIDA::IHistogram2D *aida = plot2D(col,row,hName.c_str(),hName.c_str(),-0.5,31.5,-0.5,31.5,32,m_hitMap[a][b][c].nYBins(),content);
	      TH2D* hist = Gaudi::Utils::Aida2ROOT::aida2root(aida);
	      hist->SetEntries(hist->GetEntries()+content-1);
	    }
	  }
	}
      }
    }
  }

  info()<<"Calib.Step "<<m_currentCalibrationStep<<" has U("<<nU_Peaks<<" peaks on "<<nU_HPDs<<"HPDs) & D("<<nD_Peaks<<" peaks on "<<nD_HPDs<<"HPDs)"
	<<" N(since last store)="<<m_eventsSinceEndOfPreviousPattern<<" N/CS="<<m_cache.length(m_currentCalibrationStep) 
	<<(histStored?"   HISTOS STORED  ":"      noise...    ")<<"N(steps left in cache)="<<m_cache.length()-1<<endmsg;

  if(histStored) m_eventsSinceEndOfPreviousPattern=0;
}
