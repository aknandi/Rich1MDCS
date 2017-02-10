#ifndef MDMHITMAPMANAGER_H
#define MDMHITMAPMANAGER_H 1

// base class
#include "RichKernel/RichTupleToolBase.h"

// Kernel
//#include "RichKernel/RichSmartIDSorter.h"
#include "RichUtils/RichSmartIDSorter.h"

// interface
#include "IHitMapManagerTool.h"

#include "MDMSmartIDCache.h"
#include "MDMHitMap.h"

class TH2D;

typedef boost::array<boost::array<boost::array<MDMHitMap,NHPDSPERCOLUMN>,NHPDCOLUMNS>,NHPDPANELS> Rich1HitMap;

class MDMHitMapManager : public Rich::TupleToolBase, virtual public IHitMapManagerTool {

 public:
  MDMHitMapManager( const std::string& type, const std::string& name, const IInterface* parent );
  virtual StatusCode initialize();
  virtual StatusCode finalize();
      
  void setup(int,bool,bool);
  void addHits( LHCb::RichSmartID::Vector&, int, unsigned int=0 );
  void tidy();
	
 protected:
  void dumpToNtuple();
  void fillMaps(Rich1HitMap& map,int);
  void resetMaps(Rich1HitMap& map);
  bool containsOnlyNoise(Rich1HitMap& map);

 private: 
  int m_eventsSinceEndOfPreviousPattern;
  int m_currentCalibrationStep;
  Rich1HitMap m_hitMap;
  MDMSmartIDCache m_cache;
  unsigned int m_nEvent;
  unsigned int m_pattern;
  unsigned int m_nEvtsPerStep;
  bool DEBUG;
  bool FIT;
  bool m_StoreHistos;
};

#endif // MDMHITMAPMANAGER_H
