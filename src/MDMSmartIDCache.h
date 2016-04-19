// $Id: MDMSmartIDCache.h,v 1.1.1.1 2009/04/14 16:06:30 mjohn Exp $
#ifndef MDMSMARTIDCACHE_H 
#define MDMSMARTIDCACHE_H 1

// Include files
#include "Kernel/RichSmartID.h"

#include <vector>
#include <map>

struct lessthan{  
  bool operator()(int i1, int i2) const {
    return (i1 < i2);
  }
};

typedef std::vector<LHCb::RichSmartID::Vector> SmartIDVector;
typedef std::map< int, SmartIDVector, lessthan > SmartIDMap;

class MDMSmartIDCache{
 public:
  MDMSmartIDCache(){}
  bool aStepHasCompleted(){return (m_CompleteCalibrationStep >= 0);}
  int theCompletedStepNumber(){return m_CompleteCalibrationStep;}
  void resetCompletedStep(){m_CompleteCalibrationStep=-1;}
  unsigned int firstKey(){for(SmartIDMap::iterator it=m_map.begin();it != m_map.end();it++){if((*it).second.size()>0) return (*it).first;} return 0;}
  SmartIDVector& smartIDsData(int step){return m_map[step];}
  void setMaximumLength(int max){m_maximumLength=max;}
  void removeSmartIDs(int step){m_map.erase(step);}
  void addSmartIDs(int,LHCb::RichSmartID::Vector);
  unsigned int length(int step){return m_map[step].size();}
  unsigned int length(){return m_map.size();}
  void tidy(){m_map.clear();}
 private:
  int m_CompleteCalibrationStep;
  unsigned int m_maximumLength;
  SmartIDMap m_map;
};

inline void MDMSmartIDCache::addSmartIDs(int calibrationStep, LHCb::RichSmartID::Vector v)
{
  m_map[calibrationStep].push_back(v);
  m_CompleteCalibrationStep=-1;
  for(SmartIDMap::iterator it=m_map.begin();it != m_map.end();it++){
    //std::cout << calibrationStep << "    " << (*it).second.size() << "     " << m_maximumLength << std::endl; 
    if((*it).second.size()==m_maximumLength){
      m_CompleteCalibrationStep=(*it).first;
    }
    if((*it).second.size()>m_maximumLength){
      std::cout << "Step " << (*it).first << "contains " << (*it).second.size() << " events. Expecting " << m_maximumLength << std::endl;
    }
  }
}

#endif // MDMSMARTIDCACHE_H
