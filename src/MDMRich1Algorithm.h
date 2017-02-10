// $Id: MDMRich1Algorithm.h,v 1.1.1.1 2009/04/14 16:06:30 mjohn Exp $
#ifndef MDMRICH1ALGORITHM_H 
#define MDMRICH1ALGORITHM_H 1

// from Gaudi
#include "RichKernel/RichTupleAlgBase.h"
#include "RichRecBase/RichRecTupleAlgBase.h"

#include "GaudiAlg/GaudiTupleAlg.h"

#include "RichInterfaces/IRichRawBufferToSmartIDsTool.h"

#include "MDMHitMapManager.h"

class MDMRich1Algorithm : public Rich::Rec::TupleAlgBase {
 public: 
  /// Standard constructor
  MDMRich1Algorithm( const std::string& name, ISvcLocator* pSvcLocator );
			
  virtual ~MDMRich1Algorithm( ); ///< Destructor
			
  virtual StatusCode initialize();
  virtual StatusCode execute   ();
  virtual StatusCode finalize  ();
			
 private:
  const Rich::DAQ::IRawBufferToSmartIDsTool *m_SmartIDDecoder;

  IHitMapManagerTool *m_mdmHitMapManager;
  unsigned int nEvents;
  int m_nEvtsPerStep;
  bool DEBUG;
  bool m_StoreHistos;
};

#endif // MDMRICH1ALGORITHM_H
