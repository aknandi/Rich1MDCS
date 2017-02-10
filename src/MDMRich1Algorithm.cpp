// Include files 
#include <boost/lexical_cast.hpp>

// from Gaudi
#include "GaudiKernel/DeclareFactoryEntries.h" 
#include "GaudiKernel/AlgFactory.h" 
#include "GaudiUtils/Aida2ROOT.h"

// Rich
//#include "RichKernel/RichSmartIDCnv.h"
#include "RichUtils/RichSmartIDCnv.h"
#include "RichUtils/RichDecodedData.h"

// local
#include "MDMRich1Algorithm.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MDMRich1Algorithm
//
// 2009-03-25 : Malcolm John
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MDMRich1Algorithm );

typedef Rich::PoolMap< Rich::DAQ::Level1Input, Rich::DAQ::PDInfo > HPDMap;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MDMRich1Algorithm::MDMRich1Algorithm( const std::string& name, ISvcLocator* pSvcLocator)
  : Rich::Rec::TupleAlgBase ( name , pSvcLocator )
  , nEvents(0)
{
  declareProperty("DEBUG",DEBUG=false);
  declareProperty("NumberOfEventsPerStep",m_nEvtsPerStep=5000);
  declareProperty("StoreHistos",m_StoreHistos=false);
}
//=============================================================================
// Destructor
//=============================================================================
MDMRich1Algorithm::~MDMRich1Algorithm() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MDMRich1Algorithm::initialize()
{
  StatusCode sc = Rich::Rec::TupleAlgBase::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm
  debug() << "==> Initialize" << endmsg;

  acquireTool( "RichSmartIDDecoder", m_SmartIDDecoder, 0, true );
  m_mdmHitMapManager = tool<IHitMapManagerTool>("MDMHitMapManager","MDMHitManager");
  m_mdmHitMapManager->setup(m_nEvtsPerStep,DEBUG,m_StoreHistos);
	if(!m_StoreHistos) info() << "Hitmap histograms will not be stored." << endmsg;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MDMRich1Algorithm::execute()
{
  if (!exist<LHCb::ODIN>(LHCb::ODINLocation::Default) ){
    fatal()<<"No ODIN bank found"<<endmsg;
  }
  LHCb::ODIN* odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  //std::cout << odin->calibrationStep() << std::endl;
  const Rich::DAQ::L1Map & allSmartIDs = m_SmartIDDecoder->allRichSmartIDs();    
  if (allSmartIDs.size() <= 0)	return StatusCode::SUCCESS;

  std::vector<LHCb::RichSmartID>   smartIDList; 	
  // --- Loop over L1 boards
  for ( Rich::DAQ::L1Map::const_iterator iL1 = allSmartIDs.begin();iL1 != allSmartIDs.end(); ++iL1 ){
    // --- loop over ingresses for this L1 board
    for ( Rich::DAQ::IngressMap::const_iterator iIn = (*iL1).second.begin();iIn != (*iL1).second.end(); ++iIn ){
      // --- Loop over HPDs in this ingress
      for ( HPDMap::const_iterator iAllSmartIDs = (*iIn).second.pdData().begin();
            iAllSmartIDs!= (*iIn).second.pdData().end(); ++iAllSmartIDs){
        const LHCb::RichSmartID smartIDHPD = (*iAllSmartIDs).second.pdID();
        Rich::DetectorType 	   RichNum = smartIDHPD.rich ();
        Rich::Side 	            PanNum = smartIDHPD.panel ();
        unsigned int                HPDCol = smartIDHPD.pdCol();
        unsigned int                HPDRow = smartIDHPD.pdNumInCol();
        const std::string             sPan = PanNum? "LowerBox":"UpperBox";

        if(RichNum != Rich::Rich1) continue;
        
        const Rich::DAQ::PDInfo & hpdInfo = (*iAllSmartIDs).second;
        const Rich::DAQ::RichDAQHeaderV4::RichDAQHeaderPD & header = hpdInfo.header();
        plot2D(HPDCol,HPDRow,sPan+" HPD hit count",-0.5,6.5,-0.5,13.5,7,14);
        if(header.extendedFormat()){
          plot2D(HPDCol,HPDRow,sPan+" HPDs with extended Format",-0.5,6.5,-0.5,13.5,7,14);
          continue;
        }
        
        const LHCb::RichSmartID::Vector &smartIDHits = (*iAllSmartIDs).second.smartIDs();  // SmartID array of hits in HPD
        for (LHCb::RichSmartID::Vector::const_iterator id = smartIDHits.begin();id != smartIDHits.end(); id++){ 
          smartIDList.push_back(*id);
        }
      }
    }
  } 
  m_mdmHitMapManager->addHits(smartIDList,odin->calibrationStep(),nEvents++);
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MDMRich1Algorithm::finalize()
{
  debug() << "==> Finalize" << endmsg;

  m_mdmHitMapManager->tidy();
  return Rich::Rec::TupleAlgBase::finalize();  // must be called after all other actions
}
