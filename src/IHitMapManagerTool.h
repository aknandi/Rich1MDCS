#ifndef Rich1MDCS_IHitMapManagerTool_H
#define Rich1MDCS_IHitMapManagerTool_H 1

// from Gaudi
#include "GaudiKernel/IAlgTool.h"

// Kernel
#include "Kernel/RichSmartID.h"

/// Static Interface Identification
static const InterfaceID IID_IHitMapManagerTool( "IHitMapManagerTool", 1, 0 );

class IHitMapManagerTool : public virtual IAlgTool
{
 public:

  static const InterfaceID& interfaceID() { return IID_IHitMapManagerTool; }

  virtual void setup(int,bool,bool) = 0;
  virtual void addHits( LHCb::RichSmartID::Vector&, int, unsigned int ) = 0;
  virtual void tidy() = 0;

};

#endif // Rich1MDCS_IHitMapManagerTool_H
