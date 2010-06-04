/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <init/Ionit_Initializer.h>

#include <exodusII/Ioex_IOFactory.h>
#include <heartbeat/Iohb_DatabaseIO.h>
#include <generated/Iogn_DatabaseIO.h>
#include <Ioss_ConcreteVariableType.h>
#include <Ioss_Initializer.h>
#include <transform/Iotr_Initializer.h>

namespace Ioss {
  namespace Init {
    Initializer::Initializer()
    {
      Ioex::IOFactory::factory();     // ExodusII
      Iohb::IOFactory::factory();   // HeartBeat
      Iogn::IOFactory::factory();  // Generated

      Ioss::StorageInitializer();
      Ioss::Initializer();
      Iotr::Initializer();
    }

    Initializer::~Initializer()
    {
      try {
	Ioss::IOFactory::clean();
	// Put code here that should run after sierra is finished
	// executing...
      } catch (...) {}
    }
  }
}
