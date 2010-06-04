/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#if !defined(__PUMAGON__) || (defined(__PUMAGON__) && defined(JANUS2))
#include <xdmf/Ioxf_Initializer.h>
#include <xdmf/Ioxf_IOFactory.h>

namespace Ioxf {
  int Initializer::useCount = 0;

  Initializer::Initializer()
  {
    if (useCount == 0) 
      Ioxf::IOFactory::factory();
    useCount++;
  }

  Initializer::~Initializer()
  {
    try {
      // Put code here that should run after sierra is finished
      // executing...
      useCount--;
      if (useCount <= 0)
	Ioxf::IOFactory::finalize();
    } catch (...) {
    }
  }
}

#endif
