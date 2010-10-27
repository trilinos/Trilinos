/* Copyright 2001 Sandia Corporation, Albuquerque, NM. */
#ifndef Ioi_Observer_h
#define Ioi_Observer_h

#include <Ioi_Events.h>

namespace sierra {
namespace Frio {
  class IOBase;
}
}

namespace Ioss {
namespace Interface {

  struct IOEvent {
    sierra::Frio::IOBase *		    m_ioBase;
    EventInterest   m_interest;
    EventState	    m_state;
  };
}//namespace Interface
}//namespace Ioss

#endif // Ioi_Observer_h
